import argparse
import json
import os
import sys
import tarfile
from contextlib import ExitStack
from tempfile import TemporaryDirectory

from libtiff import TIFF
from omero.gateway import BlitzGateway  # noqa
from omero.constants.namespaces import NSBULKANNOTATIONS  # noqa


def warn(message, image_identifier, warn_skip=False):
    message = message.rstrip()
    if warn_skip:
        if message[-1] in ['.', '!', '?']:
            skip_msg = ' Skipping download!'
        else:
            skip_msg = '. Skipping download!'
    else:
        skip_msg = ''
    print(
        'ImageSpecWarning for {0}: {1}{2}'
        .format(
            image_identifier,
            message,
            skip_msg
        ),
        file=sys.stderr
    )


def find_channel_index(image, channel_name):
    channel_name = channel_name.lower()
    for n, channel in enumerate(image.getChannelLabels()):
        if channel_name == channel.lower():
            return n
    # Check map annotation for information (this is necessary for some images)
    for ann in image.listAnnotations(NSBULKANNOTATIONS):
        pairs = ann.getValue()
        for p in pairs:
            if p[0] == "Channels":
                channels = p[1].replace(" ", "").split(";")
                for n, c in enumerate(channels):
                    for value in c.split(':'):
                        if channel_name == value.lower():
                            return n
    return -1


def get_clipping_region(image, x, y, w, h):
    # If the (x, y) coordinate falls outside the image boundaries, we
    # cannot just shift it because that would render the meaning of
    # w and h undefined (should width and height be decreased or the whole
    # region be shifted to keep them fixed?).
    # It may be better to abort in this situation.
    if x < 0 or y < 0:
        raise ValueError(
            'Too small upper left coordinate ({0}, {1}) for clipping region.'
            .format(x, y)
        )
    size_x = image.getSizeX()
    size_y = image.getSizeY()
    if x >= size_x or y >= size_y:
        raise ValueError(
            'Upper left coordinate ({0}, {1}) of clipping region lies '
            'outside of image.'
            .format(x, y)
        )
    # adjust width and height to the image dimensions
    if w <= 0 or x + w > size_x:
        w = size_x - x
    if h <= 0 or y + h > size_y:
        h = size_y - y
    return [x, y, w, h]


def confine_plane(image, z):
    if z < 0:
        z = 0
    else:
        max_z = image.getSizeZ() - 1
        if z > max_z:
            z = max_z
    return z


def confine_frame(image, t):
    if t < 0:
        t = 0
    else:
        max_t = image.getSizeT() - 1
        if t > max_t:
            t = max_t
    return t


def get_image_array(image, tile, z, c, t):
    pixels = image.getPrimaryPixels()
    try:
        selection = pixels.getTile(theZ=z, theT=t, theC=c, tile=tile)
    except Exception:
        warning = '{0} (ID: {1})'.format(image.getName(),
                                         image.getId())
        warn('Could not download the requested region', warning)
        return

    return selection


def download_image_data(
    image_ids,
    channel=None, z_stack=0, frame=0,
    coord=(0, 0), width=0, height=0, region_spec='rectangle',
    skip_failed=False, download_tar=False, omero_host='idr.openmicroscopy.org', omero_secured=False, config_file=None
):

    if config_file is None:  # IDR connection
        omero_username = 'public'
        omero_password = 'public'
    else:  # other omero instance
        with open(config_file) as f:
            cfg = json.load(f)
            omero_username = cfg['username']
            omero_password = cfg['password']

            if omero_username == "" or omero_password == "":
                omero_username = 'public'
                omero_password = 'public'

    # basic argument sanity checks and adjustments
    prefix = 'image-'
    # normalize image ids by stripping off prefix if it exists
    image_ids = [
        iid[len(prefix):] if iid[:len(prefix)] == prefix else iid
        for iid in image_ids
    ]

    if region_spec not in ['rectangle', 'center']:
        raise ValueError(
            'Got unknown value "{0}" as region_spec argument'
            .format(region_spec)
        )
    with ExitStack() as exit_stack:
        conn = exit_stack.enter_context(
            BlitzGateway(
                omero_username, omero_password,
                host=omero_host,
                secure=omero_secured
            )
        )
        # exit_stack.callback(conn.connect().close)
        if download_tar:
            # create an archive file to write images to
            archive = exit_stack.enter_context(
                tarfile.open('images.tar', mode='w')
            )
            tempdir = exit_stack.enter_context(
                TemporaryDirectory()
            )

        for image_id in image_ids:
            image_warning_id = 'Image-ID: {0}'.format(image_id)
            try:
                image_id = int(image_id)
            except ValueError:
                image = None
            else:
                try:
                    image = conn.getObject("Image", image_id)
                except Exception as e:
                    # respect skip_failed on unexpected errors
                    if skip_failed:
                        warn(str(e), image_warning_id, warn_skip=True)
                        continue
                    else:
                        raise

            if image is None:
                if skip_failed:
                    warn(
                        'Unable to find an image with this ID in the '
                        'database.',
                        image_warning_id,
                        warn_skip=True
                    )
                    continue
                raise ValueError(
                    '{0}: Unable to find an image with this ID in the '
                    'database. Aborting!'
                    .format(image_warning_id)
                )

            try:
                # try to extract image properties
                # if anything goes wrong here skip the image
                # or abort.
                image_name = os.path.splitext(image.getName())[0]
                image_warning_id = '{0} (ID: {1})'.format(
                    image_name, image_id
                )

                if region_spec == 'rectangle':
                    tile = get_clipping_region(image, *coord, width, height)
                elif region_spec == 'center':
                    tile = get_clipping_region(
                        image,
                        *_center_to_ul(*coord, width, height)
                    )

                ori_z, z_stack = z_stack, confine_plane(image, z_stack)
                ori_frame, frame = frame, confine_frame(image, frame)
                num_channels = image.getSizeC()
                if channel is None:
                    channel_index = 0
                else:
                    channel_index = find_channel_index(image, channel)
            except Exception as e:
                # respect skip_failed on unexpected errors
                if skip_failed:
                    warn(str(e), image_warning_id, warn_skip=True)
                    continue
                else:
                    raise

            # region sanity checks and warnings
            if tile[2] < width or tile[3] < height:
                # The downloaded image region will have smaller dimensions
                # than the specified width x height.
                warn(
                    'Downloaded image dimensions ({0} x {1}) will be smaller '
                    'than the specified width and height ({2} x {3}).'
                    .format(tile[2], tile[3], width, height),
                    image_warning_id
                )

            # z-stack sanity checks and warnings
            if z_stack != ori_z:
                warn(
                    'Specified image plane ({0}) is out of bounds. Using {1} '
                    'instead.'
                    .format(ori_z, z_stack),
                    image_warning_id
                )

            # frame sanity checks and warnings
            if frame != ori_frame:
                warn(
                    'Specified image frame ({0}) is out of bounds. Using '
                    'frame {1} instead.'
                    .format(ori_frame, frame),
                    image_warning_id
                )

            # channel index sanity checks and warnings
            if channel is None:
                if num_channels > 1:
                    warn(
                        'No specific channel selected for multi-channel '
                        'image. Using first of {0} channels.'
                        .format(num_channels),
                        image_warning_id
                    )
            else:
                if channel_index == -1 or channel_index >= num_channels:
                    if skip_failed:
                        warn(
                            str(channel)
                            + ' is not a known channel name for this image.',
                            image_warning_id,
                            warn_skip=True
                        )
                        continue
                    else:
                        raise ValueError(
                            '"{0}" is not a known channel name for image {1}. '
                            'Aborting!'
                            .format(channel, image_warning_id)
                        )

            # download and save the region as TIFF
            fname = '__'.join(
                [image_name, str(image_id)] + [str(x) for x in tile]
            )
            try:
                if fname[-5:] != '.tiff':
                    fname += '.tiff'

                fname = fname.replace(' ', '_')

                im_array = get_image_array(image, tile, z_stack, channel_index, frame)

                if download_tar:
                    fname = os.path.join(tempdir, fname)
                try:
                    tiff = TIFF.open(fname, mode='w')
                    tiff.write_image(im_array)
                finally:
                    tiff.close()
                # move image into tarball
                if download_tar:
                    archive.add(fname, os.path.basename(fname))
                    os.remove(fname)
            except Exception as e:
                if skip_failed:
                    # respect skip_failed on unexpected errors
                    warn(str(e), image_warning_id, warn_skip=True)
                    continue
                else:
                    raise


def _center_to_ul(center_x, center_y, width, height):
    if width > 0:
        ext_x = (width - 1) // 2
        ul_x = max([center_x - ext_x, 0])
        width = center_x + ext_x + 1 - ul_x
    else:
        ul_x = 0
    if height > 0:
        ext_y = (height - 1) // 2
        ul_y = max([center_y - ext_y, 0])
        height = center_y + ext_y + 1 - ul_y
    else:
        ul_y = 0
    return ul_x, ul_y, width, height


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        'image_ids', nargs='*', default=[],
        help='one or more IDR image ids for which to retrieve data (default: '
             'read ids from stdin).'
    )
    p.add_argument(
        '-c', '--channel',
        help='name of the channel to retrieve data for '
             '(note: the first channel of each image will be downloaded if '
             'left unspecified)'
    )
    region = p.add_mutually_exclusive_group()
    region.add_argument(
        '--rectangle', nargs=4, type=int, default=argparse.SUPPRESS,
        help='specify a clipping region for the image as x y width height, '
             'where x and y give the upper left coordinate of the rectangle '
             'to clip to. Set width and height to 0 to extend the rectangle '
             'to the actual size of the image.'
    )
    region.add_argument(
        '--center', nargs=4, type=int, default=argparse.SUPPRESS,
        help='specify a clipping region for the image as x y width height, '
             'where x and y define the center of a width x height rectangle. '
             'Set either width or height to 0 to extend the region to the '
             'actual size of the image along the x- or y-axis.\n'
             'Note: Even values for width and height will be rounded down to '
             'the nearest odd number.'
    )
    p.add_argument(
        '-f', '--frame', type=int, default=0
    )
    p.add_argument(
        '-z', '--z-stack', type=int, default=0
    )
    p.add_argument(
        '--skip-failed', action='store_true'
    )
    p.add_argument(
        '--download-tar', action='store_true'
    )
    p.add_argument(
        '-oh', '--omero-host', type=str, default="idr.openmicroscopy.org"
    )
    p.add_argument(
        '--omero-secured', action='store_true', default=True
    )
    p.add_argument(
        '-cf', '--config-file', dest='config_file', default=None
    )
    args = p.parse_args()
    if not args.image_ids:
        args.image_ids = sys.stdin.read().split()
    if 'center' in args:
        args.coord, args.width, args.height = (
            args.center[:2], args.center[2], args.center[3]
        )
        args.region_spec = 'center'
        del args.center
    elif 'rectangle' in args:
        args.coord, args.width, args.height = (
            args.rectangle[:2], args.rectangle[2], args.rectangle[3]
        )
        args.region_spec = 'rectangle'
        del args.rectangle
    download_image_data(**vars(args))
