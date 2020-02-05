import argparse
import os
import sys

from matplotlib import pyplot as plt
from omero.gateway import BlitzGateway
from omero.constants.namespaces import NSBULKANNOTATIONS


def find_channel_index(image, channel_name):
    index = -1
    channel_name = channel_name.lower()
    for n, channel in enumerate(image.getChannels()):
        if channel_name == channel.getLabel().lower():
            index = n
            break
    # Check map annotation for information (this is necessary for some images)
    if index == -1:
        for ann in image.listAnnotations(NSBULKANNOTATIONS):
            pairs = ann.getValue()
            for p in pairs:
                if p[0] == "Channels":
                    channels = p[1].replace(" ", "").split(";")
                    for n, c in enumerate(channels):
                        for value in c.split(':'):
                            if channel_name == value.lower():
                                index = n
    return index


def get_valid_region(image, tile):
    x, y, w, h = tile
    size_x = image.getSizeX()
    size_y = image.getSizeY()
    if x < 0 or x >= size_x:
        return None
    if y < 0 or y >= size_y:
        return None
    if w < 0 or w > size_x:
        return None
    if h < 0 or h > size_y:
        return None
    if x + w > size_x:
        return None
    if y + h > size_y:
        return None
    return [x, y, w, h]


def download_plane_as_tiff(image, tile, z, c, t):
    if z < 0 or z >= image.getSizeZ():
        z = 0
    if t < 0 or t >= image.getSizeT():
        t = 0
    pixels = image.getPrimaryPixels()
    region = get_valid_region(tile)
    selection = pixels.getTile(theZ=z, theT=t, theC=c, tile=region)

    # save the region as TIFF
    filename, file_extension = os.path.splitext(image.getName())
    # Name to be adjusted
    name = filename + "_" + str(image.getId()) + '_'.join([str(x) for x in tile]) + ".tiff"
    plt.imsave(name, selection)


def download_image_data(
    image_ids,
    channel=None, z_stack=0, frame=0,
    ul_coord=None, width=None, height=None
):

    # connect to idr
    conn = BlitzGateway('public', 'public',
                        host='idr.openmicroscopy.org',
                        secure=True)
    conn.connect()

    try:
        prefix = 'image-'
        for image_id in image_ids:
            if image_id[:len(prefix)] == prefix:
                image_id = image_id[len(prefix):]
            image_id = int(image_id)
            image = conn.getObject("Image", image_id)

            if ul_coord is None and width is None and height is None:
                tile = None
            else:
                tile = list(ul_coord) + [width, height]
                # TODO check that it is a valid region
                assert all(isinstance(d, int) for d in tile)

            # Get the channel index. If the index is not valid, skip the image
            if channel is None:
                channel_index = 0
            else:
                channel_index = find_channel_index(image, channel)
                if channel_index == -1:
                    raise ValueError(
                        '"{0}" is not a known channel name for image {1}'
                        .format(channel, image.getName())
                    )
            download_plane_as_tiff(image, tile, z_stack, channel_index, frame)
    finally:
        # Close the connection
        conn.close()


def center_to_ul(c_coord, width, height):
    # TO DO:
    # implement this
    return c_coord, width, height


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
        '--rectangle', nargs=4, type=int, default=argparse.SUPPRESS
    )
    region.add_argument(
        '--center', nargs=4, type=int, default=argparse.SUPPRESS
    )
    p.add_argument(
        '-f', '--frame', type=int, default=0
    )
    p.add_argument(
        '-z', '--z-stack', type=int, default=0
    )

    args = p.parse_args()
    if not args.image_ids:
        args.image_ids = sys.stdin.read().split()
    if 'center' in args:
        ul_coord, width, height = center_to_ul(
            args.center[:2], args.center[2], args.center[3]
        )
        del args.center
    elif 'rectangle' in args:
        ul_coord, width, height = (
            args.rectangle[:2], args.rectangle[2], args.rectangle[3]
        )
        del args.rectangle
    download_image_data(
        ul_coord=ul_coord, width=width, height=height, **vars(args)
    )
