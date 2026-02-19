import argparse
import json
import os
import sys
import tarfile
from contextlib import ExitStack
from itertools import product
from tempfile import TemporaryDirectory

import numpy
from omero.cli import cli_login
from omero.constants.namespaces import NSBULKANNOTATIONS
from omero.gateway import _ImageWrapper, BlitzGateway
from tifffile import imwrite


def warn(message: str, image_identifier: str, warn_skip: bool = False) -> None:
    """Print an error `message` to stderr and
    - prefix with the `image_identifier`
    - suffix with 'Skipping download!' if `warn_skip` is True

    Args:
        message (string): Message to print to stderr
        image_identifier (string): Image identifier
        warn_skip (bool, optional): Whether 'skipping download' should be suffix to the message. Defaults to False.
    """
    message = message.rstrip()
    if warn_skip:
        if message[-1] in [".", "!", "?"]:
            skip_msg = " Skipping download!"
        else:
            skip_msg = ". Skipping download!"
    else:
        skip_msg = ""
    print(
        "ImageSpecWarning for {0}: {1}{2}".format(image_identifier, message, skip_msg),
        file=sys.stderr,
    )


def find_channel_index(image: _ImageWrapper, channel_name: str) -> int:
    """Identify the channel index from the `image` and the `channel_name`

    Args:
        image (_ImageWrapper): image wrapper on which the channel should be identified
        channel_name (string): name of the channel to look for

    Returns:
        int: Index of the channel or -1 if was not found.
    """
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
                    for value in c.split(":"):
                        if channel_name == value.lower():
                            return n
    return -1


def get_clipping_region(
    image: _ImageWrapper, x: int, y: int, w: int, h: int
) -> list[int]:
    """Check `x`, `y` and adjust `w`, `h` to image size to be able to crop the `image` with these coordinates

    Args:
        image (_ImageWrapper): image wrapper on which region want to be cropped
        x (int): left x coordinate
        y (int): top y coordinate
        w (int): width
        h (int): height

    Raises:
        ValueError: if the x or y coordinates are negative.
        ValueError: if the x or y coordinates are larger than the width or height of the image.

    Returns:
        list[int]: new [x, y, width, height] adjusted to the image
    """
    # If the (x, y) coordinate falls outside the image boundaries, we
    # cannot just shift it because that would render the meaning of
    # w and h undefined (should width and height be decreased or the whole
    # region be shifted to keep them fixed?).
    # It may be better to abort in this situation.
    if x < 0 or y < 0:
        raise ValueError(
            "Too small upper left coordinate ({0}, {1}) for clipping region.".format(
                x, y
            )
        )
    size_x = image.getSizeX()
    size_y = image.getSizeY()
    if x >= size_x or y >= size_y:
        raise ValueError(
            "Upper left coordinate ({0}, {1}) of clipping region lies "
            "outside of image.".format(x, y)
        )
    # adjust width and height to the image dimensions
    if w <= 0 or x + w > size_x:
        w = size_x - x
    if h <= 0 or y + h > size_y:
        h = size_y - y
    return [x, y, w, h]


def confine_plane(image: _ImageWrapper, z: int) -> int:
    """Adjust/Confine `z` to be among the possible z for the `image`

    Args:
        image (_ImageWrapper): image wrapper for which the z is adjusted
        z (int): plane index that need to be confined

    Returns:
        int: confined z
    """
    if z < 0:
        z = 0
    else:
        max_z = image.getSizeZ() - 1
        if z > max_z:
            z = max_z
    return z


def confine_frame(image: _ImageWrapper, t: int) -> int:
    """Adjust/Confine `t` to be among the possible t for the `image`

    Args:
        image (_ImageWrapper): image wrapper for which the t is adjusted
        t (int): frame index that need to be confined

    Returns:
        int: confined t
    """
    if t < 0:
        t = 0
    else:
        max_t = image.getSizeT() - 1
        if t > max_t:
            t = max_t
    return t


def get_image_array(
    image: _ImageWrapper, tile: list[int], z: int, c: int, t: int
) -> numpy.ndarray:
    """Get a 2D numpy array from an `image` wrapper for a given `tile`, `z`, `c`, `t`

    Args:
        image (_ImageWrapper): image wrapper from which values are taken
        tile (list[int]): [x, y, width, height] where x,y is the top left coordinate of the region to crop
        z (int): plane index
        c (int): channel index
        t (int): frame index

    Returns:
        numpy.ndarray: image values of the selected area (2 dimensions)
    """
    pixels = image.getPrimaryPixels()
    try:
        selection = pixels.getTile(theZ=z, theT=t, theC=c, tile=tile)
    except Exception:
        warning = "{0} (ID: {1})".format(image.getName(), image.getId())
        warn("Could not download the requested region", warning)
        return

    return selection


def get_full_image_array(image: _ImageWrapper) -> numpy.ndarray:
    """Get a 5D numpy array with all values from an `image` wrapper

    Args:
        image (_ImageWrapper): image wrapper from which values are taken

    Returns:
        numpy.ndarray: image values in the TZCYX order (5 dimensions)
    """
    # The goal is to get the image in TZCYX order
    pixels = image.getPrimaryPixels()
    # Get the final tzclist in the order that will make the numpy reshape work
    tzclist = list(
        product(
            range(image.getSizeT()), range(image.getSizeZ()), range(image.getSizeC())
        )
    )
    # As getPlanes requires the indices in the zct order
    # We keep the final order but switch indices
    zctlist = [(z, c, t) for (t, z, c) in tzclist]
    try:
        all_planes = numpy.array(list(pixels.getPlanes(zctlist)))
        all_planes_reshaped = all_planes.reshape(
            image.getSizeT(),
            image.getSizeZ(),
            image.getSizeC(),
            all_planes.shape[-2],
            all_planes.shape[-1],
        )
    except Exception as e:
        warning = "{0} (ID: {1})".format(image.getName(), image.getId())
        warn(f"Could not download the full image \n {e.msg}", warning)
        return

    return all_planes_reshaped


def download_image_data(
    image_ids_or_dataset_id: str,
    dataset: bool = False,
    download_original: bool = False,
    download_full: bool = False,
    channel: str = None,
    z_stack: int = 0,
    frame: int = 0,
    coord: tuple[int, int] = (0, 0),
    width: int = 0,
    height: int = 0,
    region_spec: str = "rectangle",
    skip_failed: bool = False,
    download_tar: bool = False,
    omero_host: str = "idr.openmicroscopy.org",
    omero_secured: bool = False,
    config_file: str = None,
) -> None:
    """Download the image data of
      either a list of image ids or all images from a dataset.
      The image data can be:
       - a 2D cropped region or
       - a hyperstack written in a tiff file
       - the original image uploaded in omero
      Optionally, the final file can be in a tar

    Args:
        image_ids_or_dataset_id (list[str]): Can be either a list with a single id (int) of a dataset or a list with images ids (int) or images ids prefixed by 'image-'
        dataset (bool, optional): Whether the image_ids_or_dataset_id is a dataset id and all images from this dataset should be retrieved (true) or image_ids_or_dataset_id are individual image ids (false). Defaults to False.
        download_original (bool, optional): Whether the original file uploded to omero should be downloaded (ignored if `download_full` is set to True). Defaults to False.
        download_full (bool, optional): Whether the full image (hyperstack) on omero should be written to TIFF. Defaults to False.
        channel (string, optional): Channel name (ignored if `download_full` or `download_original` is set to True). Defaults to None.
        z_stack (int, optional): Z stack (plane) index (ignored if `download_full` or `download_original` is set to True). Defaults to 0.
        frame (int, optional): T frame index (ignored if `download_full` or `download_original` is set to True). Defaults to 0.
        coord (tuple[int, int], optional): Coordinates of the top left or center of the region to crop (ignored if `download_full` or `download_original` is set to True). Defaults to (0, 0).
        width (int, optional): Width of the region to crop (ignored if `download_full` or `download_original` is set to True). Defaults to 0.
        height (int, optional): Height of the region to crop (ignored if `download_full` or `download_original` is set to True). Defaults to 0.
        region_spec (str, optional): How the region is specified ('rectangle' = coord is top left or 'center' = the region is center, ignored if `download_full` or `download_original` is set to True). Defaults to "rectangle".
        skip_failed (bool, optional): Do not stop the downloads if one fails. Defaults to False.
        download_tar (bool, optional): Put all downloaded images into a tar file. Defaults to False.
        omero_host (str, optional): omero host url. Defaults to "idr.openmicroscopy.org".
        omero_secured (bool, optional): Whether the omero connects with secure connection. Defaults to False.
        config_file (string, optional): File path with config file with credentials to connect to OMERO. Defaults to None.

    Raises:
        ValueError: If the region_spec is not 'rectangle' nor 'center' and a cropped region is wanted.
        ValueError: If there is no dataset with this number in OMERO
        ValueError: If there is no image with this number in OMERO
        Exception: If the command to download the original image fails
        ValueError: If the channel name could not be identified
    """

    if config_file is None:  # IDR connection
        omero_username = "public"
        omero_password = "public"
    else:  # other omero instance
        with open(config_file) as f:
            cfg = json.load(f)
            omero_username = cfg["username"]
            omero_password = cfg["password"]

            if omero_username == "" or omero_password == "":
                omero_username = "public"
                omero_password = "public"

    if (
        not download_original
        and not download_full
        and region_spec not in ["rectangle", "center"]
    ):
        raise ValueError(
            'Got unknown value "{0}" as region_spec argument'.format(region_spec)
        )
    with ExitStack() as exit_stack:
        conn = exit_stack.enter_context(
            BlitzGateway(
                omero_username, omero_password, host=omero_host, secure=omero_secured
            )
        )
        if download_tar:
            # create an archive file to write images to
            archive = exit_stack.enter_context(tarfile.open("images.tar", mode="w"))
            tempdir = exit_stack.enter_context(TemporaryDirectory())

        if dataset:
            dataset_warning_id = "Dataset-ID: {0}".format(image_ids_or_dataset_id[0])
            try:
                dataset_id = int(image_ids_or_dataset_id[0])
            except ValueError:
                image_ids = None
            else:
                try:
                    dataset = conn.getObject("Dataset", dataset_id)
                except Exception as e:
                    # respect skip_failed on unexpected errors
                    if skip_failed:
                        warn(str(e), dataset_warning_id, warn_skip=True)
                    else:
                        raise
                else:
                    image_ids = [image.id for image in dataset.listChildren()]

            if image_ids is None:
                if skip_failed:
                    warn(
                        "Unable to find a dataset with this ID in the " "database.",
                        dataset_warning_id,
                        warn_skip=True,
                    )
                else:
                    raise ValueError(
                        "{0}: Unable to find a dataset with this ID in the "
                        "database. Aborting!".format(dataset_warning_id)
                    )

        else:
            # basic argument sanity checks and adjustments
            prefix = "image-"
            # normalize image ids by stripping off prefix if it exists
            image_ids = [
                iid[len(prefix):] if iid[:len(prefix)] == prefix else iid
                for iid in image_ids_or_dataset_id
            ]
        for image_id in image_ids:
            image_warning_id = "Image-ID: {0}".format(image_id)
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
                        "Unable to find an image with this ID in the database.",
                        image_warning_id,
                        warn_skip=True,
                    )
                    continue
                raise ValueError(
                    "{0}: Unable to find an image with this ID in the "
                    "database. Aborting!".format(image_warning_id)
                )
            try:
                # try to extract image name
                # if anything goes wrong here skip the image
                # or abort.
                image_name = os.path.splitext(image.getName())[0]
                image_warning_id = "{0} (ID: {1})".format(image_name, image_id)
            except Exception as e:
                # respect skip_failed on unexpected errors
                if skip_failed:
                    warn(str(e), image_warning_id, warn_skip=True)
                    continue
                else:
                    raise
            if download_full:
                fname = (
                    "__".join([image_name.replace(" ", "_"), str(image_id), "full"])
                    + ".tiff"
                )
                # download and save the region as TIFF
                try:
                    im_array = get_full_image_array(image)

                    if download_tar:
                        fname = os.path.join(tempdir, fname)

                    imwrite(fname, im_array, imagej=True)
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
            elif download_original:
                try:
                    # try to extract image properties
                    # if anything goes wrong here skip the image
                    # or abort.
                    original_image_name = image.getFileset().listFiles()[0].getName()
                    fname = (
                        image_name
                        + "__"
                        + str(image_id)
                        + os.path.splitext(original_image_name)[1]
                    )
                    fname = fname.replace(" ", "_")
                    fname = fname.replace("/", "_")
                    download_directory = "./"
                    if download_tar:
                        download_directory = tempdir
                    with cli_login(
                        "-u", omero_username, "-s", omero_host, "-w", omero_password
                    ) as cli:
                        cli.invoke(
                            ["download", f"Image:{image_id}", download_directory]
                        )
                        if cli.rv != 0:
                            raise Exception("Download failed.")
                    # This will download to download_directory/original_image_name
                    os.rename(
                        os.path.join(download_directory, original_image_name),
                        os.path.join(download_directory, fname),
                    )
                    # move image into tarball
                    if download_tar:
                        archive.add(
                            os.path.join(download_directory, fname),
                            os.path.basename(fname),
                        )
                        os.remove(os.path.join(download_directory, fname))
                except Exception as e:
                    # respect skip_failed on unexpected errors
                    if skip_failed:
                        warn(str(e), image_warning_id, warn_skip=True)
                        continue
                    else:
                        raise
            else:
                try:
                    # try to extract image properties
                    # if anything goes wrong here skip the image
                    # or abort.
                    if region_spec == "rectangle":
                        tile = get_clipping_region(image, *coord, width, height)
                    elif region_spec == "center":
                        tile = get_clipping_region(
                            image, *_center_to_ul(*coord, width, height)
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
                        "Downloaded image dimensions ({0} x {1}) will be smaller "
                        "than the specified width and height ({2} x {3}).".format(
                            tile[2], tile[3], width, height
                        ),
                        image_warning_id,
                    )

                # z-stack sanity checks and warnings
                if z_stack != ori_z:
                    warn(
                        "Specified image plane ({0}) is out of bounds. Using {1} "
                        "instead.".format(ori_z, z_stack),
                        image_warning_id,
                    )

                # frame sanity checks and warnings
                if frame != ori_frame:
                    warn(
                        "Specified image frame ({0}) is out of bounds. Using "
                        "frame {1} instead.".format(ori_frame, frame),
                        image_warning_id,
                    )

                # channel index sanity checks and warnings
                if channel is None:
                    if num_channels > 1:
                        warn(
                            "No specific channel selected for multi-channel "
                            "image. Using first of {0} channels.".format(num_channels),
                            image_warning_id,
                        )
                else:
                    if channel_index == -1 or channel_index >= num_channels:
                        if skip_failed:
                            warn(
                                str(channel)
                                + " is not a known channel name for this image.",
                                image_warning_id,
                                warn_skip=True,
                            )
                            continue
                        else:
                            raise ValueError(
                                '"{0}" is not a known channel name for image {1}. '
                                "Aborting!".format(channel, image_warning_id)
                            )

                fname = "__".join([image_name, str(image_id)] + [str(x) for x in tile])
                fname += ".tiff"
                fname = fname.replace(" ", "_")
                # download and save the region as TIFF
                try:
                    im_array = get_image_array(
                        image, tile, z_stack, channel_index, frame
                    )

                    if download_tar:
                        fname = os.path.join(tempdir, fname)
                    imwrite(fname, im_array)
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


def _center_to_ul(center_x: int, center_y: int, width: int, height: int) -> list[int]:
    """Convert the center coordinates (`center_x`, `center_y`), `width`, `height` to upper left coordinates, width, height

    Args:
        center_x (int): x coordinate of center
        center_y (int): y coordinate of center
        width (int): width
        height (int): height

    Returns:
        list[int]: [x, y, width, height] where x,y are the upper left coordinates
    """
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
        "image_ids_or_dataset_id",
        nargs="*",
        default=[],
        help="one or more IDR image ids or a single dataset id"
        "for which to retrieve data (default: "
        "read ids from stdin).",
    )
    region = p.add_mutually_exclusive_group()
    region.add_argument(
        "--rectangle",
        nargs=4,
        type=int,
        default=argparse.SUPPRESS,
        help="specify a clipping region for the image as x y width height, "
        "where x and y give the upper left coordinate of the rectangle "
        "to clip to. Set width and height to 0 to extend the rectangle "
        "to the actual size of the image.",
    )
    region.add_argument(
        "--center",
        nargs=4,
        type=int,
        default=argparse.SUPPRESS,
        help="specify a clipping region for the image as x y width height, "
        "where x and y define the center of a width x height rectangle. "
        "Set either width or height to 0 to extend the region to the "
        "actual size of the image along the x- or y-axis.\n"
        "Note: Even values for width and height will be rounded down to "
        "the nearest odd number.",
    )
    region.add_argument(
        "--download-original",
        dest="download_original",
        action="store_true",
        help="download the original file uploaded to omero",
    )
    region.add_argument(
        "--download-full",
        dest="download_full",
        action="store_true",
        help="download the full image on omero",
    )
    p.add_argument(
        "-c",
        "--channel",
        help="name of the channel to retrieve data for "
        "(note: the first channel of each image will be downloaded if "
        "left unspecified), ignored with `--download-original` and "
        "`--download-full`",
    )
    p.add_argument(
        "-f",
        "--frame",
        type=int,
        default=0,
        help="index of the frame to retrive data for (first frame is 0),"
        " ignored with `--download-original` and `--download-full`",
    )
    p.add_argument(
        "-z",
        "--z-stack",
        type=int,
        default=0,
        help="index of the slice to retrive data for (first slice is 0),"
        " ignored with `--download-original` and `--download-full`",
    )
    p.add_argument("--skip-failed", action="store_true")
    p.add_argument("--download-tar", action="store_true")
    p.add_argument("-oh", "--omero-host", type=str, default="idr.openmicroscopy.org")
    p.add_argument("--omero-secured", action="store_true", default=True)
    p.add_argument("-cf", "--config-file", dest="config_file", default=None)
    p.add_argument("--dataset", action="store_true")
    args = p.parse_args()
    if not args.image_ids_or_dataset_id:
        args.image_ids_or_dataset_id = sys.stdin.read().split()
    if args.dataset and len(args.image_ids_or_dataset_id) > 1:
        warn("Multiple dataset ids provided. Only the first one will be used.")
    if "center" in args:
        args.coord, args.width, args.height = (
            args.center[:2],
            args.center[2],
            args.center[3],
        )
        args.region_spec = "center"
        del args.center
    elif "rectangle" in args:
        args.coord, args.width, args.height = (
            args.rectangle[:2],
            args.rectangle[2],
            args.rectangle[3],
        )
        args.region_spec = "rectangle"
        del args.rectangle
    download_image_data(**vars(args))
