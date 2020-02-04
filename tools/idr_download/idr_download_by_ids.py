import argparse
import os
import sys

from matplotlib import pyplot as plt
from omero.gateway import BlitzGateway


def download_image_data(
    image_ids,
    channel, z_stack=0, frame=0,
    ul_coord=None, width=None, height=None
):

    # connect to idr
    conn = BlitzGateway('public', 'public',
                        host='ws://idr.openmicroscopy.org/omero-ws',
                        secure=True)
    conn.connect()

    prefix = 'image-'
    for image_id in image_ids:
        # get the image
        #image_id = 1884807
        if image_id[:len(prefix)] == prefix:
            image_id = image_id[len(prefix):]
        image_id = int(image_id)
        image = conn.getObject("Image", image_id)
        image_name = image.getName()
        print(image_name)

        # example loading one single plane and saving it 
        pixels = image.getPrimaryPixels()

        if ul_coord is None and width is None and height is None:
            tile = None
        else:
            tile = list(ul_coord) + [width, height]
            assert all(isinstance(d, int) for d in tile)
        # TO DO:
        # Get the channel information
        # Check if it is the channel metadata if not in the map annotation

        selection = pixels.getTile(theZ=z_stack, theT=frame, theC=0, tile=tile)

        # save the crop image as TIFF
        filename, file_extension = os.path.splitext(image_name)
        name = filename + "_" + str(image.getId()) + '_'.join(str(x) for x in tile) + ".tiff"
        plt.imsave(name, selection)

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
        '-c', '--channel', required=True,
        help='name of the channel to retrieve data for'
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
