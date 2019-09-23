#!/usr/bin/python

import argparse

from PIL import Image


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Combines two PNG files into a single PNG by averaging the pixel RGB values from input PNGs. PNGs need to be of same dimensions ')
    parser.add_argument('senseImage', help='PNG file')
    parser.add_argument('antiSenseImage', help='PNG file')
    args = parser.parse_args()

    # open first image
    image1 = Image.open(args.senseImage)

    # open second image
    image2 = Image.open(args.antiSenseImage)

    # retrieve the image dimensions.
    width, height = image1.size
    width2, height2 = image2.size

    if [width, height] != [width2, height2]:
        print("Image dimensions do not match! The Two inputs must have equal dimensions")
        exit()
    else:
        print(image1.size)
        print(image2.size)
        # Create a new image object.
        merged = Image.new('RGB', image1.size)

        total = 0
        for i in range(0, width):
            for j in range(0, height):
                ima1 = list(image1.getpixel((i, j)))
                ima2 = list(image2.getpixel((i, j)))
                if ima1 == ima2:
                    r, g, b, a = ima1
                elif [ima1[0], ima1[1], ima1[2]] == [0, 0, 0] and [ima2[0], ima2[1], ima2[2]] != [0, 0, 0]:
                    r, g, b, a = ima2
                elif [ima1[0], ima1[1], ima1[2]] != [0, 0, 0] and [ima2[0], ima2[1], ima2[2]] == [0, 0, 0]:
                    r, g, b, a = ima1
                elif [ima1[0], ima1[1], ima1[2]] != [0, 0, 0] and ima2 == [255, 255, 255, 255]:
                    r, g, b, a = ima1
                elif [ima2[0], ima2[1], ima2[2]] != [0, 0, 0] and ima1 == [255, 255, 255, 255]:
                    r, g, b, a = ima2
                else:
                    # print ima1,ima2
                    r = (ima1[0] + ima2[0]) // 2
                    g = (ima1[1] + ima2[1]) // 2
                    b = (ima1[2] + ima2[2]) // 2
                    a = 255
                    # print [r,g,b,a]

                merged.putpixel((i, j), (r, g, b, a))
        merged.save("merge.png")
