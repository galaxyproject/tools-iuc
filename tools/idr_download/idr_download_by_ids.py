# import libraries
import argparse
import csv
import os
import subprocess
import sys


import numpy as np
from pandas import DataFrame
from connect_omero_IDR import establish_IDR_connection, establish_OMERO_connection
from tifffile import imwrite


IDR_URL = "https://idr.openmicroscopy.org"
ZARR_ENTITY_TYPE = "com.glencoesoftware.ngff:multiscales"
MAX_IDR_FILES = 20


# import environmental variables for OMERO login
usr = os.getenv("OMERO_USER")
psw = os.getenv("OMERO_PASSWORD")
uuid_key = os.getenv("UUID_SESSION_KEY")


# helper functions
def write_tsv(path, header, rows):
    with open(path, "w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def get_json(session, url):
    response = session.get(url)
    response.raise_for_status()
    return response.json()

# IDR block


# get image data Zarr/Images
def get_image_from_IDR(image_id, session, out_dir):
    r = get_json(session, f"{IDR_URL}/api/v0/m/images/{image_id}/")
    external_info = r["data"]["omero:details"].get("externalInfo") or {}
    is_zarr = external_info.get("EntityType") == ZARR_ENTITY_TYPE
    print(f"Is Zarr? {is_zarr}")
    # get the Zarr url if the dataset is a ZARR
    if is_zarr:
        zarr_url = external_info["Lsid"].replace("s3:", "https:").replace("?anonymous=true", "").rstrip("/")
        validator_url = f"https://ome.github.io/ome-ngff-validator/?source={zarr_url}"
        # write a TSV file with different informations
        out_tsv = os.path.join(out_dir, f"zarr_info_{image_id}.tsv")
        write_tsv(out_tsv, ["image_id", "zarr_url", "validator_url"],
                  [[image_id, zarr_url, validator_url]])
        print(f"Zarr info saved to {out_tsv}")
        # wite a HTML file with the link to the NGFF-validator
        out_html = os.path.join(out_dir, f"NGFF_validator_{image_id}.html")
        with open(out_html, "w", encoding="utf-8") as f:
            f.write(f"Open <a href='{validator_url}'>image in NGFF-validator</a>\n")
        print(f"HTML saved to {out_html}")
        return
    # if not Zarr get the original files only if there are at most MAX_IDR_FILES
    # (to avoid fetching too much data from IDR)
    download_urls = get_json(session, f"{IDR_URL}/image/{image_id}/download_urls/")["download_urls"]
    if len(download_urls) > MAX_IDR_FILES:
        out_tsv = os.path.join(out_dir, f"download_urls_{image_id}.tsv")
        write_tsv(out_tsv, ["image_name", "url"],
                  [[url.split("/")[-1], url] for url in download_urls])
        print(f"More than {MAX_IDR_FILES} files: URLs saved to {out_tsv}")
        return
    for tif_url in download_urls:
        image_path = os.path.join(out_dir, tif_url.split("/")[-1])
        # stream to disk instead of holding the whole file in memory
        with session.get(tif_url, stream=True) as response:
            response.raise_for_status()
            with open(image_path, mode="wb") as file:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    file.write(chunk)
        print(f"Saved {image_path}")


# get annotation (key-value pairs) from IDR
def get_annotation_from_IDR(obj_id, obj_type, session, out_dir):
    url = f"{IDR_URL}/webclient/api/annotations/?type=map&{obj_type}={obj_id}"
    rows = [[a["ns"], key, value]
            for a in get_json(session, url)["annotations"]
            for key, value in a["values"]]
    out_tsv = os.path.join(out_dir, f"annotations_{obj_type}_{obj_id}.tsv")
    write_tsv(out_tsv, ["namespace", "key", "value"], rows)
    print(f"{len(rows)} annotations saved to {out_tsv}")


# get bulk (tables) from IDR
def get_bulk_annotation_from_IDR(screen_id, session, out_dir):
    url = f"{IDR_URL}/webclient/api/annotations/?type=file&screen={screen_id}"
    json_data = get_json(session, url)
    # look for name = "bulk_annotations"
    file_id = None
    for ann in json_data["annotations"]:
        if ann.get("file", {}).get("name") == "bulk_annotations":
            file_id = ann["file"]["id"]
            break
    if file_id is None:
        print(f"No bulk_annotations file found for screen {screen_id}")
        return
    print(f"Found bulk_annotations file {file_id}")
    limit = 1000
    offset = 0
    row_count = None
    data = []
    while row_count is None or offset < row_count:
        url = f"{IDR_URL}/webclient/omero_table/{file_id}/json/?limit={limit}&offset={offset}"
        json_data = get_json(session, url)
        row_count = json_data["meta"]["totalCount"]
        column_names = json_data["data"]["columns"]
        data.extend(json_data["data"]["rows"])
        offset += limit
        print(f"loaded {min(offset, row_count)}/{row_count} rows")
    df_ann = DataFrame(data, columns=column_names)
    out_tsv = os.path.join(out_dir, f"bulk_annotations_{screen_id}.tsv")
    df_ann.to_csv(out_tsv, sep="\t", index=False)
    print(f"Bulk annotations saved to {out_tsv} ({len(df_ann)} rows)")


# OMERO BLOCK

# fetch single image function
def fetch_image(image, channel, stack, time, tile, mode, out_dir):
    image_id = image.getId()
    pixels = image.getPrimaryPixels()
    size_z, size_c, size_t = image.getSizeZ(), image.getSizeC(), image.getSizeT()
    # the function can get planes, stacks or hypercubes which will be saved as TIFF
    if mode == "plane":
        planes = pixels.getPlane(stack, channel, time)
        axes = "YX"
        out_name = f"image_{image_id}_z{stack}_c{channel}_t{time}.tiff"
    elif mode == "stack":
        if tile is None:
            tile = (0, 0, pixels.getSizeX(), pixels.getSizeY())   # x, y, width, height of tile
        zct_list = [(iz, channel, time, tile) for iz in range(size_z)]
        planes = np.array(list(pixels.getTiles(zct_list)))
        axes = "ZYX"
        out_name = f"image_{image_id}_stack_c{channel}_t{time}.tiff"
    else:  # hypercube
        zct_list = [(z, c, t) for z in range(size_z) for c in range(size_c) for t in range(size_t)]
        planes = np.array(list(pixels.getPlanes(zct_list)))
        # restore the Z, C, T dimensions (same order as zct_list)
        planes = planes.reshape(size_z, size_c, size_t, *planes.shape[-2:])
        axes = "ZCTYX"
        out_name = f"image_{image_id}_hypercube.tiff"
    out = os.path.join(out_dir, out_name)
    imwrite(out, planes, metadata={"axes": axes})
    print(f"Saved {out} (shape: {planes.shape}, axes: {axes})")


# function to fetch one or more images calling fetch_image
def get_images_from_OMERO(image_ids, channel, stack, time, tile, mode, conn, out_dir):
    missing = []
    for i, image_id in enumerate(image_ids, start=1):
        print(f"[{i}/{len(image_ids)}] Image {image_id}")
        image = conn.getObject("Image", image_id)
        if image is None:
            # skip missing ids instead of aborting the whole batch
            print(f"Warning: no image with id {image_id}, skipping", file=sys.stderr)
            missing.append(image_id)
            continue
        fetch_image(image, channel, stack, time, tile, mode, out_dir)
    if missing:
        raise ValueError(f"{len(missing)} image id(s) not found: {missing}")


# function to fetch a dataset with all images
def get_images_from_dataset(dataset_id, channel, stack, time, tile, mode, conn, out_dir):
    dataset = conn.getObject("Dataset", dataset_id)
    if dataset is None:
        raise ValueError(f"No dataset with id {dataset_id}")
    for image in dataset.listChildren():
        fetch_image(image, channel, stack, time, tile, mode, out_dir)


# function to get the original data
def get_original_data_from_OMERO(conn, host, port, image_ids, dataset_id, out_dir):
    if image_ids:
        targets = [("Image", image_id) for image_id in image_ids]
    else:
        targets = [("Dataset", dataset_id)]
    # reuse the session of the open connection: no password on the command line
    session_key = conn.getEventContext().sessionUuid
    for obj_type, obj_id in targets:
        # one subfolder per object to avoid file name collisions between filesets
        target_dir = os.path.join(out_dir, f"{obj_type.lower()}_{obj_id}")
        os.makedirs(target_dir, exist_ok=True)
        subprocess.run(["omero", "download", "-s", host, "-p", str(port), "-k", session_key,
                        f"{obj_type}:{obj_id}", target_dir], check=True)
        print(f"Original files of {obj_type} {obj_id} saved to {target_dir}")


def parse_args():
    examples = """examples:
  image        python %(prog)s image 12345 -o results
  annotation   python %(prog)s annotation 12345 image
  bulk         python %(prog)s bulk 777
  omero plane  python %(prog)s omero 54321 --mode plane -z 0 -c 1 -t 0
  omero multi  python %(prog)s omero 54321 54322 54323 --mode plane -o results
  omero stack  python %(prog)s omero 54321 --mode stack --tile 50 50 10 10
  omero cube   python %(prog)s omero 54321 --mode hypercube
  omero ds     python %(prog)s omero --dataset 777 --mode plane
  omero orig   python %(prog)s omero 54321 54322 --original -o originals
"""
    parser = argparse.ArgumentParser(
        description="Fetch images and map annotations from IDR, or pixel data from OMERO.",
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # options shared by all subcommands (can be given after the subcommand)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--output-dir", default=".", help="Destination folder for all output files (created if missing, default: .)")

    sub = parser.add_subparsers(dest="command", required=True)

    p_img = sub.add_parser("image", parents=[common], help="Download an IDR image (Zarr images: info saved to a TSV)")
    p_img.add_argument("image_id", type=int, help="OMERO image id")
    p_anno = sub.add_parser("annotation", parents=[common], help="Save IDR map annotations to a TSV file")
    p_anno.add_argument("obj_id", type=int, help="Object id")
    p_anno.add_argument("obj_type", help="Object type for the API query (e.g. image, project)")
    p_bulk = sub.add_parser("bulk", parents=[common], help="Save a screen's bulk_annotations table to a TSV file")
    p_bulk.add_argument("screen_id", type=int, help="Screen id with a bulk_annotations file annotation")
    p_omero = sub.add_parser("omero", parents=[common], help="Fetch pixel data from one or more OMERO images")
    p_omero.add_argument("--original", action="store_true",
                         help="Download the original files with 'omero download'")
    p_omero.add_argument("image_ids", type=int, nargs="*", default=[],
                         help="One or more OMERO image ids (omit when using --dataset)")
    p_omero.add_argument("--dataset", type=int, default=None,
                         help="Dataset id: fetch all images in the dataset")
    p_omero.add_argument("--mode", choices=["plane", "stack", "hypercube"], default="plane",
                         help="plane: one ZCT plane | stack: tiled Z-stack | hypercube: all ZCT planes")
    p_omero.add_argument("-z", "--stack", type=int, default=0, help="Z index (plane mode)")
    p_omero.add_argument("-c", "--channel", type=int, default=0, help="Channel index (plane and stack modes)")
    p_omero.add_argument("-t", "--time", type=int, default=0, help="Time index (plane and stack modes)")
    p_omero.add_argument("--tile", type=int, nargs=4, metavar=("X", "Y", "W", "H"), help="Tile region x y width height (stack mode, default: full frame)")
    p_omero.add_argument("--host", default="idr.openmicroscopy.org", help="OMERO server host")
    p_omero.add_argument("--port", type=int, default=4064, help="OMERO server port")

    return parser.parse_args()


def main(args):
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    if args.command == "omero":
        if not (uuid_key or (usr and psw)):
            raise EnvironmentError("Set UUID_SESSION_KEY, or both OMERO_USER and OMERO_PASSWORD")
        if bool(args.image_ids) == (args.dataset is not None):
            raise ValueError("Provide either one or more image ids or --dataset (not both)")
        # drop duplicate ids while keeping the given order
        image_ids = list(dict.fromkeys(args.image_ids))
        tile = tuple(args.tile) if args.tile is not None else None
        conn = establish_OMERO_connection(uuid_key, usr, psw, args.host, args.port)
        try:
            if args.original:
                get_original_data_from_OMERO(conn, args.host, args.port, image_ids,
                                             args.dataset, out_dir)
            elif args.dataset is not None:
                get_images_from_dataset(args.dataset, args.channel, args.stack, args.time,
                                        tile, args.mode, conn, out_dir)
            else:
                get_images_from_OMERO(image_ids, args.channel, args.stack, args.time,
                                      tile, args.mode, conn, out_dir)
        finally:
            conn.close()
    else:
        session = establish_IDR_connection()
        try:
            if args.command == "image":
                get_image_from_IDR(args.image_id, session, out_dir)
            elif args.command == "annotation":
                get_annotation_from_IDR(args.obj_id, args.obj_type, session, out_dir)
            else:  # bulk
                get_bulk_annotation_from_IDR(args.screen_id, session, out_dir)
        finally:
            session.close()


if __name__ == "__main__":
    args = parse_args()
    try:
        main(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
