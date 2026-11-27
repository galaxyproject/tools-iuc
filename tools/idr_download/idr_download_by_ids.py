# import libraries
import argparse
import csv
import os
import sys
from pathlib import Path


import numpy as np
from connect_omero_IDR import establish_IDR_connection, establish_OMERO_connection, IDR_URL
from tifffile import imwrite, memmap


ZARR_ENTITY_TYPE = "com.glencoesoftware.ngff:multiscales"
MAX_IDR_FILES = 20
BIGTIFF_THRESHOLD = 2**32 - 2**25


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


# (getFileInChunks) had the same "open file, write every chunk" loop.
def write_chunks(path, chunks):
    with open(path, "wb") as file:
        for chunk in chunks:
            file.write(chunk)


def run_batch(items, label, func, skip_failed):
    failed = []
    for i, item in enumerate(items, start=1):
        print(f"[{i}/{len(items)}] {label} {item}", flush=True)
        try:
            func(item)
        except Exception as e:
            print(f"Warning: {label} {item} failed: {e}", file=sys.stderr, flush=True)
            failed.append(item)
    if failed:
        message = f"{len(failed)} of {len(items)} {label}(s) failed: {failed}"
        if skip_failed:
            print(f"Warning: {message} (skipped, --skip-failed)", file=sys.stderr)
        else:
            raise RuntimeError(message)


# IDR block


# get image data Zarr/Images
def get_image_from_IDR(image_id, session, out_dir):
    r = get_json(session, f"{IDR_URL}/api/v0/m/images/{image_id}/")
    external_info = r["data"]["omero:details"].get("externalInfo") or {}
    is_zarr = external_info.get("EntityType") == ZARR_ENTITY_TYPE
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
            write_chunks(image_path, response.iter_content(chunk_size=1024 * 1024))
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
    out_tsv = os.path.join(out_dir, f"bulk_annotations_{screen_id}.tsv")
    write_tsv(out_tsv, column_names, data)
    print(f"Bulk annotations saved to {out_tsv} {row_count} rows)")


# OMERO BLOCK


def write_planes(out, planes, shape, axes):
    tiff = None
    for index, plane in zip(np.ndindex(*shape[:-2]), planes):
        if tiff is None:
            # dtype is only known once the first plane arrives
            nbytes = int(np.prod(shape, dtype=np.int64)) * plane.dtype.itemsize
            print(f"Writing {out}: shape {shape}, {plane.dtype}, {nbytes / 2**20:.1f} MiB", flush=True)
            tiff = memmap(out, shape=shape, dtype=plane.dtype, metadata={"axes": axes},
                          bigtiff=nbytes > BIGTIFF_THRESHOLD)
        tiff[index] = plane
    if tiff is None:
        raise ValueError(f"No planes returned for {out}")
    tiff.flush()
    del tiff


# fetch single image function
def fetch_image(image, channel, stack, time, tile, mode, out_dir):
    image_id = image.getId()
    pixels = image.getPrimaryPixels()
    size_z, size_c, size_t = image.getSizeZ(), image.getSizeC(), image.getSizeT()
    # the function can get planes, stacks or hypercubes which will be saved as TIFF
    if mode == "plane":
        out = os.path.join(out_dir, f"image_{image_id}_z{stack}_c{channel}_t{time}.tiff")
        plane = pixels.getPlane(stack, channel, time)
        imwrite(out, plane, metadata={"axes": "YX"})
        print(f"Saved {out} (shape: {plane.shape}, axes: YX)")
    elif mode == "stack":
        if tile is None:
            tile = (0, 0, pixels.getSizeX(), pixels.getSizeY())   # x, y, width, height of tile
        zct_list = [(iz, channel, time, tile) for iz in range(size_z)]
        out = os.path.join(out_dir, f"image_{image_id}_stack_c{channel}_t{time}.tiff")
        write_planes(out, pixels.getTiles(zct_list), (size_z, tile[3], tile[2]), "ZYX")
    else:  # hypercube
        zct_list = [(z, c, t) for z in range(size_z) for c in range(size_c) for t in range(size_t)]
        out = os.path.join(out_dir, f"image_{image_id}_hypercube.tiff")
        write_planes(out, pixels.getPlanes(zct_list),
                     (size_z, size_c, size_t, pixels.getSizeY(), pixels.getSizeX()), "ZCTYX")


def get_image(conn, image_id):
    image = conn.getObject("Image", image_id)
    if image is None:
        raise ValueError(f"no image with id {image_id}")
    return image


# function to fetch one or more images calling fetch_image
def get_dataset_images(conn, dataset_id):
    dataset = conn.getObject("Dataset", dataset_id)
    if dataset is None:
        raise ValueError(f"No dataset with id {dataset_id}")
    return {image.getId(): image for image in dataset.listChildren()}


def resolve_images(conn, image_ids, dataset_id):
    if dataset_id is None:
        return dict.fromkeys(image_ids)
    images = get_dataset_images(conn, dataset_id)
    print(f"Dataset {dataset_id}: {len(images)} image(s)")
    return images


# function to fetch a dataset with all images
def download_original_files(image, fileset, out_dir):
    folder = f"fileset_{fileset.getId()}" if fileset is not None else f"image_{image.getId()}"
    target = Path(out_dir) / folder
    target.mkdir(exist_ok=True)
    for original_file in image.getImportedImageFiles():
        # name/stem/suffix contain no folders, so nothing is written outside the folder
        name = Path(original_file.getName())
        dest = target / name.name
        # same name in different subfolders: add _1, _2, ... before the extension
        n = 1
        while dest.exists():
            dest = target / f"{name.stem}_{n}{name.suffix}"
            n += 1
        write_chunks(dest, original_file.getFileInChunks())
        print(f"Downloaded {dest}")


# function to get the original data
def download_original(image, out_dir, seen_filesets):
    fileset = image.getFileset()
    if fileset is not None:
        if fileset.getId() in seen_filesets:
            print(f"Fileset {fileset.getId()} already downloaded with image "
                  f"{seen_filesets[fileset.getId()]}, skipping image {image.getId()}")
            return
        seen_filesets[fileset.getId()] = image.getId()
    download_original_files(image, fileset, out_dir)


def parse_args():
    examples = """examples:
  image        python %(prog)s image 12345 -o results
  image multi  python %(prog)s image 12345 12346 --skip-failed -o results
  annotation   python %(prog)s annotation 12345 image
  bulk         python %(prog)s bulk 777
  omero plane  python %(prog)s omero 54321 --mode plane -z 0 -c 1 -t 0
  omero multi  python %(prog)s omero 54321 54322 54323 --mode plane -o results
  omero stack  python %(prog)s omero 54321 --mode stack --tile 50 50 10 10
  omero cube   python %(prog)s omero 54321 --mode hypercube
  omero ds     python %(prog)s omero --dataset 777 --mode plane
  omero orig   python %(prog)s omero 54321 54322 --original -o originals
  omero ds-org python %(prog)s omero --dataset 777 --original -o originals
"""
    parser = argparse.ArgumentParser(
        description="Fetch images and map annotations from IDR, or pixel data from OMERO.",
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # options shared by all subcommands (can be given after the subcommand)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--output-dir", default=".", help="Destination folder for all output files (created if missing, default: .)")

    batch = argparse.ArgumentParser(add_help=False)
    batch.add_argument("--skip-failed", action="store_true",
                       help="Warn about ids that fail instead of exiting with an error")

    sub = parser.add_subparsers(dest="command", required=True)

    p_img = sub.add_parser("image", parents=[common, batch], help="Download IDR images (Zarr images: info saved to a TSV)")
    p_img.add_argument("image_ids", type=int, nargs="+", help="One or more IDR image ids")
    p_anno = sub.add_parser("annotation", parents=[common], help="Save IDR map annotations to a TSV file")
    p_anno.add_argument("obj_id", type=int, help="Object id")
    p_anno.add_argument("obj_type", help="Object type for the API query (e.g. image, project)")
    p_bulk = sub.add_parser("bulk", parents=[common], help="Save a screen's bulk_annotations table to a TSV file")
    p_bulk.add_argument("screen_id", type=int, help="Screen id with a bulk_annotations file annotation")
    p_omero = sub.add_parser("omero", parents=[common, batch], help="Fetch pixel data from one or more OMERO images")
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
    image_ids = list(dict.fromkeys(getattr(args, "image_ids", [])))

    if args.command == "omero":
        if not (uuid_key or (usr and psw)):
            raise EnvironmentError("Set UUID_SESSION_KEY, or both OMERO_USER and OMERO_PASSWORD")
        if bool(args.image_ids) == (args.dataset is not None):
            raise ValueError("Provide either one or more image ids or --dataset (not both)")
        tile = tuple(args.tile) if args.tile is not None else None
        conn = establish_OMERO_connection(uuid_key, usr, psw, args.host, args.port)
        try:
            images = resolve_images(conn, image_ids, args.dataset)
            if args.original:
                label = "Image (original)"
                seen_filesets = {}

                def task(image):
                    download_original(image, out_dir, seen_filesets)
            else:
                label = "Image"

                def task(image):
                    fetch_image(image, args.channel, args.stack, args.time, tile,
                                args.mode, out_dir)

            def process(image_id):
                image = images[image_id]
                task(image if image is not None else get_image(conn, image_id))
            run_batch(list(images), label, process, args.skip_failed)
        finally:
            conn.close(hard=not uuid_key)
    else:
        session = establish_IDR_connection()
        try:
            if args.command == "image":
                run_batch(image_ids, "IDR image",
                          lambda image_id: get_image_from_IDR(image_id, session, out_dir),
                          args.skip_failed)
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
