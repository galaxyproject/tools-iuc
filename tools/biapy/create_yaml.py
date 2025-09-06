import argparse

import requests
import yaml


def download_yaml_template(workflow, dims, biapy_version=""):
    """
    Download a YAML template for a specific workflow and dimensions.

    Parameters:
        workflow (str): The workflow type.
        dims (str): The dimensions (e.g., 2d, 3d).
        biapy_version (str): The BiaPy version to use.

    Returns:
        dict: The YAML template as a dictionary.
    """
    template_dir_map = {
        "SEMANTIC_SEG": "semantic_segmentation",
        "INSTANCE_SEG": "instance_segmentation",
        "DETECTION": "detection",
        "DENOISING": "denoising",
        "SUPER_RESOLUTION": "super-resolution",
        "CLASSIFICATION": "classification",
        "SELF_SUPERVISED": "self-supervised",
        "IMAGE_TO_IMAGE": "image-to-image",
    }
    template_name = (
        template_dir_map[workflow]
        + "/"
        + dims.lower() + "_" + template_dir_map[workflow] + ".yaml"
    )

    url = (
        f"https://raw.githubusercontent.com/BiaPyX/BiaPy/"
        f"refs/tags/v{biapy_version}/templates/{template_name}"
    )
    print(f"Downloading YAML template from {url}")
    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(
            f"Failed to download YAML template: {response.status_code}"
        )
    return yaml.safe_load(response.text)


def tuple_to_list(obj):
    """Convert tuples to lists recursively."""
    if isinstance(obj, tuple):
        return list(obj)
    if isinstance(obj, dict):
        return {k: tuple_to_list(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [tuple_to_list(v) for v in obj]
    return obj


def main():
    parser = argparse.ArgumentParser(
        description="Generate a YAML configuration from given arguments."
    )
    parser.add_argument(
        '--input_config_path', default='', type=str,
        help="Input configuration file to reuse"
    )
    parser.add_argument(
        '--new_config', action='store_true',
        help="Whether to create a new config or reuse an existing one."
    )
    parser.add_argument(
        '--out_config_path', required=True, type=str,
        help="Path to save the generated YAML configuration."
    )
    parser.add_argument(
        '--workflow', default='semantic', type=str,
        choices=['semantic', 'instance', 'detection', 'denoising',
                 'sr', 'cls', 'sr2', 'i2i'],
    )
    parser.add_argument(
        '--dims', default='2d', type=str,
        choices=['2d_stack', '2d', '3d'],
        help="Number of dimensions for the problem"
    )
    parser.add_argument(
        '--obj_slices', default='', type=str,
        choices=['', '1-5', '5-10', '10-20', '20-60', '60+'],
        help="Number of slices for the objects in the images"
    )
    parser.add_argument(
        '--obj_size', default='0-25', type=str,
        choices=['0-25', '25-100', '100-200', '200-500', '500+'],
        help="Size of the objects in the images"
    )
    parser.add_argument(
        '--img_channel', default=1, type=int,
        help="Number of channels in the input images"
    )
    parser.add_argument(
        '--model_source', default='biapy',
        choices=['biapy', 'bmz', 'torchvision'],
        help="Source of the model."
    )
    parser.add_argument(
        '--model', default='', type=str,
        help=("Path to the model file if using a pre-trained model "
              "from BiaPy or name of the model within BioImage "
              "Model Zoo or TorchVision.")
    )
    parser.add_argument(
        '--raw_train', default='', type=str,
        help="Path to the training raw data."
    )
    parser.add_argument(
        '--gt_train', default='', type=str,
        help="Path to the training ground truth data."
    )
    parser.add_argument(
        '--test_raw_path', default='', type=str,
        help="Path to the testing raw data."
    )
    parser.add_argument(
        '--test_gt_path', default='', type=str,
        help="Path to the testing ground truth data."
    )
    parser.add_argument(
        '--biapy_version', default='', type=str,
        help="BiaPy version to use."
    )

    args = parser.parse_args()

    if args.new_config:
        workflow_map = {
            "semantic": "SEMANTIC_SEG",
            "instance": "INSTANCE_SEG",
            "detection": "DETECTION",
            "denoising": "DENOISING",
            "sr": "SUPER_RESOLUTION",
            "cls": "CLASSIFICATION",
            "sr2": "SELF_SUPERVISED",
            "i2i": "IMAGE_TO_IMAGE",
        }
        workflow_type = workflow_map[args.workflow]

        if args.dims == "2d_stack":
            ndim = "2D"
            as_stack = True
        elif args.dims == "2d":
            ndim = "2D"
            as_stack = True
        elif args.dims == "3d":
            ndim = "3D"
            as_stack = False

        config = download_yaml_template(
            workflow_type, ndim, biapy_version=args.biapy_version
        )

        config["PROBLEM"]["TYPE"] = workflow_type
        config["PROBLEM"]["NDIM"] = ndim
        config["TEST"]["ANALIZE_2D_IMGS_AS_3D_STACK"] = as_stack

        if args.model_source == "biapy":
            config["MODEL"]["SOURCE"] = "biapy"
            if args.model:
                config["MODEL"]["LOAD_CHECKPOINT"] = True
                config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = True
                config.setdefault("PATHS", {})
                config["PATHS"]["CHECKPOINT_FILE"] = args.model
            else:
                config["MODEL"]["LOAD_CHECKPOINT"] = False
                config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
        elif args.model_source == "bmz":
            config["MODEL"]["SOURCE"] = "bmz"
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
            config.setdefault("MODEL", {}).setdefault("BMZ", {})
            config["MODEL"]["BMZ"]["SOURCE_MODEL_ID"] = args.model
        elif args.model_source == "torchvision":
            config["MODEL"]["SOURCE"] = "torchvision"
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
            config["MODEL"]["TORCHVISION_MODEL_NAME"] = args.model

        obj_size_map = {
            "0-25": (256, 256),
            "25-100": (256, 256),
            "100-200": (512, 512),
            "200-500": (512, 512),
            "500+": (1024, 1024),
        }
        obj_size = obj_size_map[args.obj_size]

        obj_slices_map = {
            "": -1,
            "1-5": 5,
            "5-10": 10,
            "10-20": 20,
            "20-60": 40,
            "60+": 80,
        }
        obj_slices = obj_slices_map[args.obj_slices]
        if config["PROBLEM"]["NDIM"] == "2D":
            config["DATA"]["PATCH_SIZE"] = obj_size + (args.img_channel,)
        else:
            assert obj_slices != -1, (
                "For 3D problems, obj_slices must be specified."
            )
            config["DATA"]["PATCH_SIZE"] = (
                (obj_slices,) + obj_size + (args.img_channel,)
            )
        config["DATA"]["PATCH_SIZE"] = str(config["DATA"]["PATCH_SIZE"])
    else:
        assert args.input_config_path, (
            "Input configuration path must be specified when not "
            "creating a new config."
        )
        with open(args.input_config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        if args.model:
            config["MODEL"]["SOURCE"] = "biapy"
            config["MODEL"]["LOAD_CHECKPOINT"] = True
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = True
            config.setdefault("PATHS", {})
            config["PATHS"]["CHECKPOINT_FILE"] = args.model
        else:
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False

    if args.raw_train:
        config["TRAIN"]["ENABLE"] = True
        config["DATA"]["TRAIN"]["PATH"] = args.raw_train
        config["DATA"]["TRAIN"]["GT_PATH"] = args.gt_train
    else:
        config["TRAIN"]["ENABLE"] = False

    if args.test_raw_path:
        config["TEST"]["ENABLE"] = True
        config["DATA"]["TEST"]["PATH"] = args.test_raw_path
        if args.test_gt_path:
            config["DATA"]["TEST"]["LOAD_GT"] = True
            config["DATA"]["TEST"]["GT_PATH"] = args.test_gt_path
        else:
            config["DATA"]["TEST"]["LOAD_GT"] = False
    else:
        config["TEST"]["ENABLE"] = False

    # Always use safetensors in Galaxy
    config["MODEL"]["OUT_CHECKPOINT_FORMAT"] = "safetensors"

    config = tuple_to_list(config)

    with open(args.out_config_path, 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"YAML configuration written to {args.out_config_path}")


if __name__ == "__main__":
    main()
