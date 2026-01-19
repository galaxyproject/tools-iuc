import argparse
import requests
import yaml
import sys
from collections import defaultdict

def download_yaml_template(workflow, dims, biapy_version=""):
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
    
    # Use .get() to avoid KeyError if workflow is unexpected
    dir_name = template_dir_map.get(workflow)
    if not dir_name:
        raise ValueError(f"Unknown workflow: {workflow}")

    template_name = f"{dir_name}/{dims.lower()}_{dir_name}.yaml"
    url = f"https://raw.githubusercontent.com/BiaPyX/BiaPy/refs/tags/v{biapy_version}/templates/{template_name}"
    
    print(f"Downloading YAML template from {url}")
    try:
        response = requests.get(url, timeout=10) # Added timeout
        response.raise_for_status() # Automatically raises HTTPError for 4xx/5xx
        return yaml.safe_load(response.text) or {}
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download template. {e}")
        sys.exit(1) # Exit gracefully rather than crashing with a stack trace

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
        
        ndim = "3D" if args.dims == "3d" else "2D"
        as_stack = args.dims in ["2d_stack", "2d"]

        config = download_yaml_template(workflow_type, ndim, biapy_version=args.biapy_version)
        
        # Initialization using setdefault to prevent KeyErrors
        config.setdefault("PROBLEM", {})
        config["PROBLEM"].update({"TYPE": workflow_type, "NDIM": ndim})
        
        config.setdefault("TEST", {})["ANALIZE_2D_IMGS_AS_3D_STACK"] = as_stack

        # Handle MODEL and PATHS
        model_cfg = config.setdefault("MODEL", {})
        if args.model_source == "biapy":
            model_cfg["SOURCE"] = "biapy"
            is_loading = bool(args.model)
            model_cfg["LOAD_CHECKPOINT"] = is_loading
            model_cfg["LOAD_MODEL_FROM_CHECKPOINT"] = is_loading
            if is_loading:
                config.setdefault("PATHS", {})["CHECKPOINT_FILE"] = args.model
        elif args.model_source == "bmz":
            model_cfg["SOURCE"] = "bmz"
            model_cfg.setdefault("BMZ", {})["SOURCE_MODEL_ID"] = args.model
        elif args.model_source == "torchvision":
            model_cfg["SOURCE"] = "torchvision"
            model_cfg["TORCHVISION_MODEL_NAME"] = args.model

        # PATCH_SIZE Logic
        obj_size_map = {
            "0-25": (256, 256), "25-100": (256, 256),
            "100-200": (512, 512), "200-500": (512, 512), "500+": (1024, 1024),
        }
        obj_size = obj_size_map[args.obj_size]
        
        obj_slices_map = {"": -1, "1-5": 5, "5-10": 10, "10-20": 20, "20-60": 40, "60+": 80}
        obj_slices = obj_slices_map.get(args.obj_slices, -1)

        if ndim == "2D":
            patch_size = obj_size + (args.img_channel,)
        else:
            if obj_slices == -1:
                print("Error: For 3D problems, obj_slices must be specified.")
                sys.exit(1)
            patch_size = (obj_slices,) + obj_size + (args.img_channel,)
        
        config.setdefault("DATA", {})["PATCH_SIZE"] = str(patch_size)
        config["DATA"]["REFLECT_TO_COMPLETE_SHAPE"] = True

    else:
        if not args.input_config_path:
            print("Error: Input configuration path must be specified.")
            sys.exit(1)
        try:
            with open(args.input_config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f) or {}
        except FileNotFoundError:
            print(f"Error: File {args.input_config_path} not found.")
            sys.exit(1)

    # Global overrides (Train/Test)
    config.setdefault("TRAIN", {})
    config.setdefault("DATA", {})
    
    if args.raw_train:
        config["TRAIN"]["ENABLE"] = True
        config["DATA"].setdefault("TRAIN", {}).update({
            "PATH": args.raw_train,
            "GT_PATH": args.gt_train
        })
    else:
        config["TRAIN"]["ENABLE"] = False

    test_cfg = config.setdefault("TEST", {})
    if args.test_raw_path:
        test_cfg["ENABLE"] = True
        data_test = config["DATA"].setdefault("TEST", {})
        data_test["PATH"] = args.test_raw_path
        data_test["LOAD_GT"] = bool(args.test_gt_path)
        if args.test_gt_path:
            data_test["GT_PATH"] = args.test_gt_path
    else:
        test_cfg["ENABLE"] = False

    config.setdefault("MODEL", {})["OUT_CHECKPOINT_FORMAT"] = "safetensors"
    
    # Final cleanup and save
    config = tuple_to_list(config)
    with open(args.out_config_path, 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"Success: YAML configuration written to {args.out_config_path}")

if __name__ == "__main__":
    main()