import argparse
import yaml
import requests

def download_yaml_template(template_name="template_biax.yaml"):
    url = f"https://raw.githubusercontent.com/BiaPyX/BiaPy/master/templates/{template_name}"
    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download YAML template: {response.status_code}")
    return yaml.safe_load(response.text)

def main():
    parser = argparse.ArgumentParser(description="Generate a YAML configuration from given arguments.")
    parser.add_argument('--input_config_path', default='', type=str, help="Input configuration file to reuse")
    parser.add_argument('--new_config', action='store_true', help="Whether to create a new config or reuse an existing one.")
    parser.add_argument('--out_config_path', required=True, type=str, help="Path to save the generated YAML configuration.")
    parser.add_argument('--workflow', default='', type=str, choices=['semantic', 'instance', 'detection', 'denoising', 'sr', 'cls', 'sr2', 'i2i'],)
    parser.add_argument('--dims', default='', type=str, help="Number of dimensions for the problem", choices=['2d_stack', '2d', '3d'])
    parser.add_argument('--obj_slices', default='', type=str, help="Number of slices for the objects in the images", choices=['1-5', '5-10', '10-20', '20-60', '60+'])
    parser.add_argument('--obj_size', default='', type=str, help="Size of the objects in the images", choices=['0-25', '25-100', '100-200', '200-500', '500+'])
    parser.add_argument('--model_source', default='biapy', choices=['biapy', 'bmz', 'torchvision'], help="Source of the model.")
    parser.add_argument('--model', default='', type=str, help="Path to the model file if using a pre-trained model from BiaPy or name of the model within BioImage Model Zoo or TorchVision.")
    parser.add_argument('--train_raw_path', default='', type=str, help="Path to the training raw data.")
    parser.add_argument('--train_gt_path', default='', type=str, help="Path to the training ground truth data.")
    parser.add_argument('--test_raw_path', default='', type=str, help="Path to the testing raw data.")
    parser.add_argument('--test_gt_path', default='', type=str, help="Path to the testing ground truth data.")

    args = parser.parse_args()

    if args.new_config:
        # Load the template
        config = download_yaml_template()
    else:
        assert args.input_config_path != '', "Input configuration path must be specified when not creating a new config."
        # Load the existing configuration file
        with open(args.input_config_path, 'r') as f:
            config = yaml.safe_load(f)

    # Q1
    if args.dims == "2d_stack":
        config["PROBLEM"]["NDIM"] = "2D"
        config["TEST"]["ANALIZE_2D_IMGS_AS_3D_STACK"] = True
    elif args.dims == "2d":
        config["PROBLEM"]["NDIM"] = "2D"
        config["TEST"]["ANALIZE_2D_IMGS_AS_3D_STACK"] = True
    elif args.dims == "3d":
        config["PROBLEM"]["NDIM"] = "3D"
        config["TEST"]["ANALIZE_2D_IMGS_AS_3D_STACK"] = False

    if args.new_config:
        assert args.workflow != "", "Workflow must be specified when creating a new config."
        # Q2
        # Map input workflow values from UI to BiaPy workflow constants
        workflow_map = {
            "semantic": "SEMANTIC_SEG",
            "instance": "INSTANCE_SEG",
            "detection": "DETECTION",
            "denoising": "DENOISING",
            "sr": "SUPER_RESOLUTION",
            "cls": "CLASSIFICATION",
            "sr2": "SELF_SUPERVISED",  # assuming this is meant to represent self-supervised restoration
            "i2i": "IMAGE_TO_IMAGE"
        }
        workflow_type = workflow_map[args.workflow]
        config["PROBLEM"]["TYPE"] = workflow_type

        # Q3, Q4 and Q5
        assert args.model_source != "", "Model source must be specified."
        if args.model_source == "biapy":
            config["MODEL"]["SOURCE"] = "biapy"
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
            if "PATHS" not in config:
                config["PATHS"] = {}
            config["PATHS"]["CHECKPOINT_FILE"] = args.model
        elif args.model_source == "bmz":
            config["MODEL"]["SOURCE"] = "bmz"
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
            if "BMZ" not in config["MODEL"]:
                config["MODEL"]["BMZ"] = {}
            config["MODEL"]["BMZ"]["SOURCE_MODEL_ID"] = args.model
        elif args.model_source == "torchvision":
            config["MODEL"]["SOURCE"] = "torchvision"
            config["MODEL"]["LOAD_CHECKPOINT"] = False
            config["MODEL"]["LOAD_MODEL_FROM_CHECKPOINT"] = False
            config["MODEL"]["TORCHVISION_MODEL_NAME"] = args.model

        # Q6
        assert args.obj_size != "", "Object size must be specified."
        obj_size_map = {
            "0-25": (256, 256),
            "25-100": (256, 256),
            "100-200": (512, 512),
            "200-500": (512, 512),
            "500+": (1024, 1024),
        }
        obj_size = obj_size_map[args.obj_size]

        # Q7
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
            config["DATA"]["PATCH_SIZE"] = obj_size
        else:
            assert obj_slices != -1, "For 3D problems, obj_slices must be specified."
            config["DATA"]["PATCH_SIZE"] = (obj_slices,) + obj_size

    # Q8, Q9, Q10, Q11 and Q12
    if args.train_raw_path != "":
        config["TRAIN"]["ENABLE"] = True
        config["DATA"]["TRAIN"]["PATH"] = args.train_raw_path
        config["DATA"]["TRAIN"]["GT_PATH"] = args.train_gt_path
    else:
        config["TRAIN"]["ENABLE"] = False
    if args.test_raw_path != "":
        config["TEST"]["ENABLE"] = True
        config["DATA"]["TEST"]["PATH"] = args.test_raw_path
        if args.test_gt_path != "":
            config["DATA"]["TEST"]["LOAD_GT"] = True
            config["DATA"]["TEST"]["GT_PATH"] = args.test_gt_path
    else:
        config["TEST"]["ENABLE"] = False
    

    # Save the YAML to out_config_path
    with open(args.out_config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"YAML configuration written to {args.out_config_path}")

if __name__ == "__main__":
    main()