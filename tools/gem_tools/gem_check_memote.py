import argparse
import subprocess


def __main__():
    parser = argparse.ArgumentParser(
        prog="MemoteReport",
        description="This program takes a memote snapshot report of a GEM",
        epilog="Adding an epilog, but doubt it's needed.",
    )
    parser.add_argument(
        "-m",
        "--cb_model_location",
        dest="cb_model_location",
        action="store", type=str,
        default=None,
        required=True,
        help="The model to use."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        action="store",
        type=str,
        default=None,
        required=True,
        help="The output file."
    )

    output_location = parser.parse_args().out_file
    cb_model_location = parser.parse_args().cb_model_location

    command = (
        f"memote report snapshot --filename '{output_location}' "
        f"'{cb_model_location}'"
    )
    try:
        # Execute the command
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Output the result
        print("Command executed successfully.")
        print("Output:\n", result.stdout)
    except subprocess.CalledProcessError as e:
        # Handle errors in the subprocess
        if "Error:" in e.stderr or "error:" in e.stderr:
            raise Exception(e.stderr)
        else:
            print(e.stderr)


if __name__ == "__main__":
    __main__()
