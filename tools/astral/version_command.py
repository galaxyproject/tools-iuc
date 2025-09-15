#!/usr/bin/env python
"""Modules to execute shell commands (through child process) and regex"""
import re
import subprocess


def get_astral_version():
    """Function that will parse the Astral version from `astral --help`."""
    try:
        # run the `astral --help` command
        result = subprocess.run(
            ['astral', '--help'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        # save the output
        output = result.stdout + result.stderr
        # Regex pattern that matches the version string from the help message
        version_pattern = re.compile(r'This is ASTRAL version (\d+\.\d+\.\d+)')
        # search for the version pattern in the output
        match = version_pattern.search(output)
        if match:
            # extract and return the version string
            return match.group(1)
        print("Version information not found in `astral --help` output.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        return None
    except FileNotFoundError as e:
        print(f"Command not found: {e}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# call the function and print the version
version = get_astral_version()
if version:
    print(version)
