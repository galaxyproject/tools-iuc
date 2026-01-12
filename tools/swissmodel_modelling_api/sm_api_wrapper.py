"""Wrapper for the SWISS-MODEL API."""

import argparse
import json
import os
import sys
import time
from urllib.parse import urlsplit

import requests


class _SmApiWhisperer:
    """Parent class for talking to the SWISS-MODEL API."""

    PROJECT_TYPE = ""

    def __init__(self, targets, token, project_title="Untitled Project"):
        self.project_id = None
        self.project_title = project_title
        self.targets = targets
        self.token = token

    def get_json_payload(self):
        """Needs to be implemented per project type."""
        raise NotImplementedError

    def submit_request(self):
        """Send off a request to the SM API."""
        json_payload = self.get_json_payload()
        json_payload["project_title"] = self.project_title
        try:
            response = requests.post(
                f"https://swissmodel.expasy.org/{self.PROJECT_TYPE}",
                headers={"Authorization": f"Token {self.token}"},
                json=json_payload,
                timeout=60,
            )
        except requests.exceptions.ConnectTimeout:
            print(
                "SWISS-MODEL seems to temporarily unavailable",
                file=sys.stderr,
            )
            sys.exit(3)
        if response.ok is not True:
            raise RuntimeError(
                f"Submitting modelling job failed ({response.status_code})"
            )
        self.project_id = response.json()["project_id"]

        return response.status_code

    def wait(self):
        """Poll the API for job to be finished."""
        response = None
        # Wait at the end, there is a chance that this project is already
        # available from cache.
        while True:
            # Update the status from the server
            # response = requests.get(
            #     f"https://swissmodel.expasy.org/project/{self.project_id}/"
            #     + "models/summary/",
            #     headers={"Authorization": f"Token {self.token}"},
            #     timeout=360,
            # )
            response = requests.get(
                f"https://swissmodel.expasy.org/project/{self.project_id}/"
                + "models/full-details/",
                headers={"Authorization": f"Token {self.token}"},
                timeout=360,
            )
            # Update the status
            status = response.json()["status"]
            if status.upper() in ["COMPLETED", "FAILED"]:
                break
            # Wait for some time before the next request
            time.sleep(17)

        return response.json()

    def fetch_results(
        self, response_object, output_dir, fetch_modelcif=True, fetch_pdb=True
    ):
        """Get results of the modelling job."""

        def _store_model_json(model_json, outdir):
            fname = f"model_{model_json['model_id']}.json"
            with open(
                os.path.join(outdir, "JSON", fname), "w", encoding="utf8"
            ) as jfh:
                json.dump(model_json, jfh)

        def _fetch_file(url, file_type, outdir):
            response = requests.get(url, timeout=360)
            if response.ok is not True:
                raise RuntimeError(
                    f"Fetching {file_type} output failed ("
                    + f"{response.status_code})."
                )
            try:
                os.mkdir(os.path.join(outdir, file_type))
            except FileExistsError:
                pass
            fname = f"model_{os.path.basename(urlsplit(url).path)}"
            with open(os.path.join(outdir, file_type, fname), "wb") as mfh:
                for chunk in response.iter_content(chunk_size=8192):
                    mfh.write(chunk)

        # make sure a JSON directory exists
        os.mkdir(os.path.join(output_dir, "JSON"))
        if response_object["status"] == "COMPLETED":
            for model in response_object["models"]:
                _store_model_json(model, output_dir)
                if fetch_modelcif:
                    _fetch_file(model["modelcif_url"], "ModelCIF", output_dir)
                if fetch_pdb:
                    _fetch_file(model["coordinates_url"], "PDB", output_dir)


class _AutoModelWhisperer(_SmApiWhisperer):
    """SM automodel project."""

    PROJECT_TYPE = "automodel"

    def get_json_payload(self):
        """Payload for automodel mode."""
        return {"target_sequences": self.targets}


class _AlignmentWhisperer(_SmApiWhisperer):
    """SM alignemt project."""

    PROJECT_TYPE = "alignment"

    def __init__(
        self,
        targets,
        token,
        template_sequence,
        template_seqres_offset,
        pdb_id,
        auth_asym_id,
        assembly_id,
        project_title="Untitled Project",
    ):
        # Not sure how to reduce the number of arguments as they are required
        # by the API, so make an exception in Pylint.
        # pylint: disable=too-many-arguments,too-many-positional-arguments
        """Initialise alignment mode, add mode-specific info to the method."""
        super().__init__(targets, token, project_title=project_title)
        self.assembly_id = assembly_id
        self.auth_asym_id = auth_asym_id
        self.pdb_id = pdb_id.lower()
        self.template_seqres_offset = template_seqres_offset
        self.template_sequence = template_sequence

    def get_json_payload(self):
        """Payload for alignment mode."""

        return {
            "assembly_id": self.assembly_id,
            "auth_asym_id": self.auth_asym_id,
            "pdb_id": self.pdb_id,
            "target_sequences": self.targets,
            "template_seqres_offset": self.template_seqres_offset,
            "template_sequence": self.template_sequence,
        }


class _UserTemplateWhisperer(_SmApiWhisperer):
    """SM user-template project."""

    PROJECT_TYPE = "user_template"

    def __init__(
        self,
        targets,
        token,
        template_file,
        project_title="Untitled Project",
    ):
        """Initialise user template mode."""
        super().__init__(targets, token, project_title=project_title)
        self.template_file = template_file

    def get_json_payload(self):
        """Payload for user upload mode."""
        with open(self.template_file, encoding="utf8") as tfh:
            template_coordinates = tfh.read()

        return {
            "project_title": self.project_title,
            "target_sequences": self.targets,
            "template_coordinates": template_coordinates,
        }


def _defastarise_targets(sequences):
    """In case some of the targets carry FastA headers, remove them."""
    targets = []
    for seq in sequences:
        seq = seq.split(" ")
        if len(seq) > 1:
            if seq[0].strip().startswith((">", "__gt__")):
                targets.append("".join(seq[1:]))
            else:
                targets.append("".join(seq))
        else:
            targets.extend(seq)

    return targets


def _parse_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-d",
        "--project-title",
        help="Title for the modelling project",
        metavar="<TITLE>",
    )
    parser.add_argument(
        "-m",
        "--no-modelcif",
        help="Do not download models in ModelCIF format.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-l",
        "--fetch-pdb",
        help="Download models in PDB legacy format.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--template-sequence",
        help="The template sequence used for alignment mode",
        metavar="<SEQUENCE>",
    )
    # ToDo: do we need the offset from the user? Doesn't interactive alignment
    #       mode compute it?
    parser.add_argument(
        "-o",
        "--template-seqres-offset",
        help="Offset of the template sequence segment compared to the full "
        + "template sequence",
        metavar="<NUMBER>",
        type=int,
    )
    parser.add_argument(
        "-p",
        "--pdb-id",
        help="PDB ID (SMTL ID) for the template used in alignment mode",
        metavar="<PDB ID>",
    )
    parser.add_argument(
        "-c",
        "--auth-asym-id",
        help="The chain name to be used in alignment mode",
        metavar="<CHAIN NAME>",
    )
    parser.add_argument(
        "-a",
        "--assembly-id",
        help="ID of the assembly of the SMTL template to be used in alignment "
        + "mode",
        metavar="<NUMBER>",
        type=int,
    )
    parser.add_argument(
        "-f",
        "--template-file",
        help="PDB formatted file to serve as template for modelling",
        metavar="<PDB FILE>",
    )
    parser.add_argument(
        "project_type",
        choices=("alignment", "automodel", "usertemplate"),
        help="Kind of project ('alignmet', 'automodel', 'usertemplate')",
        metavar="<PROJECT TYPE>",
    )
    metas = {
        "outdir": "<OUTPUT DIRECTORY>",
        "target_sequences": "<SEQUENCE[S]>",
    }
    parser.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar=metas["outdir"],
    )
    parser.add_argument(
        "target_sequences",
        help="Target sequence to be modelled; to add multiple sequences, "
        + "delimit with a space",
        metavar=metas["target_sequences"],
        nargs=argparse.REMAINDER,
    )

    opts = parser.parse_args()

    # Make sure arguments for the different modelling modes are there
    req_opts = {
        "alignment": [
            "assembly_id",
            "auth_asym_id",
            "pdb_id",
            "template_seqres_offset",
            "template_sequence",
        ],
        "automodel": [],
        "usertemplate": ["template_file"],
    }
    # check mandatory arguments
    for req in req_opts[opts.project_type]:
        value = getattr(opts, req)
        if value is None:
            print(
                f"Option '--{req.replace('_', '-')}' missing for "
                + f"'{opts.project_type}' mode",
                file=sys.stderr,
            )
            sys.exit(2)
        if isinstance(value, str) and len(value) == 0:
            print(
                f"Option '--{req.replace('_', '-')}' can not be an empty "
                + "string",
                file=sys.stderr,
            )
            sys.exit(2)
    # check positional arguments
    for req, mta in metas.items():
        value = getattr(opts, req)
        if isinstance(value, str):
            if len(value) == 0:
                print(
                    f"Argument of '{mta}' can not be an empty string",
                    file=sys.stderr,
                )
                sys.exit(2)
        elif isinstance(value, list):
            if len(value) == 0 or not all(value):
                print(
                    f"Argument of '{mta}' can not be an empty",
                    file=sys.stderr,
                )
                sys.exit(2)
        else:
            raise RuntimeError(
                f"Value with unknown type '{type(value).__name__}' found for "
                + f"'{mta}'"
            )
    # check optional & positional arguments
    for opt in ["project_title"]:
        value = getattr(opts, opt)
        if value is not None and len(value) == 0:
            print(
                f"Option '--{opt.replace('_', '-')}' can not have an empty "
                + "string as value",
                file=sys.stderr,
            )
            sys.exit(2)

    return opts


def _main():
    """Run as script."""
    opts = _parse_args()

    token = os.getenv("SWISSMODEL_API_TOKEN")
    if not token:
        print(
            "SWISS-MODEL token is not provided in credentials!",
            file=sys.stderr,
        )
        sys.exit(1)
    target_sequences = _defastarise_targets(opts.target_sequences)
    # determine class
    whsprr = None
    if opts.project_type.lower() == "automodel":
        whsprr = _AutoModelWhisperer(
            target_sequences, token, project_title=opts.project_title
        )
    elif opts.project_type.lower() == "alignment":
        template_sequence = _defastarise_targets([opts.template_sequence])
        assert len(template_sequence) == 1
        template_sequence = template_sequence[0]
        whsprr = _AlignmentWhisperer(
            target_sequences,
            token,
            template_sequence,
            opts.template_seqres_offset,
            opts.pdb_id,
            opts.auth_asym_id,
            opts.assembly_id,
            project_title=opts.project_title,
        )
    elif opts.project_type.lower() == "usertemplate":
        whsprr = _UserTemplateWhisperer(
            target_sequences,
            token,
            opts.template_file,
            project_title=opts.project_title,
        )
    else:
        raise RuntimeError(
            f"Not a suitable project type: '{opts.project_type}'"
        )
    # run the modelling job and wait for it to finish
    whsprr.submit_request()
    response = whsprr.wait()
    whsprr.fetch_results(
        response,
        opts.outdir,
        fetch_modelcif=not opts.no_modelcif,
        fetch_pdb=opts.fetch_pdb,
    )

    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint
