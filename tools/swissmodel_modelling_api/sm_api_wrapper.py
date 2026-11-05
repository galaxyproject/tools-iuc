"""Wrapper for the SWISS-MODEL API."""

import argparse
import functools
import json
import os
import re
import sys
import time
from urllib.parse import urlsplit

import requests


def _api_error_detail(response, max_len=2000):
    """Return the error message the SWISS-MODEL API sent with a response.

    The API explains why a request was rejected (e.g. an invalid target
    sequence) in the response body; without it, users only see a status code.
    """
    try:
        detail = json.dumps(response.json(), indent=2)
    except ValueError:
        detail = response.text.strip()
    if not detail:
        return ""
    if len(detail) > max_len:
        detail = detail[:max_len] + "... [truncated]"
    return f"\nResponse from SWISS-MODEL:\n{detail}"


class _TransientError(Exception):
    """A request failed in a way that may go away when trying again."""


def _send(func, *args, **kwargs):
    """Send a request, raise _TransientError for failures worth retrying."""
    # HTTP status codes worth retrying, usually the server being overloaded
    # or down for maintenance
    retry_status_codes = (429, 500, 502, 503, 504)
    try:
        response = func(*args, **kwargs)
    except (
        requests.exceptions.ConnectionError,
        requests.exceptions.Timeout,
        requests.exceptions.ChunkedEncodingError,
        requests.exceptions.ContentDecodingError,
    ) as exc:
        raise _TransientError(type(exc).__name__) from exc
    if response.status_code in retry_status_codes:
        raise _TransientError(f"HTTP {response.status_code}")

    return response


class _SmApiWhisperer:
    """Parent class for talking to the SWISS-MODEL API."""

    PROJECT_TYPE = ""
    # Polling intervals in seconds
    WAIT_QUEUEING = 31
    WAIT_RUNNING = 17
    WAIT_RUNNING_LONG = 127
    # After this many seconds in RUNNING, switch to the long interval
    RUNNING_LONG_AFTER = 60 * 60
    # Stop waiting after this many seconds in RUNNING
    RUNNING_MAX = 24 * 60 * 60
    # Seconds to wait before trying a failed request again
    WAIT_RETRY = 31
    # Give up when the server can not be reached for this many seconds
    OUTAGE_MAX = 60 * 60

    def __init__(self, targets, token, project_title="Untitled Project"):
        self.project_id = None
        self.project_title = project_title
        self.targets = targets
        self.token = token

    def get_json_payload(self):
        """Needs to be implemented per project type."""
        raise NotImplementedError

    def _retry(self, send, what):
        """Call send() again and again while it raises _TransientError.

        Gives up when the server could not be reached for OUTAGE_MAX seconds.
        Messages must not contain "error:" or "exception:", Galaxy would mark
        the job as failed even if a retry succeeds.
        """
        failing_since = None
        while True:
            try:
                return send()
            except _TransientError as exc:
                if failing_since is None:
                    failing_since = time.monotonic()
                if time.monotonic() - failing_since > self.OUTAGE_MAX:
                    print(
                        "SWISS-MODEL seems to be temporarily unavailable, "
                        + f"could not reach it while {what} for "
                        + f"{self.OUTAGE_MAX // 60} minutes ({exc})",
                        file=sys.stderr,
                    )
                    sys.exit(3)
                print(
                    f"Could not reach SWISS-MODEL while {what} ({exc}), "
                    + f"trying again in {self.WAIT_RETRY} s",
                    file=sys.stderr,
                )
                time.sleep(self.WAIT_RETRY)

    def submit_request(self):
        """Send off a request to the SM API."""
        json_payload = self.get_json_payload()
        json_payload["project_title"] = self.project_title
        # Submitting the same job twice should reuse the existing project, so
        # this can be retried even if the first request reached the server.
        response = self._retry(
            functools.partial(
                _send,
                requests.post,
                f"https://swissmodel.expasy.org/{self.PROJECT_TYPE}",
                headers={"Authorization": f"Token {self.token}"},
                json=json_payload,
                timeout=60,
            ),
            "submitting the modelling job",
        )
        if response.ok is not True:
            raise RuntimeError(
                f"Submitting modelling job failed ({response.status_code})"
                + _api_error_detail(response)
            )
        self.project_id = response.json()["project_id"]
        print(f"SWISS-MODEL project ID: {self.project_id}")

        return response.status_code

    def _fetch_details(self):
        """Get the current state of the project from the server."""
        response = _send(
            requests.get,
            f"https://swissmodel.expasy.org/project/{self.project_id}/"
            + "models/full-details/",
            headers={"Authorization": f"Token {self.token}"},
            timeout=360,
        )
        if response.ok is not True:
            raise RuntimeError(
                f"Polling project '{self.project_id}' failed "
                + f"({response.status_code})"
                + _api_error_detail(response)
            )
        try:
            return response.json()
        except ValueError as exc:
            # e.g. a maintenance page instead of JSON
            raise _TransientError("response is not JSON") from exc

    def wait(self):
        """Poll the API for job to be finished."""
        details = None
        running_since = None
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
            details = self._retry(
                self._fetch_details, f"polling project '{self.project_id}'"
            )
            # Update the status
            status = details["status"].upper()
            if status in ["COMPLETED", "FAILED"]:
                break
            # Wait for some time before the next request, depending on the
            # status. Lots of jobs finish within ~30s, so the first wait in
            # the queue may already be the last one.
            if status in ["INITIALISED", "QUEUEING"]:
                time.sleep(self.WAIT_QUEUEING)
                continue
            if status != "RUNNING":
                raise RuntimeError(
                    f"Unknown status '{status}' for project "
                    + f"'{self.project_id}'"
                )
            # The time limit starts when the job starts running.
            if running_since is None:
                running_since = time.monotonic()
            running_for = time.monotonic() - running_since
            if running_for > self.RUNNING_MAX:
                print(
                    "SWISS-MODEL modelling job for project "
                    + f"'{self.project_id}' is still running after "
                    + f"{self.RUNNING_MAX // 3600} hours, giving up waiting.",
                    file=sys.stderr,
                )
                sys.exit(5)
            if running_for < self.RUNNING_LONG_AFTER:
                time.sleep(self.WAIT_RUNNING)
            else:
                time.sleep(self.WAIT_RUNNING_LONG)

        return details

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
            fname = f"model_{os.path.basename(urlsplit(url).path)}"
            response = self._retry(
                functools.partial(_send, requests.get, url, timeout=360),
                f"downloading {fname}",
            )
            if response.ok is not True:
                raise RuntimeError(
                    f"Fetching {file_type} output failed ("
                    + f"{response.status_code})."
                    + _api_error_detail(response)
                )
            try:
                os.mkdir(os.path.join(outdir, file_type))
            except FileExistsError:
                pass
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
    """SM alignment project."""

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


def _defastarise_sequences(sequences):
    """In case some of the targets carry FastA headers, remove them.

    Sequences come from single line text fields of the Galaxy form, where
    Galaxy replaces '>' by '__gt__'.
    """
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

    return ["".join(seq.split()) for seq in targets]


def _read_fasta(fasta_file):
    """Read a Multi-FastA file into (header, sequence) records.

    Whitespace is removed from sequences, records without sequence lines are
    kept with an empty sequence.
    """
    records = []

    # utf-8-sig drops a byte order mark, as written by some Windows editors
    with open(fasta_file, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                records.append([line[1:], []])
            elif line:
                if not records:
                    records.append(["", []])
                records[-1][1].append("".join(line.split()))

    return [(header, "".join(seq)) for header, seq in records]


def _defastarise_targets(fasta_file):
    """Parse sequences from a Multi-FastA file, removing headers.

    A record without sequence raises ValueError, dropping it would silently
    turn a heteromer into a different complex.
    """
    targets = []
    for header, seq in _read_fasta(fasta_file):
        if not seq:
            raise ValueError(f"Target FASTA record '{header}' has no sequence")
        targets.append(seq)

    return targets


def _parse_alignment_file(fasta_file):
    """Read a target-template alignment for alignment mode.

    The file holds two records, target first, then template. The template
    header starts with the SMTL ID '<PDB ID>.<assembly ID>.<chain name>' and
    may carry 'offset=X' for the template sequence offset.
    """
    records = _read_fasta(fasta_file)
    if len(records) != 2:
        raise ValueError(
            "Alignment file must contain exactly 2 sequences, found "
            + f"{len(records)}"
        )
    (_, target_sequence), (tpl_header, template_sequence) = records
    if not target_sequence or not template_sequence:
        raise ValueError("Alignment file has an empty sequence")
    tpl_words = tpl_header.split()
    smtl_id = tpl_words[0] if tpl_words else ""
    match = re.match(r"^(\d[A-Za-z0-9]{3})\.(\d+)\.([A-Za-z0-9]+)$", smtl_id)
    if match is None:
        raise ValueError(
            "Template header must start with an SMTL ID like '1crn.1.A', "
            + f"found '{smtl_id}'"
        )
    offset = 0
    for word in tpl_words[1:]:
        if word.startswith("offset="):
            value = word.removeprefix("offset=")
            if not re.match(r"^\d+$", value):
                raise ValueError(
                    "Template offset must be a non-negative integer, found "
                    + f"'{word}'"
                )
            offset = int(value)
            break
        if word.lower().startswith("offset"):
            raise ValueError(
                "Template offset must be written as 'offset=X', found "
                + f"'{word}'"
            )

    return {
        "target_sequence": target_sequence,
        "template_sequence": template_sequence,
        "pdb_id": match.group(1),
        "assembly_id": int(match.group(2)),
        "auth_asym_id": match.group(3),
        "template_seqres_offset": offset,
    }


def _add_target_sequences(subparser):
    """Add positional target sequences to a single job subparser."""
    subparser.add_argument(
        "target_sequences",
        help="Target sequence(s) to be modelled; multiple sequences for "
        + "hetero modelling",
        metavar="<SEQUENCE>",
        nargs="+",
    )


def _add_target_sequences_file(subparser):
    """Add the target FastA file to a batch subparser."""
    subparser.add_argument(
        "target_sequences_file",
        help="FastA file with target sequence to be modelled; "
        + "Multi-FastA file for multiple target sequences",
        metavar="<[Multi-]FASTA FILE>",
    )


def _add_template_file(subparser):
    """Add the template file option to a user template subparser."""
    subparser.add_argument(
        "-f",
        "--template-file",
        help="PDB formatted file to serve as template for modelling",
        metavar="<PDB FILE>",
        required=True,
    )


def _parse_args(argv=None):
    """Get command line arguments."""
    # Options shared by all subcommands
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar="<OUTPUT DIRECTORY>",
    )
    common.add_argument(
        "-d",
        "--project-title",
        help="Title for the modelling project",
        metavar="<TITLE>",
    )
    common.add_argument(
        "-m",
        "--no-modelcif",
        help="Do not download models in ModelCIF format.",
        default=False,
        action="store_true",
    )
    common.add_argument(
        "-l",
        "--fetch-pdb",
        help="Download models in PDB legacy format.",
        default=False,
        action="store_true",
    )

    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True, metavar="<MODE>")
    # both modes offer the same project types
    project_types = {
        "dest": "project_type",
        "required": True,
        "metavar": "<PROJECT TYPE>",
    }
    single = modes.add_parser(
        "single", help="Single job, sequences given on the command line"
    ).add_subparsers(**project_types)
    batch = modes.add_parser(
        "batch", help="Batch job, sequences given as files"
    ).add_subparsers(**project_types)

    # single job mode
    _add_target_sequences(
        single.add_parser("automodel", parents=[common], help="Automodel mode")
    )
    sgl_aln = single.add_parser(
        "alignment", parents=[common], help="Alignment mode"
    )
    # Sequences are options, aligned sequences may start with a gap ('-')
    sgl_aln.add_argument(
        "--target-sequence",
        help="Target sequence of the alignment, give as "
        + "'--target-sequence=<SEQUENCE>'",
        metavar="<SEQUENCE>",
        required=True,
    )
    sgl_aln.add_argument(
        "-t",
        "--template-sequence",
        help="Template sequence of the alignment, give as "
        + "'--template-sequence=<SEQUENCE>'",
        metavar="<SEQUENCE>",
        required=True,
    )
    # ToDo: do we need the offset from the user? Doesn't interactive alignment
    #       mode compute it?
    sgl_aln.add_argument(
        "-o",
        "--template-seqres-offset",
        help="Offset of the template sequence segment compared to the full "
        + "template sequence",
        metavar="<NUMBER>",
        type=int,
        required=True,
    )
    sgl_aln.add_argument(
        "-p",
        "--pdb-id",
        help="PDB ID (SMTL ID) for the template used in alignment mode",
        metavar="<PDB ID>",
        required=True,
    )
    sgl_aln.add_argument(
        "-c",
        "--auth-asym-id",
        help="The chain name to be used in alignment mode",
        metavar="<CHAIN NAME>",
        required=True,
    )
    sgl_aln.add_argument(
        "-a",
        "--assembly-id",
        help="ID of the assembly of the SMTL template to be used in alignment "
        + "mode",
        metavar="<NUMBER>",
        type=int,
        required=True,
    )
    sgl_usr = single.add_parser(
        "usertemplate", parents=[common], help="User template mode"
    )
    _add_target_sequences(sgl_usr)
    _add_template_file(sgl_usr)

    # batch mode
    _add_target_sequences_file(
        batch.add_parser("automodel", parents=[common], help="Automodel mode")
    )
    batch.add_parser(
        "alignment", parents=[common], help="Alignment mode"
    ).add_argument(
        "alignment_file",
        help="FastA file with the target-template alignment, target first, "
        + "template header starting with the SMTL ID, e.g. '1crn.1.A', "
        + "optionally followed by 'offset=X'",
        metavar="<ALIGNMENT FILE>",
    )
    bch_usr = batch.add_parser(
        "usertemplate", parents=[common], help="User template mode"
    )
    _add_target_sequences_file(bch_usr)
    _add_template_file(bch_usr)

    opts = parser.parse_args(argv)

    # string arguments can not be empty
    for name in [
        "outdir",
        "pdb_id",
        "auth_asym_id",
        "template_file",
        "target_sequences_file",
        "alignment_file",
    ]:
        value = getattr(opts, name, None)
        if value is not None and len(value) == 0:
            print(
                f"Argument '{name}' can not be an empty string",
                file=sys.stderr,
            )
            sys.exit(2)
    if opts.project_title is not None and len(opts.project_title) == 0:
        print(
            "Option '--project-title' can not have an empty string as value",
            file=sys.stderr,
        )
        sys.exit(2)

    return opts


def _load_inputs(opts):
    """Read the sequences for the modelling job.

    Returns the keyword arguments for the whisperer classes besides token,
    project title and template file. Raises ValueError for unusable input.
    """
    if opts.mode == "batch" and opts.project_type == "alignment":
        aln = _parse_alignment_file(opts.alignment_file)
        aln["targets"] = [aln.pop("target_sequence")]
        return aln

    if opts.mode == "single" and opts.project_type == "alignment":
        targets = _defastarise_sequences([opts.target_sequence])
    elif opts.mode == "single":
        targets = _defastarise_sequences(opts.target_sequences)
    else:
        targets = _defastarise_targets(opts.target_sequences_file)
    if len(targets) == 0 or not all(targets):
        raise ValueError("No target sequence given")
    inputs = {"targets": targets}

    if opts.project_type == "alignment":
        template_sequence = _defastarise_sequences([opts.template_sequence])
        if not template_sequence[0]:
            raise ValueError("No template sequence given")
        inputs.update(
            {
                "template_sequence": template_sequence[0],
                "template_seqres_offset": opts.template_seqres_offset,
                "pdb_id": opts.pdb_id,
                "auth_asym_id": opts.auth_asym_id,
                "assembly_id": opts.assembly_id,
            }
        )

    return inputs


def _main():
    """Run as script."""
    opts = _parse_args()
    # check the input before anything else, so the user does not need to
    # wait for the API to learn about broken input
    try:
        inputs = _load_inputs(opts)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(2)

    token = os.getenv("SWISSMODEL_API_TOKEN")
    if not token:
        print(
            "SWISS-MODEL token is not provided in credentials!",
            file=sys.stderr,
        )
        sys.exit(1)
    # determine class
    whsprr = None
    if opts.project_type.lower() == "automodel":
        whsprr = _AutoModelWhisperer(
            inputs["targets"], token, project_title=opts.project_title
        )
    elif opts.project_type.lower() == "alignment":
        whsprr = _AlignmentWhisperer(
            inputs["targets"],
            token,
            inputs["template_sequence"],
            inputs["template_seqres_offset"],
            inputs["pdb_id"],
            inputs["auth_asym_id"],
            inputs["assembly_id"],
            project_title=opts.project_title,
        )
    elif opts.project_type.lower() == "usertemplate":
        whsprr = _UserTemplateWhisperer(
            inputs["targets"],
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
    # A job without models is still COMPLETED, FAILED means the modelling
    # pipeline crashed on the server.
    if response["status"] == "FAILED":
        print(
            "SWISS-MODEL modelling job failed for project "
            + f"'{whsprr.project_id}'. This is an error on the SWISS-MODEL "
            + "server, not a problem with your input. Please try again "
            + "later; if it keeps failing, report the project ID to "
            + "SWISS-MODEL.",
            file=sys.stderr,
        )
        sys.exit(4)
    whsprr.fetch_results(
        response,
        opts.outdir,
        fetch_modelcif=not opts.no_modelcif,
        fetch_pdb=opts.fetch_pdb,
    )

    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint FASTA
