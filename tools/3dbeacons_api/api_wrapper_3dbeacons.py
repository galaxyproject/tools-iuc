"""Wrapper for the 3D-Beacons API."""

import argparse
import collections
import csv
import functools
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import urlsplit

import requests

# Known providers (lowercase)
PROVIDERS = [
    "alphafill",
    "alphafold",
    "hegelab",
    "isoformio",
    "levylab",
    "modelarchive",
    "pdbe",
    "ped",
    "sasbdb",
    "swissmodel",
]

# Known parameter columns for batch mode
KNOWN_PARAMS = {
    "provider",
    "exclude_provider",
    "template",
    "sequence_range",
    "uniprot_checksum",
}


def _validate_provider(value):
    """Validate provider name. Returns lowercase value."""
    value_lower = value.lower()
    if value_lower not in PROVIDERS:
        raise ValueError(
            f"invalid provider '{value}' - must be one of: "
            + f"{', '.join(PROVIDERS)}"
        )
    return value_lower


def _validate_template(value):
    """Validate template PDB ID format. Returns lowercase value."""
    value_lower = value.lower()
    if not re.match(r"^[a-z0-9]{4}([._-]?\d+[._-]?[A-Za-z])?$", value_lower):
        raise ValueError(
            f"invalid template '{value}' - must be 4-character PDB ID with "
            + "optional chain/assembly (e.g., '1ake', '1ake.1.A', '1ake1A')"
        )
    return value_lower


def _validate_sequence_range(value):
    """Validate sequence range format."""
    if not re.match(r"^\d+-\d+$", value):
        raise ValueError(
            f"invalid sequence_range '{value}' - must be in format "
            + "'start-end' (e.g., '1-100')"
        )
    return value


def _validate_checksum(value):
    """Validate CRC64 checksum format. Returns uppercase value."""
    value_upper = value.upper()
    if not re.match(r"^[A-F0-9]{16}$", value_upper):
        raise ValueError(
            f"invalid uniprot_checksum '{value}' - must be 16-character CRC64 "
            + "hex value (e.g., '1A2B3C4D5E6F7890')"
        )
    return value_upper


def _parse_args():
    """Get command line arguments."""
    # pylint: disable=too-many-branches
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Single query mode
    single_parser = subparsers.add_parser(
        "single",
        help="Query a single entry",
    )
    single_parser.add_argument(
        "qualifier",
        help="ID/AC to search for",
        metavar="<QUALIFIER>",
    )
    single_parser.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar="<OUTPUT DIRECTORY>",
        type=Path,
    )
    single_parser.add_argument(
        "--provider",
        choices=PROVIDERS,
        help="Request from provider only: '"
        + "', '".join(PROVIDERS)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--exclude-provider",
        choices=PROVIDERS,
        help="Exclude provider: '"
        + "', '".join(PROVIDERS)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--template",
        help="Template PDB ID (4-character code, e.g. '1ake') or SMTL ID with "
        + "assembly and chain name (e.g. 1ake.1.A)",
        metavar="<TEMPLATE>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--sequence-range",
        help="Specify a UniProt sequence residue range (format: 'start-end', "
        + "e.g., '1-100')",
        metavar="<RANGE>",
    )
    single_parser.add_argument(
        "--uniprot-checksum",
        help="CRC64 checksum of the UniProt sequence (16-character hex value)",
        metavar="<CHECKSUM>",
        type=str.upper,
    )

    # Batch mode
    batch_parser = subparsers.add_parser(
        "batch", help="Process multiple queries from a TSV file"
    )
    batch_parser.add_argument(
        "tsv_file",
        help="Path to TSV query file with columns: qualifier and optional "
        + "parameters (provider, exclude_provider, template, sequence_range, "
        + "uniprot_checksum)",
        metavar="<TSV_FILE>",
        type=Path,
    )
    batch_parser.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar="<OUTPUT_DIRECTORY>",
        type=Path,
    )
    opts = parser.parse_args()

    # Validate single-query mode
    if opts.mode == "single":
        # check qualifier
        if len(opts.qualifier) == 0:
            single_parser.error(
                "Argument of '<QUALIFIER>' can not be an empty string"
            )
        # provider takes precedence over exclude-provider (warning if conflict)
        if opts.provider is not None and opts.exclude_provider is not None:
            if opts.provider == opts.exclude_provider:
                print(
                    "Warning: --provider takes precedence over "
                    + f"--exclude-provider (both set to '{opts.provider}')",
                    file=sys.stderr,
                )
                opts.exclude_provider = None
        # validate template format (PDB ID with optional chain/assembly)
        if opts.template is not None:
            try:
                opts.template = _validate_template(opts.template)
            except ValueError as e:
                single_parser.error(str(e))
        # validate sequence range format
        if opts.sequence_range is not None:
            try:
                opts.sequence_range = _validate_sequence_range(
                    opts.sequence_range
                )
            except ValueError as e:
                single_parser.error(str(e))
        # validate CRC64 checksum format
        if opts.uniprot_checksum is not None:
            try:
                opts.uniprot_checksum = _validate_checksum(
                    opts.uniprot_checksum
                )
            except ValueError as e:
                single_parser.error(str(e))
    # Validate batch mode
    elif opts.mode == "batch":
        if not opts.tsv_file.is_file():
            batch_parser.error(f"Query file not found: '{opts.tsv_file}'")

    return opts


def _parse_query_file(filepath):
    """Parse TSV query file with flexible column order (streaming).

    Validates parameter values and warns on unknown columns.

    Args:
        filepath: Path to TSV query file

    Yields:
        Tuples of (qualifier, params_dict)

    Raises:
        ValueError: If file format is invalid or parameter validation fails
    """

    # pylint: disable=too-many-branches
    def _filtered_line(file_handle):
        """Generator that yields non-comment, non-empty lines."""
        for line in file_handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                yield stripped

    try:
        with open(filepath, encoding="ascii") as fh:
            reader = csv.DictReader(_filtered_line(fh), delimiter="\t")

            # Validate header before processing rows
            required = {"qualifier"}
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"Query file must have {required} column in header. "
                    f"Found: {', '.join(reader.fieldnames)}"
                )

            # Process data rows
            for row_num, row in enumerate(reader, start=1):
                qualifier = row.get("qualifier", "").strip()
                if not qualifier:
                    raise ValueError(f"Row {row_num}: empty qualifier value")

                # Build params dict with validation
                params = {}
                for key, value in row.items():
                    if key == "qualifier":
                        continue
                    if not value or not value.strip():
                        continue

                    # Warn on unknown columns
                    if key not in KNOWN_PARAMS:
                        print(
                            f"Warning: Row {row_num}: unknown column '{key}' "
                            + f"(value: '{value}') - ignored",
                            file=sys.stderr,
                        )
                        continue

                    # Validate known parameters
                    try:
                        if key == "provider":
                            params[key] = _validate_provider(value)
                        elif key == "exclude_provider":
                            params[key] = _validate_provider(value)
                        elif key == "template":
                            params[key] = _validate_template(value)
                        elif key == "sequence_range":
                            params[key] = _validate_sequence_range(value)
                        elif key == "uniprot_checksum":
                            params[key] = _validate_checksum(value)
                    except ValueError as e:
                        raise ValueError(f"Row {row_num}: {e}") from e

                # provider takes precedence over exclude-provider (warning if
                # conflict)
                if (
                    "provider" in params
                    and "exclude_provider" in params
                    and params["provider"] == params["exclude_provider"]
                ):
                    print(
                        f"Warning: Row {row_num}: provider takes precedence "
                        + "over exclude_provider (both set to "
                        + f"'{params['provider']}')",
                        file=sys.stderr,
                    )
                    del params["exclude_provider"]

                yield (qualifier, params)

    except FileNotFoundError as err:
        raise ValueError(f"Query file not found: '{filepath}'") from err
    except csv.Error as err:
        raise ValueError(f"Invalid TSV format in '{filepath}': {err}") from err


def _assemble_call(
    qualifier,
    provider,
    exclude_provider,
    template,
    sequence_range,
    uniprot_checksum,
):
    """Get the URL for the request assembled."""
    # pylint: disable=too-many-positional-arguments,too-many-arguments
    url = (
        "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/uniprot"
        + f"/summary/{qualifier}.json"
    )
    params = {}
    if provider is not None:
        params["provider"] = provider
    if exclude_provider is not None:
        params["exclude_provider"] = exclude_provider
    if template is not None:
        params["template"] = template
    if sequence_range is not None:
        params["range"] = sequence_range
    if uniprot_checksum is not None:
        params["uniprot_checksum"] = uniprot_checksum

    return url, params


class _TransientError(Exception):
    """A request failed in a way that may go away when trying again."""


class _UnreachableError(Exception):
    """A server could not be reached for too long."""


class _DownloadError(Exception):
    """A model file could not be downloaded."""


def _send(url, **kwargs):
    """Send a GET request, raise _TransientError for failures to retry."""
    # HTTP status codes worth retrying, usually the server being overloaded
    # or down for maintenance
    retry_status_codes = (429, 500, 502, 503, 504)
    try:
        # callers pass the timeout in kwargs
        response = requests.get(  # pylint: disable=missing-timeout
            url, **kwargs
        )
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


def _retry(send, what, outage_max):
    """Call send() again and again while it raises _TransientError.

    Raises _UnreachableError when the server could not be reached for
    outage_max seconds. Messages must not contain "error:" or "exception:",
    Galaxy would mark the job as failed even if a retry succeeds.
    """
    # Seconds to wait before trying a failed request again
    wait_retry = 31
    failing_since = None
    while True:
        try:
            return send()
        except _TransientError as exc:
            if failing_since is None:
                failing_since = time.monotonic()
            if time.monotonic() - failing_since > outage_max:
                raise _UnreachableError(
                    f"could not reach the server while {what} for "
                    + f"{outage_max // 60} minutes ({exc})"
                ) from exc
            print(
                f"Could not reach the server while {what} ({exc}), "
                + f"trying again in {wait_retry} s",
                file=sys.stderr,
            )
            time.sleep(wait_retry)


class _FailedDownloads:
    """Keep track of models that could not be downloaded."""

    # pylint: disable=too-few-public-methods

    def __init__(self):
        # providers whose server could not be reached, skipped from then on
        self.unreachable = set()
        self.per_provider = collections.Counter()

    def report(self):
        """Print a summary of the models not downloaded."""
        if not self.per_provider:
            return
        print(
            f"{sum(self.per_provider.values())} model(s) could not be "
            + "downloaded: "
            + ", ".join(
                f"{provider} ({count})"
                for provider, count in sorted(self.per_provider.items())
            ),
            file=sys.stderr,
        )


_TRANSLATE_TABLE = str.maketrans(".:", "__")


def _make_coordfile_name(mdl_info):
    """Assemble a name for a coordinate file."""
    # some providers deliver non-unique file name, use model identifier instead
    if mdl_info["provider"].upper() in ("SWISS-MODEL", "PED"):
        mdl_id = mdl_info["model_identifier"].translate(_TRANSLATE_TABLE)
    else:
        mdl_id = str(Path(urlsplit(mdl_info["model_url"]).path).name)
    fname = Path(
        f"{mdl_info['provider']}__{mdl_info['model_format'].upper()}"
        + f"__{mdl_id}"
    )
    # make sure file has a lowercase extension, Galaxy discovers files by it
    if fname.suffix:
        fname = fname.with_suffix(fname.suffix.lower())
    else:
        if mdl_info["model_format"].upper() == "MMCIF":
            fname = fname.with_suffix(".cif")
        else:
            raise NotImplementedError(
                f"Unknown model format '{mdl_info['model_format']}' found, "
                + "needs to be added to the wrapper script."
            )

    return fname, mdl_id


def _fetch_file(mdl_info, outdir):
    """Download a model file, raise _DownloadError if that does not work."""
    # Give up on a provider when its server can not be reached for this many
    # seconds
    outage_max = 5 * 60
    fname, mid = _make_coordfile_name(mdl_info)
    response = _retry(
        functools.partial(_send, mdl_info["model_url"], timeout=360),
        f"downloading {mdl_info['model_url']}",
        outage_max,
    )
    if not response.ok:
        raise _DownloadError(f"HTTP {response.status_code} {response.reason}")
    with open(outdir / fname, "wb") as mfh:
        for chunk in response.iter_content(chunk_size=8192):
            mfh.write(chunk)
    return Path(mid).stem


def _collect_results(response, output_dir, failed):
    """Get JSON data and coordinate files."""
    # write entry info as JSON
    out_path = output_dir / "INFO"
    out_path.mkdir(parents=True, exist_ok=True)
    with open(
        out_path / f"uniprot_{response['uniprot_entry']['ac']}.json",
        "w",
        encoding="utf8",
    ) as jfh:
        json.dump(response["uniprot_entry"], jfh)
    # get coordinate files & per file a summary JSON, a failing provider must
    # not stop fetching models from the others
    for mdl in response["structures"]:
        mdl = mdl["summary"]
        provider = mdl["provider"]
        if provider in failed.unreachable:
            print(
                f"Skipping model '{mdl['model_identifier']}' from {provider}, "
                + "its server could not be reached before",
                file=sys.stderr,
            )
            failed.per_provider[provider] += 1
            continue
        # coordinate file
        out_path = output_dir / "MODELS"
        out_path.mkdir(parents=True, exist_ok=True)
        try:
            mid = _fetch_file(mdl, out_path)
        except _UnreachableError as exc:
            print(
                f"Giving up on {provider}, {exc}, skipping its remaining "
                + "models",
                file=sys.stderr,
            )
            failed.unreachable.add(provider)
            failed.per_provider[provider] += 1
            continue
        except _DownloadError as exc:
            print(
                f"Could not download model '{mdl['model_identifier']}' from "
                + f"{provider} ({exc}), {mdl['model_url']}; exclude "
                + f"{provider} if this keeps happening",
                file=sys.stderr,
            )
            failed.per_provider[provider] += 1
            continue
        except NotImplementedError as exc:
            print(
                f"Skipping model '{mdl['model_identifier']}' from "
                + f"{provider}: {exc}",
                file=sys.stderr,
            )
            failed.per_provider[provider] += 1
            continue
        # summary file in SUMMARIES/provider/
        summaries_path = output_dir / "SUMMARIES"
        summaries_path.mkdir(parents=True, exist_ok=True)
        with open(
            summaries_path / f"{mdl['provider']}__{mid}.json",
            "w",
            encoding="utf8",
        ) as jfh:
            json.dump(mdl, jfh)


def _run_query(
    qualifier,
    out_dir,
    provider,
    exclude_provider,
    template,
    sequence_range,
    uniprot_checksum,
    failed,
):
    """Try to fetch entry data"""
    # pylint: disable=too-many-positional-arguments,too-many-arguments
    # Give up when the 3D-Beacons API can not be reached for this many seconds
    api_outage_max = 60 * 60
    call, api_params = _assemble_call(
        qualifier,
        provider,
        exclude_provider,
        template,
        sequence_range,
        uniprot_checksum,
    )
    try:
        response = _retry(
            functools.partial(
                _send,
                call,
                headers={"Accept": "application/json"},
                params=api_params,
                timeout=360,
            ),
            f"querying 3D-Beacons for '{qualifier}'",
            api_outage_max,
        )
    except _UnreachableError as exc:
        print(
            f"3D-Beacons seems to be temporarily unavailable, {exc}",
            file=sys.stderr,
        )
        sys.exit(3)
    # 3D-Beacons answers 404 for unknown entries and when no model matches the
    # filters, both are a valid result, not a failure
    if response.status_code == 404:
        print(f"No models found for '{qualifier}' with the given filters")
        return
    if not response.ok:
        # leave out the response text, it may contain "error:" which makes
        # Galaxy fail the job
        print(
            f"3D-Beacons rejected the query for '{qualifier}' (HTTP "
            + f"{response.status_code} {response.reason}), skipping it",
            file=sys.stderr,
        )
        return
    _collect_results(response.json(), out_dir, failed)


def _main():
    """Run as script."""
    opts = _parse_args()
    failed = _FailedDownloads()

    if opts.mode == "batch":
        # Batch mode: process TSV file (streaming, memory efficient)
        # Process each query
        queries_processed = False
        for qualifier, params in _parse_query_file(opts.tsv_file):
            queries_processed = True
            qual_outdir = opts.outdir / qualifier
            _run_query(
                qualifier,
                qual_outdir,
                params.get("provider"),
                params.get("exclude_provider"),
                params.get("template"),
                params.get("sequence_range"),
                params.get("uniprot_checksum"),
                failed,
            )

        if not queries_processed:
            print("Error: Query file contains no data rows", file=sys.stderr)
            sys.exit(1)
    else:
        # Single query mode
        _run_query(
            opts.qualifier,
            opts.outdir,
            opts.provider,
            opts.exclude_provider,
            opts.template,
            opts.sequence_range,
            opts.uniprot_checksum,
            failed,
        )

    failed.report()
    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint FASTA
