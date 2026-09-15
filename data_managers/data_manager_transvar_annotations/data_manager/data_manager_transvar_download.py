#!/usr/bin/env python
"""Download pre-built TransVar transcript annotation databases.

TransVar's own ``transvar config --download_anno`` writes both the databases
and its configuration file into the installed Python package directory, which
is read-only (and shared between jobs) on a Galaxy installation. This data
manager instead downloads the same pre-built ``.transvardb`` files into a
Galaxy managed directory and registers them in the ``transvar_annotations``
data table.

The download URLs are taken from the installed TransVar release itself, so
they stay in sync with the requirement version rather than being duplicated
here. Files are stored under the name of the ``transvar {g,c,p}anno`` command
line flag that consumes them (``ucsc.transvardb``, ``refseq.transvardb``, ...)
so that the wrapper does not have to know TransVar's internal file naming.
"""

import argparse
import json
import os
import shutil
import sys
import tempfile
import time
from urllib.error import HTTPError
from urllib.request import urlopen

# A full annotation set is ~35 files and a few hundred megabytes from a single
# academic server, so a transient failure part way through is expected.
ATTEMPTS = 3
TIMEOUT = 60

PROXY_HINT = (
    "If this Galaxy server reaches the internet through an HTTP proxy, note "
    "that Galaxy strips the environment of a job down to HOME, LC_CTYPE, PATH "
    "and TMPDIR, so http_proxy/https_proxy have to be set for the job "
    "destination in job_conf.yml."
)

# TransVar config key -> the "transvar {g,c,p}anno" flag that takes the
# database. Only annotation sets that are openly licensed and freely
# redistributable are offered. AceView is deliberately excluded: it is
# distributed for use with attribution but states no redistribution licence.
FLAG_BY_CONFIG_KEY = {
    "ensembl": "ensembl",
    "gencode": "gencode",
    "refseq": "refseq",
    "ccds": "ccds",
    "ucsc": "ucsc",
    "known_gene": "kg",
}

SOURCE_LABELS = {
    "ensembl": "Ensembl",
    "gencode": "GENCODE",
    "refseq": "RefSeq",
    "ccds": "CCDS",
    "ucsc": "UCSC RefGene",
    "kg": "UCSC knownGene",
}


def files_to_download(refversion, wanted_flags):
    """Return [(target_basename, url), ...] for the requested sources.

    Every database consists of one ``.transvardb`` file plus a set of index
    files that TransVar locates by appending suffixes to that name, so the
    whole group has to be renamed consistently.
    """
    from transvar import config as transvar_config

    try:
        entries = transvar_config.fns[(refversion, "anno")]
    except KeyError:
        sys.exit(
            "The installed TransVar release provides no annotation databases "
            "for reference version '%s'." % refversion
        )

    # base .transvardb file name -> flag, for the requested sources only
    bases = {}
    for config_key, filename, _url in entries:
        flag = FLAG_BY_CONFIG_KEY.get(config_key)
        if flag is not None and flag in wanted_flags:
            bases[filename] = flag

    missing = sorted(set(wanted_flags) - set(bases.values()))
    if missing:
        sys.exit(
            "Reference version '%s' has no %s annotation database in this "
            "TransVar release." % (refversion, ", ".join(missing))
        )

    downloads = []
    for _config_key, filename, url in entries:
        for base, flag in bases.items():
            if filename == base or filename.startswith(base + "."):
                target = flag + ".transvardb" + filename[len(base):]
                downloads.append((target, url))
                break
    return downloads


def download(url, target_path, test=False):
    """Download url to target_path, atomically. Returns False if absent."""
    if test:
        # Exercise the file naming without fetching hundreds of megabytes
        # from an external server.
        open(target_path, "w").close()
        return True

    # The per-database ID mapping indices only exist for some sources and
    # TransVar works without them, so their absence is not an error.
    optional = target_path.endswith(".idmap_idx")
    last_error = None

    for attempt in range(1, ATTEMPTS + 1):
        tmp_fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(target_path))
        os.close(tmp_fd)
        try:
            with urlopen(url, timeout=TIMEOUT) as response, open(tmp_path, "wb") as out:
                shutil.copyfileobj(response, out)
        except HTTPError as e:
            os.remove(tmp_path)
            last_error = e
            if e.code < 500:
                break  # the file is not there, retrying will not help
        except Exception as e:
            os.remove(tmp_path)
            last_error = e
        else:
            # mkstemp creates the file 0600, but the databases end up in a
            # shared data directory every job user has to be able to read.
            os.chmod(tmp_path, 0o644)
            os.rename(tmp_path, target_path)
            return True
        if attempt < ATTEMPTS:
            time.sleep(2 ** attempt)

    if optional:
        return False
    sys.exit(
        "Failed to download %s after %d attempt(s): %s\n%s"
        % (url, attempt, last_error, PROXY_HINT)
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out_file")
    parser.add_argument("--refversion", required=True)
    parser.add_argument(
        "--sources",
        required=True,
        help="comma separated list of annotation source flags",
    )
    parser.add_argument("--name", default="")
    parser.add_argument(
        "--test",
        action="store_true",
        help="create empty placeholder files instead of downloading",
    )
    args = parser.parse_args()

    wanted_flags = [s for s in args.sources.split(",") if s]
    if not wanted_flags:
        sys.exit("Select at least one transcript annotation source.")
    unknown = sorted(set(wanted_flags) - set(FLAG_BY_CONFIG_KEY.values()))
    if unknown:
        sys.exit("Unsupported annotation source(s): %s" % ", ".join(unknown))

    with open(args.out_file) as fh:
        params = json.load(fh)
    target_directory = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_directory, exist_ok=True)

    downloaded_flags = set()
    for basename, url in files_to_download(args.refversion, wanted_flags):
        if download(url, os.path.join(target_directory, basename), args.test):
            downloaded_flags.add(basename.split(".", 1)[0])

    sources = sorted(downloaded_flags)
    value = "%s_%s" % (args.refversion, "-".join(sources))
    name = args.name or "%s (%s)" % (
        args.refversion,
        ", ".join(SOURCE_LABELS[flag] for flag in sources),
    )

    data_manager_dict = {
        "data_tables": {
            "transvar_annotations": [
                {
                    "value": value,
                    "dbkey": args.refversion,
                    "refversion": args.refversion,
                    "name": name,
                    "sources": ",".join(sources),
                    "path": target_directory,
                }
            ]
        }
    }
    with open(args.out_file, "w") as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
