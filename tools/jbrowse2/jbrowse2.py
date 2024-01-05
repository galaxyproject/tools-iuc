#!/usr/bin/env python
# change to accumulating all configuration for config.json based on the default from the clone
import argparse
import datetime
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("jbrowse")
TODAY = datetime.datetime.now().strftime("%Y-%m-%d")
GALAXY_INFRASTRUCTURE_URL = None
mapped_chars = {
    ">": "__gt__",
    "<": "__lt__",
    "'": "__sq__",
    '"': "__dq__",
    "[": "__ob__",
    "]": "__cb__",
    "{": "__oc__",
    "}": "__cc__",
    "@": "__at__",
    "#": "__pd__",
    "": "__cn__",
}


def etree_to_dict(t):
    if t is None:
        return {}

    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d


INSTALLED_TO = os.path.dirname(os.path.realpath(__file__))


def metadata_from_node(node):
    metadata = {}
    try:
        if len(node.findall("dataset")) != 1:
            # exit early
            return metadata
    except Exception:
        return {}

    for (key, value) in node.findall("dataset")[0].attrib.items():
        metadata["dataset_%s" % key] = value

    for (key, value) in node.findall("history")[0].attrib.items():
        metadata["history_%s" % key] = value

    for (key, value) in node.findall("metadata")[0].attrib.items():
        metadata["metadata_%s" % key] = value

    for (key, value) in node.findall("tool")[0].attrib.items():
        metadata["tool_%s" % key] = value

    # Additional Mappings applied:
    metadata[
        "dataset_edam_format"
    ] = '<a target="_blank" href="http://edamontology.org/{0}">{1}</a>'.format(
        metadata["dataset_edam_format"], metadata["dataset_file_ext"]
    )
    metadata["history_user_email"] = '<a href="mailto:{0}">{0}</a>'.format(
        metadata["history_user_email"]
    )
    metadata["hist_name"] = metadata["history_display_name"]
    metadata[
        "history_display_name"
    ] = '<a target="_blank" href="{galaxy}/history/view/{encoded_hist_id}">{hist_name}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_hist_id=metadata["history_id"],
        hist_name=metadata["history_display_name"],
    )
    metadata[
        "tool_tool"
    ] = '<a target="_blank" href="{galaxy}/datasets/{encoded_id}/show_params">{tool_id}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_id=metadata["dataset_id"],
        tool_id=metadata["tool_tool_id"],
        # tool_version=metadata['tool_tool_version'],
    )
    return metadata


class JbrowseConnector(object):
    def __init__(self, jbrowse, outdir, genomes):
        self.debug = False
        self.usejson = True
        self.giURL = GALAXY_INFRASTRUCTURE_URL
        self.jbrowse = jbrowse
        self.outdir = outdir
        os.makedirs(self.outdir, exist_ok=True)
        self.genome_paths = genomes
        self.trackIdlist = []
        self.tracksToAdd = []
        self.config_json = {}
        self.config_json_file = os.path.join(outdir, "config.json")
        self.clone_jbrowse(self.jbrowse, self.outdir)

    def subprocess_check_call(self, command, output=None):
        if output:
            if self.debug:
                log.debug("cd %s && %s >  %s", self.outdir, " ".join(command), output)
            subprocess.check_call(command, cwd=self.outdir, stdout=output)
        else:
            log.debug("cd %s && %s", self.outdir, " ".join(command))
            subprocess.check_call(command, cwd=self.outdir)

    def subprocess_popen(self, command):
        if self.debug:
            log.debug(command)
        p = subprocess.Popen(
            command,
            cwd=self.outdir,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, err = p.communicate()
        retcode = p.returncode
        if retcode != 0:
            log.error(command)
            log.error(output)
            log.error(err)
            raise RuntimeError("Command failed with exit code %s" % (retcode))

    def subprocess_check_output(self, command):
        if self.debug:
            log.debug(" ".join(command))
        return subprocess.check_output(command, cwd=self.outdir)

    def _jbrowse_bin(self, command):
        return os.path.join(self.jbrowse, "bin", command)

    def symlink_or_copy(self, src, dest):
        if "GALAXY_JBROWSE_SYMLINKS" in os.environ and bool(
            os.environ["GALAXY_JBROWSE_SYMLINKS"]
        ):
            cmd = ["ln", "-s", src, dest]
        else:
            cmd = ["cp", src, dest]

        return self.subprocess_check_call(cmd)

    def process_genomes(self):
        assemblies = []
        for i, genome_node in enumerate(self.genome_paths):
            if self.debug:
                log.info("genome_node=%s" % str(genome_node))
            genome_name = genome_node["meta"]["dataset_dname"].strip().split()[0]
            fapath = genome_node["path"]
            faname = genome_name + ".fa.gz"
            fadest = os.path.join(self.outdir, faname)
            # fadest = os.path.realpath(os.path.join(self.outdir, faname))
            cmd = "bgzip -i -c %s -I %s.gzi > %s && samtools faidx %s" % (fapath, fadest, fadest, fadest)
            if self.debug:
                log.info("### cmd = %s" % cmd)
            self.subprocess_popen(cmd)
            adapter = {
                "type": "BgzipFastaAdapter",
                "fastaLocation": {
                    "uri": faname,
                },
                "faiLocation": {
                    "uri": faname + ".fai",
                },
                "gziLocation": {
                    "uri": faname + ".gzi",
                },
            }
            trackDict = {
                "name": genome_name,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": genome_name,
                    "adapter": adapter,
                },
                "rendering": {"type": "DivSequenceRenderer"},
            }
            assemblies.append(trackDict)
        self.genome_name = genome_name
        if self.usejson:
            self.config_json["assemblies"] = assemblies
        else:
            cmd = [
                "jbrowse",
                "add-assembly",
                faname,
                "-t",
                "bgzipFasta",
                "-n",
                genome_name,
                "--load",
                "inPlace",
                "--faiLocation",
                faname + ".fai",
                "--gziLocation",
                faname + ".gzi",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def add_default_view(self):
        cmd = [
            "jbrowse",
            "set-default-session",
            "-s",
            self.config_json_file,
            "-t",
            ",".join(self.trackIdlist),
            "-n",
            "JBrowse2 in Galaxy",
            "--target",
            self.config_json_file,
            "-v",
            " LinearGenomeView",
        ]
        if self.debug:
            log.info("### calling set-default-session with cmd=%s" % "  ".join(cmd))
        self.subprocess_check_call(cmd)

    def write_config(self):
        with open(self.config_json_file, "w") as fp:
            json.dump(self.config_json, fp)

    def add_hic(self, data, trackData):
        """
        HiC adapter.
        https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md
        for testing locally, these work:
        HiC data is from https://s3.amazonaws.com/igv.broadinstitute.org/data/hic/intra_nofrag_30.hic
        using hg19 reference track as a
        'BgzipFastaAdapter'
            fastaLocation:
            uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz',
            faiLocation:
            uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz.fai',
            gziLocation:
            uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/fasta/GRCh38.fa.gz.gzi',
        Cool will not be likely to be a good fit - see discussion at https://github.com/GMOD/jbrowse-components/issues/2438
        """
        log.info("#### trackData=%s" % trackData)
        tId = trackData["label"]
        dsId = trackData["metadata"]["dataset_id"]
        url = "%s/api/datasets/%s/display?to_ext=hic " % (
            self.giURL,
            dsId,
        )
        hname = trackData["name"]
        dest = os.path.join(self.outdir, hname)
        url = hname
        cmd = ["cp", data, dest]
        self.subprocess_check_call(cmd)
        floc = {
            "uri": hname,
        }
        trackDict = {
            "type": "HicTrack",
            "trackId": tId,
            "name": hname,
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "HicAdapter",
                "hicLocation": floc,
            },
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "HicTrack",
                "-a",
                self.genome_name,
                "-n",
                hname,
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def add_maf(self, data, trackData):
        """
        from https://github.com/cmdcolin/maf2bed
        Note: Both formats start with a MAF as input, and note that your MAF file should contain the species name and chromosome name
        e.g. hg38.chr1 in the sequence identifiers.
        need the reference id - eg hg18, for maf2bed.pl as the first parameter
        """
        mafPlugin = {
            "plugins": [
                {
                    "name": "MafViewer",
                    "url": "https://unpkg.com/jbrowse-plugin-mafviewer/dist/jbrowse-plugin-mafviewer.umd.production.min.js",
                }
            ]
        }
        tId = trackData["label"]
        fname = "%s.bed" % tId
        dest = "%s/%s" % (self.outdir, fname)
        # self.symlink_or_copy(data, dest)
        # Process MAF to bed-like. Need build to munge chromosomes
        gname = self.genome_name
        cmd = [
            "bash",
            os.path.join(INSTALLED_TO, "convertMAF.sh"),
            data,
            gname,
            INSTALLED_TO,
            dest,
        ]
        self.subprocess_check_call(cmd)
        if True or self.debug:
            log.info("### convertMAF.sh called as %s" % " ".join(cmd))
        # Construct samples list
        # We could get this from galaxy metadata, not sure how easily.
        ps = subprocess.Popen(["grep", "^s [^ ]*", "-o", data], stdout=subprocess.PIPE)
        output = subprocess.check_output(("sort", "-u"), stdin=ps.stdout)
        ps.wait()
        outp = output.decode("ascii")
        soutp = outp.split("\n")
        samp = [x.split("s ")[1] for x in soutp if x.startswith("s ")]
        samples = [x.split(".")[0] for x in samp]
        if self.debug:
            log.info("### got samples = %s " % (samples))
        trackDict = {
            "type": "MafTrack",
            "trackId": tId,
            "name": trackData["name"],
            "adapter": {
                "type": "MafTabixAdapter",
                "samples": samples,
                "bedGzLocation": {
                    "uri": fname + ".sorted.bed.gz",
                },
                "index": {
                    "location": {
                        "uri": fname + ".sorted.bed.gz.tbi",
                    },
                },
            },
            "assemblyNames": [self.genome_name],
        }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)
        if self.config_json.get("plugins", None):
            self.config_json["plugins"].append(mafPlugin[0])
        else:
            self.config_json.update(mafPlugin)

    def _blastxml_to_gff3(self, xml, min_gap=10):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = [
            "python",
            os.path.join(INSTALLED_TO, "blastxml_to_gapped_gff3.py"),
            "--trim",
            "--trim_end",
            "--include_seq",
            "--min_gap",
            str(min_gap),
            xml,
        ]
        subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_unrebased)
        gff3_unrebased.close()
        return gff3_unrebased.name

    def add_blastxml(self, data, trackData, blastOpts, **kwargs):
        gff3 = self._blastxml_to_gff3(data, min_gap=blastOpts["min_gap"])

        if "parent" in blastOpts and blastOpts["parent"] != "None":
            gff3_rebased = tempfile.NamedTemporaryFile(delete=False)
            cmd = ["python", os.path.join(INSTALLED_TO, "gff3_rebase.py")]
            if blastOpts.get("protein", "false") == "true":
                cmd.append("--protein2dna")
            cmd.extend([os.path.realpath(blastOpts["parent"]), gff3])
            subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_rebased)
            gff3_rebased.close()

            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)
        url = "%s.gff3" % trackData["label"]
        dest = "%s/%s" % (self.outdir, url)
        self._sort_gff(gff3, dest)
        url = url + ".gz"
        tId = trackData["label"]
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                    "uri": url,
                },
                "index": {
                    "location": {
                        "uri": url + ".tbi",
                    }
                },
            },
            "displays": [
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearBasicDisplay" % tId,
                },
                {"type": "LinearArcDisplay", "displayId": "%s-LinearArcDisplay" % tId},
            ],
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "FeatureTrack",
                "-a",
                self.genome_name,
                "--indexFile",
                url + ".tbi",
                "-n",
                trackData["name"],
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)
        os.unlink(gff3)

    def add_bigwig(self, data, trackData):
        url = "%s.bw" % trackData["name"]
        dest = os.path.join(self.outdir, url)
        cmd = ["cp", data, dest]
        self.subprocess_check_call(cmd)
        bwloc = {"uri": url}
        tId = trackData["label"]
        trackDict = {
            "type": "QuantitativeTrack",
            "trackId": tId,
            "name": url,
            "assemblyNames": [
                self.genome_name,
            ],
            "adapter": {
                "type": "BigWigAdapter",
                "bigWigLocation": bwloc,
            },
            "displays": [
                {
                    "type": "LinearWiggleDisplay",
                    "displayId": "%s-LinearWiggleDisplay" % tId,
                }
            ],
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "QuantitativeTrack",
                "-a",
                self.genome_name,
                "-n",
                trackData["name"],
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def add_bam(self, data, trackData, bamOpts, bam_index=None, **kwargs):
        tId = trackData["label"]
        fname = "%s.bam" % trackData["label"]
        dest = "%s/%s" % (self.outdir, fname)
        url = fname
        self.subprocess_check_call(["cp", data, dest])
        log.info("### copied %s to %s" % (data, dest))
        bloc = {"uri": url}
        if bam_index is not None and os.path.exists(os.path.realpath(bam_index)):
            # bai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
            self.subprocess_check_call(
                ["cp", os.path.realpath(bam_index), dest + ".bai"]
            )
        else:
            # Can happen in exotic condition
            # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
            #      => no index generated by galaxy, but there might be one next to the symlink target
            #      this trick allows to skip the bam sorting made by galaxy if already done outside
            if os.path.exists(os.path.realpath(data) + ".bai"):
                self.symlink_or_copy(os.path.realpath(data) + ".bai", dest + ".bai")
            else:
                log.warn("Could not find a bam index (.bai file) for %s", data)
        trackDict = {
            "type": "AlignmentsTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "BamAdapter",
                "bamLocation": bloc,
                "index": {
                    "location": {
                        "uri": fname + ".bai",
                    }
                },
            },
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                fname,
                "-t",
                "AlignmentsTrack",
                "-l",
                "inPlace",
                "-a",
                self.genome_name,
                "--indexFile",
                fname + ".bai",
                "-n",
                trackData["name"],
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def add_vcf(self, data, trackData):
        tId = trackData["label"]
        url = "%s/api/datasets/%s/display" % (
            self.giURL,
            trackData["metadata"]["dataset_id"],
        )
        url = "%s.vcf.gz" % tId
        dest = "%s/%s" % (self.outdir, url)
        cmd = "bgzip -c %s  > %s" % (data, dest)
        self.subprocess_popen(cmd)
        cmd = ["tabix", "-p", "vcf", dest]
        self.subprocess_check_call(cmd)
        trackDict = {
            "type": "VariantTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "VcfTabixAdapter",
                "vcfGzLocation": {
                    "uri": url,
                },
                "index": {
                    "location": {
                        "uri": url + ".tbi",
                    }
                },
            },
            "displays": [
                {
                    "type": "LinearVariantDisplay",
                    "displayId": "%s-LinearVariantDisplay" % tId,
                },
                {
                    "type": "ChordVariantDisplay",
                    "displayId": "%s-ChordVariantDisplay" % tId,
                },
                {
                    "type": "LinearPairedArcDisplay",
                    "displayId": "%s-LinearPairedArcDisplay" % tId,
                },
            ],
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "VariantTrack",
                "-a",
                self.genome_name,
                "--indexFile",
                url + ".tbi",
                "-n",
                trackData["name"],
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest + ".gz"):
            cmd = "jbrowse sort-gff %s | bgzip -c > %s.gz" % (
                data,
                dest,
            )  # "gff3sort.pl --precise '%s' | grep -v \"^$\" > '%s'"
            self.subprocess_popen(cmd)
            self.subprocess_check_call(["tabix", "-f", "-p", "gff", dest + ".gz"])

    def _sort_bed(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "sort -k1,1 -k2,2n %s | bgzip -c > %s" % (data, dest)
            self.subprocess_popen(cmd)
            cmd = ["tabix", "-f", "-p", "bed", dest]
            self.subprocess_check_call(cmd)

    def add_gff(self, data, ext, trackData):
        url = "%s.%s" % (trackData["label"], ext)
        dest = "%s/%s" % (self.outdir, url)
        self._sort_gff(data, dest)
        url = url + ".gz"
        tId = trackData["label"]
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                    "uri": url,
                },
                "index": {
                    "location": {
                        "uri": url + ".tbi",
                    }
                },
            },
            "displays": [
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearBasicDisplay" % tId,
                },
                {"type": "LinearArcDisplay", "displayId": "%s-LinearArcDisplay" % tId},
            ],
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "FeatureTrack",
                "-a",
                self.genome_name,
                "-n",
                trackData["name"],
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def add_bed(self, data, ext, trackData):
        url = "%s.%s" % (trackData["label"], ext)
        dest = "%s/%s.gz" % (self.outdir, url)
        self._sort_bed(data, dest)
        tId = trackData["label"]
        url = url + ".gz"
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "BedTabixAdapter",
                "bedGzLocation": {
                    "uri": url,
                },
                "index": {
                    "location": {
                        "uri": url + ".tbi",
                    }
                },
            },
            "displays": [
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearBasicDisplay" % tId,
                },
                {"type": "LinearArcDisplay", "displayId": "%s-LinearArcDisplay" % tId},
            ],
        }
        if self.usejson:
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            cmd = [
                "jbrowse",
                "add-track",
                url,
                "-t",
                "FeatureTrack",
                "-a",
                self.genome_name,
                "--indexFile",
                url + ".tbi",
                "-n",
                trackData["name"],
                "--load",
                "inPlace",
                "--target",
                self.outdir,
            ]
            self.subprocess_check_call(cmd)

    def process_annotations(self, track):
        category = track["category"].replace("__pd__date__pd__", TODAY)
        for i, (
            dataset_path,
            dataset_ext,
            track_human_label,
            extra_metadata,
        ) in enumerate(track["trackfiles"]):
            # Unsanitize labels (element_identifiers are always sanitized by Galaxy)
            for key, value in mapped_chars.items():
                track_human_label = track_human_label.replace(value, key)
            outputTrackConfig = {
                "category": category,
            }
            if self.debug:
                log.info(
                    "Processing category = %s, track_human_label = %s",
                    category,
                    track_human_label,
                )
            # We add extra data to hash for the case of REST + SPARQL.
            if (
                "conf" in track
                and "options" in track["conf"]
                and "url" in track["conf"]["options"]
            ):
                rest_url = track["conf"]["options"]["url"]
            else:
                rest_url = ""

            # I chose to use track['category'] instead of 'category' here. This
            # is intentional. This way re-running the tool on a different date
            # will not generate different hashes and make comparison of outputs
            # much simpler.
            hashData = [
                str(dataset_path),
                track_human_label,
                track["category"],
                rest_url,
            ]
            hashData = "|".join(hashData).encode("utf-8")
            outputTrackConfig["label"] = hashlib.md5(hashData).hexdigest() + "_%s" % i
            outputTrackConfig["metadata"] = extra_metadata
            outputTrackConfig["name"] = track_human_label

            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                )
            elif dataset_ext in ("hic",):
                self.add_hic(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext in ("bed",):
                self.add_bed(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                )
            elif dataset_ext in ("maf",):
                self.add_maf(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext == "bigwig":
                self.add_bigwig(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext == "bam":
                real_indexes = track["conf"]["options"]["pileup"]["bam_indices"][
                    "bam_index"
                ]
                if not isinstance(real_indexes, list):
                    # <bam_indices>
                    #  <bam_index>/path/to/a.bam.bai</bam_index>
                    # </bam_indices>
                    #
                    # The above will result in the 'bam_index' key containing a
                    # string. If there are two or more indices, the container
                    # becomes a list. Fun!
                    real_indexes = [real_indexes]

                self.add_bam(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["pileup"],
                    bam_index=real_indexes[i],
                )
            elif dataset_ext == "blastxml":
                self.add_blastxml(
                    dataset_path, outputTrackConfig, track["conf"]["options"]["blast"]
                )
            elif dataset_ext == "vcf":
                self.add_vcf(dataset_path, outputTrackConfig)
            else:
                log.warn("Do not know how to handle %s", dataset_ext)

    def clone_jbrowse(self, jbrowse_dir, destination):
        """Clone a JBrowse directory into a destination directory."""
        cmd = ["jbrowse", "create", "-f", os.path.realpath(self.outdir)]
        self.subprocess_check_call(cmd)
        for fn in [
            "asset-manifest.json",
            "favicon.ico",
            "robots.txt",
            "umd_plugin.js",
            "version.txt",
            "test_data",
        ]:
            cmd = ["rm", "-rf", os.path.join(self.outdir, fn)]
            self.subprocess_check_call(cmd)
        cmd = ["cp", os.path.join(INSTALLED_TO, "servejb2.py"), self.outdir]
        self.subprocess_check_call(cmd)

    def clone_jbrowse2(self, jbrowse_dir, destination):
        """Clone a JBrowse directory into a destination directory."""
        cmd = ["cp", "-rv", jbrowse_dir + "/*", self.outdir]
        self.subprocess_check_call(cmd)
        cmd = ["cp", os.path.join(INSTALLED_TO, "servejb2.py"), self.outdir]
        self.subprocess_check_call(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("xml", type=argparse.FileType("r"), help="Track Configuration")

    parser.add_argument("--jbrowse", help="Folder containing a jbrowse release")
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument("--version", "-V", action="version", version="%(prog)s 0.8.0")
    args = parser.parse_args()

    tree = ET.parse(args.xml.name)
    root = tree.getroot()

    # This should be done ASAP
    GALAXY_INFRASTRUCTURE_URL = root.find("metadata/galaxyUrl").text
    # Sometimes this comes as `localhost` without a protocol
    if not GALAXY_INFRASTRUCTURE_URL.startswith("http"):
        # so we'll prepend `http://` and hope for the best. Requests *should*
        # be GET and not POST so it should redirect OK
        GALAXY_INFRASTRUCTURE_URL = "http://" + GALAXY_INFRASTRUCTURE_URL
    jb = args.jbrowse
    jb1, one = os.path.split(jb)
    jb1 += "/opt/jbrowse2"  # /../opt/jbrowse
    jb2, two = os.path.split(jb1)
    jb2 += "/opt/jbrowse2"  # /../../opt/jbrowse for container
    if os.path.exists(jb1) and "manifest.json" in os.listdir(jb1):
        jbdir = jb1
    elif os.path.exists(jb2) and "manifest.json" in os.listdir(jb2):
        jbdir = jb2
    else:
        log.error(
            "unable to find the jbrowse2 directory for cloning jb1= %s, jb2 = %s - args.jbrowse = %s"
            % (jb1, jb2, args.jbrowse)
        )
        sys.exit(10)
    jc = JbrowseConnector(
        jbrowse=jbdir,
        outdir=args.outdir,
        genomes=[
            {
                "path": os.path.realpath(x.attrib["path"]),
                "meta": metadata_from_node(x.find("metadata")),
            }
            for x in root.findall("metadata/genomes/genome")
        ],
    )
    jc.process_genomes()

    for track in root.findall("tracks/track"):
        track_conf = {}
        track_conf["trackfiles"] = []

        is_multi_bigwig = False
        try:
            if track.find("options/wiggle/multibigwig") and (
                track.find("options/wiggle/multibigwig").text == "True"
            ):
                is_multi_bigwig = True
                multi_bigwig_paths = []
        except KeyError:
            pass

        trackfiles = track.findall("files/trackFile")
        if trackfiles:
            for x in track.findall("files/trackFile"):
                if is_multi_bigwig:
                    multi_bigwig_paths.append(
                        (x.attrib["label"], os.path.realpath(x.attrib["path"]))
                    )
                else:
                    if trackfiles:
                        metadata = metadata_from_node(x.find("metadata"))
                        track_conf["dataset_id"] = metadata["dataset_id"]
                        track_conf["trackfiles"].append(
                            (
                                os.path.realpath(x.attrib["path"]),
                                x.attrib["ext"],
                                x.attrib["label"],
                                metadata,
                            )
                        )
        else:
            # For tracks without files (rest, sparql)
            track_conf["trackfiles"].append(
                (
                    "",  # N/A, no path for rest or sparql
                    track.attrib["format"],
                    track.find("options/label").text,
                    {},
                )
            )

        if is_multi_bigwig:
            metadata = metadata_from_node(x.find("metadata"))

            track_conf["trackfiles"].append(
                (
                    multi_bigwig_paths,  # Passing an array of paths to represent as one track
                    "bigwig_multiple",
                    "MultiBigWig",  # Giving an hardcoded name for now
                    {},  # No metadata for multiple bigwig
                )
            )

        track_conf["category"] = track.attrib["cat"]
        track_conf["format"] = track.attrib["format"]
        try:
            # Only pertains to gff3 + blastxml. TODO?
            track_conf["style"] = {t.tag: t.text for t in track.find("options/style")}
        except TypeError:
            track_conf["style"] = {}
            pass
        track_conf["conf"] = etree_to_dict(track.find("options"))
        jc.process_annotations(track_conf)
        print("## processed", str(track_conf), "trackIdlist", jc.trackIdlist)
    print(
        "###done processing, trackIdlist=",
        jc.trackIdlist,
        "config=",
        str(jc.config_json),
    )
    jc.config_json["tracks"] = jc.tracksToAdd
    if jc.usejson:
        jc.write_config()
    jc.add_default_view()
