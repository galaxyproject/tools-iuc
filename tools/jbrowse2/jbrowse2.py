#!/usr/bin/env python
import argparse
import datetime
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("jbrowse")
TODAY = datetime.datetime.now().strftime("%Y-%m-%d")
SELF_LOCATION = os.path.dirname(os.path.realpath(__file__))
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

    for key, value in node.findall("dataset")[0].attrib.items():
        metadata[f"dataset_{key}"] = value

    for key, value in node.findall("history")[0].attrib.items():
        metadata[f"history_{key}"] = value

    for key, value in node.findall("metadata")[0].attrib.items():
        metadata[f"metadata_{key}"] = value

    for key, value in node.findall("tool")[0].attrib.items():
        metadata[f"tool_{key}"] = value

    # Additional Mappings applied:
    metadata[
        "dataset_edam_format"
    ] = '<a target="_blank" href="http://edamontology.org/{0}">{1}</a>'.format(
        metadata["dataset_edam_format"], metadata["dataset_file_ext"]
    )
    metadata["history_user_email"] = '<a href="mailto:{0}">{0}</a>'.format(
        metadata["history_user_email"]
    )
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
    def __init__(self, jbrowse, outdir, update):
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.update = update

        self.tracksToIndex = {}

        # This is the id of the current assembly
        self.assembly_ids = {}

        self.default_views = {}

        # If upgrading, look at the existing data
        self.check_existing(self.outdir)

        self.clone_jbrowse(self.jbrowse, self.outdir)

    def get_cwd(self, cwd):
        if cwd:
            return self.outdir
        else:
            return subprocess.check_output(['pwd']).decode('utf-8').strip()
            # return None

    def subprocess_check_call(self, command, output=None, cwd=True):
        if output:
            log.debug(f"cd {self.get_cwd(cwd)} && {' '.join(command)} > {output}")
            subprocess.check_call(command, cwd=self.get_cwd(cwd), stdout=output)
        else:
            log.debug(f"cd {self.get_cwd(cwd)} && {' '.join(command)}")
            subprocess.check_call(command, cwd=self.get_cwd(cwd))

    def subprocess_popen(self, command, cwd=True):
        log.debug(f"cd {self.get_cwd(cwd)} && {command}")
        p = subprocess.Popen(
            command,
            cwd=self.get_cwd(cwd),
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, err = p.communicate()
        retcode = p.returncode
        if retcode != 0:
            log.error(f"cd {self.get_cwd(cwd)} && {command}")
            log.error(output)
            log.error(err)
            raise RuntimeError(f"Command failed with exit code {retcode}")

    def subprocess_check_output(self, command, cwd=True):
        log.debug(f"cd {self.get_cwd(cwd)} && {' '.join(command)}")
        return subprocess.check_output(command, cwd=self.get_cwd(cwd))

    def symlink_or_copy(self, src, dest):
        if "GALAXY_JBROWSE_SYMLINKS" in os.environ and bool(
            os.environ["GALAXY_JBROWSE_SYMLINKS"]
        ):
            cmd = ["ln", "-s", src, dest]
        else:
            cmd = ["cp", src, dest]

        return self.subprocess_check_call(cmd)

    def _prepare_track_style(self, xml_conf):
        style_data = {
            "type": "LinearBasicDisplay",  # TODO choose a better default?
        }

        if "display" in xml_conf["style"]:
            style_data["type"] = xml_conf["style"]["display"]

        style_data["displayId"] = f"{xml_conf['label']}_{style_data['type']}"

        style_data.update(self._prepare_renderer_config(style_data["type"], xml_conf["style"]))

        return {"displays": [style_data]}

    def _prepare_renderer_config(self, display_type, xml_conf):

        style_data = {}

        # if display_type in ("LinearBasicDisplay", "LinearVariantDisplay"):
        # TODO LinearVariantDisplay does not understand these options when written in config.json
        if display_type in ("LinearBasicDisplay"):

            # Doc: https://jbrowse.org/jb2/docs/config/svgfeaturerenderer/
            style_data["renderer"] = {
                "type": "SvgFeatureRenderer",
                "showLabels": xml_conf.get("show_labels", True),
                "showDescriptions": xml_conf.get("show_descriptions", True),
                "labels": {
                    "name": xml_conf.get("labels_name", "jexl:get(feature,'name') || get(feature,'id')"),
                    "description": xml_conf.get("descriptions_name", "jexl:get(feature,'note') || get(feature,'description')")
                },
                "displayMode": xml_conf.get("display_mode", "normal"),
                "maxHeight": xml_conf.get("max_height", 1200),
            }

        elif display_type == "LinearArcDisplay":

            # Doc: https://jbrowse.org/jb2/docs/config/arcrenderer/
            style_data["renderer"] = {
                "type": "ArcRenderer",
                "label": xml_conf.get("labels_name", "jexl:get(feature,'score')"),
                "displayMode": xml_conf.get("display_mode", "arcs"),
            }

        elif display_type == "LinearWiggleDisplay":

            wig_renderer = xml_conf.get("renderer", "xyplot")
            style_data["defaultRendering"] = wig_renderer

            # Doc: https://jbrowse.org/jb2/docs/config/snpcoveragerenderer/

        elif display_type == "LinearSNPCoverageDisplay":

            # Does not work
            # style_data["renderer"] = {
            #     "type": "SNPCoverageRenderer",
            #     "displayCrossHatches": xml_conf.get("display_cross_hatches", True),
            # }

            style_data["scaleType"] = xml_conf.get("scale_type", "linear")
            if "min_score" in xml_conf:
                style_data["minScore"] = xml_conf["min_score"]

            if "max_score" in xml_conf:
                style_data["maxScore"] = xml_conf["max_score"]

            # Doc: https://jbrowse.org/jb2/docs/config/snpcoveragerenderer/

        return style_data

    def _prepare_format_details(self, xml_conf):
        formatDetails = {
        }

        if "feature" in xml_conf["formatdetails"]:
            feat_jexl = xml_conf["formatdetails"]["feature"]
            for key, value in mapped_chars.items():
                feat_jexl = feat_jexl.replace(value, key)
            formatDetails["feature"] = feat_jexl

        if "subfeature" in xml_conf["formatdetails"]:
            sfeat_jexl = xml_conf["formatdetails"]["subfeature"]
            for key, value in mapped_chars.items():
                sfeat_jexl = sfeat_jexl.replace(value, key)
            formatDetails["subfeatures"] = sfeat_jexl

        if "depth" in xml_conf["formatdetails"]:
            formatDetails["depth"] = int(xml_conf["formatdetails"]["depth"])

        return {"formatDetails": formatDetails}

    def check_existing(self, destination):
        existing = os.path.join(destination, "config.json")
        if os.path.exists(existing):
            with open(existing, "r") as existing_conf:
                conf = json.load(existing_conf)
                if "assemblies" in conf:
                    for assembly in conf["assemblies"]:
                        if "name" in assembly:
                            self.assembly_ids[assembly["name"]] = None

    def add_assembly(self, path, label):
        # Find a non-existing filename for the new genome
        # (to avoid colision when upgrading an existing instance)
        rel_seq_path = os.path.join("data", label)
        seq_path = os.path.join(self.outdir, rel_seq_path)
        fn_try = 1
        while (
            os.path.exists(seq_path + ".fasta")
            or os.path.exists(seq_path + ".fasta.gz")
            or os.path.exists(seq_path + ".fasta.gz.fai")
            or os.path.exists(seq_path + ".fasta.gz.gzi")
        ):
            rel_seq_path = os.path.join("data", f"{label}{fn_try}")
            seq_path = os.path.join(self.outdir, rel_seq_path)
            fn_try += 1

        # Check if the assembly already exists from a previous run (--update mode)
        if self.update:

            config_path = os.path.join(self.outdir, "config.json")
            with open(config_path, "r") as config_file:
                config_json = json.load(config_file)

                for asby in config_json['assemblies']:
                    if asby['name'] == label:
                        log.info("Found existing assembly from existing JBrowse2 instance, preserving it")

                        self.assembly_ids[label] = ""

                        # Find default views existing for this assembly
                        if 'defaultSession' in config_json and 'views' in config_json['defaultSession']:
                            for view in config_json['defaultSession']['views']:
                                if 'init' in view and 'assembly' in view['init']:
                                    if view['init']['assembly'] == label:
                                        self.default_views[view['init']['assembly']] = view

                        return label

        # Find a non-existing label for the new genome
        # (to avoid colision when upgrading an existing instance)
        lab_try = 1
        uniq_label = label
        while uniq_label in self.assembly_ids:
            uniq_label = label + str(lab_try)
            lab_try += 1

        # Find a default scaffold to display
        with open(path, "r") as fa_handle:
            fa_header = fa_handle.readline()[1:].strip().split(" ")[0]

        self.assembly_ids[uniq_label] = fa_header

        copied_genome = seq_path + ".fasta"
        shutil.copy(path, copied_genome)

        # Compress with bgzip
        cmd = ["bgzip", copied_genome]
        self.subprocess_check_call(cmd)

        # FAI Index
        cmd = ["samtools", "faidx", copied_genome + ".gz"]
        self.subprocess_check_call(cmd)

        cmd_jb = [
            "jbrowse",
            "add-assembly",
            "--load",
            "inPlace",
            "--name",
            uniq_label,
            "--type",
            "bgzipFasta",
            "--out",
            self.outdir,
            "--skipCheck",
            rel_seq_path + ".fasta.gz",
        ]

        self.subprocess_check_call(cmd_jb)

        return uniq_label

    def text_index(self):

        for ass in self.tracksToIndex:
            tracks = self.tracksToIndex[ass]
            args = [
                "jbrowse",
                "text-index",
                "--target",
                self.outdir,
                "--assemblies",
                ass,
            ]

            tracks = ",".join(tracks)
            if tracks:
                args += ["--tracks", tracks]

                log.info(f"-----> Running text-index on assembly {ass} and tracks {tracks}")

                # Only run index if we want to index at least one
                # If --tracks is not specified, it will index everything
                self.subprocess_check_call(args)

    def add_bigwig(self, parent, data, trackData, wiggleOpts, **kwargs):
        rel_dest = os.path.join("data", trackData["label"] + ".bw")
        dest = os.path.join(self.outdir, rel_dest)
        self.symlink_or_copy(os.path.realpath(data), dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
        )

    # Anything ending in "am" (Bam or Cram)
    def add_xam(self, parent, data, trackData, xamOpts, index=None, ext="bam", **kwargs):
        index_ext = "bai"
        if ext == "cram":
            index_ext = "crai"

        rel_dest = os.path.join("data", trackData["label"] + f".{ext}")
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        if index is not None and os.path.exists(os.path.realpath(index)):
            # xai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
            self.subprocess_check_call(
                ["cp", os.path.realpath(index), dest + f".{index_ext}"]
            )
        else:
            # Can happen in exotic condition
            # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
            #      => no index generated by galaxy, but there might be one next to the symlink target
            #      this trick allows to skip the bam sorting made by galaxy if already done outside
            if os.path.exists(os.path.realpath(data) + f".{index_ext}"):
                self.symlink_or_copy(
                    os.path.realpath(data) + f".{index_ext}", dest + f".{index_ext}"
                )
            else:
                log.warn(
                    f"Could not find a bam index (.{index_ext} file) for {data}"
                )

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
        )

    def add_vcf(self, parent, data, trackData, vcfOpts={}, zipped=False, **kwargs):
        if zipped:
            rel_dest = os.path.join("data", trackData["label"] + ".vcf.gz")
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), dest)
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".vcf")
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), dest)

            cmd = ["bgzip", dest]
            self.subprocess_check_call(cmd)
            cmd = ["tabix", dest + ".gz"]
            self.subprocess_check_call(cmd)

            rel_dest = os.path.join("data", trackData["label"] + ".vcf.gz")

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
        )

    def add_gff(self, parent, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join("data", trackData["label"] + ".gff")
        dest = os.path.join(self.outdir, rel_dest)

        self._sort_gff(data, dest)

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        if gffOpts.get('index', 'false') in ("yes", "true", "True"):
            if parent['uniq_id'] not in self.tracksToIndex:
                self.tracksToIndex[parent['uniq_id']] = []
            self.tracksToIndex[parent['uniq_id']].append(trackData["label"])

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest + ".gz",
            parent,
            config=style_json,
        )

    def add_bed(self, parent, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join("data", trackData["label"] + ".bed")
        dest = os.path.join(self.outdir, rel_dest)

        self._sort_bed(data, dest)

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        if gffOpts.get('index', 'false') in ("yes", "true", "True"):
            if parent['uniq_id'] not in self.tracksToIndex:
                self.tracksToIndex[parent['uniq_id']] = []
            self.tracksToIndex[parent['uniq_id']].append(trackData["label"])

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest + ".gz",
            parent,
            config=style_json,
        )

    def add_paf(self, parent, data, trackData, pafOpts, parentgenome, **kwargs):
        # TODO: how to get both parent genomes.

        # print(trackData)
        rel_dest = os.path.join("data", trackData["label"] + ".paf")
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        # TODO: this was disabled because it was adding spurious assemblies,
        # when one of that name already exists.
        # added_assembly = self.add_assembly(
        #     pafOpts["genome"], pafOpts["genome_label"], default=False
        # )

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            f"{parentgenome} v {trackData['key']}",
            trackData["category"],
            rel_dest,
            parent,
            assemblies=[trackData['key'], parentgenome],
            config=style_json,
            trackType="SyntenyTrack",
        )

    def add_hic(self, parent, data, trackData, hicOpts, **kwargs):
        rel_dest = os.path.join("data", trackData["label"] + ".hic")
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
        )

    def add_sparql(self, parent, url, query, query_refnames, trackData):
        json_track_data = {
            "type": "FeatureTrack",
            "trackId": id,
            "name": trackData["label"],
            "adapter": {
                "type": "SPARQLAdapter",
                "endpoint": {"uri": url, "locationType": "UriLocation"},
                "queryTemplate": query,
            },
            "category": [trackData["category"]],
            "assemblyNames": [parent['uniq_id']],
        }

        if query_refnames:
            json_track_data["adapter"]["refNamesQueryTemplate"]: query_refnames

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json_track_data,
            ]
        )

        # TODO Doesn't work as of 1.6.4, might work in the future
        # self.subprocess_check_call([
        #     'jbrowse', 'add-track',
        #     '--trackType', 'sparql',
        #     '--name', trackData["label"],
        #     '--category', trackData['category'],
        #     '--target', os.path.join(self.outdir, 'data'),
        #     '--trackId', id,
        #     '--config', '{"queryTemplate": "%s"}' % query,
        #     url])

    def _add_track(self, id, label, category, path, assembly, config=None, trackType=None, load_action="inPlace"):
        """
        Adds a track to config.json using Jbrowse add-track cli

        By default, using `--load inPlace`: the file is supposed to be already placed at the `path` relative to
        the outdir, `jbrowse add-track` will not touch it and trust us that the file is there and ready to use.

        With `load_action` parameter, you can ask `jbrowse add-track` to copy/move/symlink the file for you.
        Not done by default because we often need more control on file copying/symlink for specific cases (indexes, symlinks of symlinks, ...)
        """

        cmd = [
            "jbrowse",
            "add-track",
            "--load",
            load_action,
            "--name",
            label,
            "--category",
            category,
            "--target",
            self.outdir,
            "--trackId",
            id,
            "--assemblyNames",
            assembly['uniq_id'],
        ]

        if config:
            cmd.append("--config")
            cmd.append(json.dumps(config))

        if trackType:
            cmd.append("--trackType")
            cmd.append(trackType)

        cmd.append(path)

        self.subprocess_check_call(cmd)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            # TODO: replace with jbrowse sort-gff
            cmd = f"gff3sort.pl --precise '{data}' | grep -v \"^$\" > '{dest}'"
            self.subprocess_popen(cmd, cwd=False)

            self.subprocess_check_call(["bgzip", "-f", dest], cwd=False)
            self.subprocess_check_call(["tabix", "-f", "-p", "gff", dest + ".gz"], cwd=False)

    def _sort_bed(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = ["sort", "-k1,1", "-k2,2n", data]
            with open(dest, "w") as handle:
                self.subprocess_check_call(cmd, output=handle)

            self.subprocess_check_call(["bgzip", "-f", dest])
            self.subprocess_check_call(["tabix", "-f", "-p", "bed", dest + ".gz"])

    def process_annotations(self, track, parent):
        _parent_genome = parent["label"]
        category = track["category"].replace("__pd__date__pd__", TODAY)

        track_labels = []

        for i, (
            dataset_path,
            dataset_ext,
            track_human_label,
            extra_metadata,
        ) in enumerate(track["trackfiles"]):
            # Unsanitize labels (element_identifiers are always sanitized by Galaxy)
            for key, value in mapped_chars.items():
                track_human_label = track_human_label.replace(value, key)

            log.info(
                f"-----> Processing track {category} / {track_human_label} ({dataset_ext}, {len(dataset_path) if is_multi else 1} files)"
            )

            outputTrackConfig = {
                "category": category,
            }

            outputTrackConfig["key"] = track_human_label
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
                parent["uniq_id"],
            ]
            hashData = "|".join(hashData).encode("utf-8")
            outputTrackConfig["label"] = hashlib.md5(hashData).hexdigest() + f"_{track['track_num']}_{i}"
            outputTrackConfig["metadata"] = extra_metadata

            outputTrackConfig["style"] = track["style"]

            outputTrackConfig["formatdetails"] = track["formatdetails"]

            print(f"Adding track {dataset_ext}")
            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    parent,
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                )
            elif dataset_ext == "bed":
                self.add_bed(
                    parent,
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                )
            elif dataset_ext == "bigwig":
                self.add_bigwig(
                    parent,
                    dataset_path, outputTrackConfig, track["conf"]["options"]["wiggle"]
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

                self.add_xam(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["pileup"],
                    index=real_indexes[i],
                    ext="bam",
                )
            elif dataset_ext == "cram":
                real_indexes = track["conf"]["options"]["cram"]["cram_indices"][
                    "cram_index"
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

                self.add_xam(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["cram"],
                    index=real_indexes[i],
                    ext="cram",
                )
            elif dataset_ext == "vcf":
                self.add_vcf(parent, dataset_path, outputTrackConfig)
            elif dataset_ext == "vcf_bgzip":
                self.add_vcf(parent, dataset_path, outputTrackConfig, zipped=True)
            elif dataset_ext == "rest":
                self.add_rest(
                    parent,
                    track["conf"]["options"]["rest"]["url"], outputTrackConfig
                )
            elif dataset_ext == "paf":
                log.debug("===== PAF =====")
                self.add_paf(
                    parent,
                    dataset_path, outputTrackConfig, track["conf"]["options"]["synteny"],
                    _parent_genome
                )
            elif dataset_ext in ("hic", "juicebox_hic"):
                self.add_hic(
                    parent,
                    dataset_path, outputTrackConfig, track["conf"]["options"]["hic"]
                )
            elif dataset_ext == "sparql":
                sparql_query = track["conf"]["options"]["sparql"]["query"]
                for key, value in mapped_chars.items():
                    sparql_query = sparql_query.replace(value, key)
                sparql_query_refnames = track["conf"]["options"]["sparql"][
                    "query_refnames"
                ]
                for key, value in mapped_chars.items():
                    sparql_query_refnames = sparql_query_refnames.replace(value, key)
                self.add_sparql(
                    parent,
                    track["conf"]["options"]["sparql"]["url"],
                    sparql_query,
                    sparql_query_refnames,
                    outputTrackConfig,
                )
            else:
                log.error(f"Do not know how to handle {dataset_ext}")

            track_labels.append(outputTrackConfig["label"])

        # Return non-human label for use in other fields
        return track_labels

    def add_default_view(self, genome, default_loc, tracks_on):

        refName = ""
        start = end = None
        if default_loc:
            loc_match = re.search(r"^(\w+):(\d+)\.+(\d+)$", default_loc)
            if loc_match:
                refName = loc_match.group(1)
                start = int(loc_match.group(2))
                end = int(loc_match.group(3))

        if not refName and self.assembly_ids[genome['uniq_id']] is not None:
            refName = self.assembly_ids[genome['uniq_id']]

        if start and end:
            loc_str = f"{refName}:{start}-{end}"
        else:
            loc_str = refName

        # Updating an existing jbrowse instance, merge with pre-existing view
        if self.update:
            for existing in self.default_views.values():
                if existing['init']['assembly'] == genome['uniq_id']:
                    view_specs = existing
                    if loc_str:
                        view_specs['init']['loc'] = loc_str
                    view_specs['init']['tracks'].extend(tracks_on)

        else:
            view_specs = {
                "type": "LinearGenomeView",
                "init": {
                    "assembly": genome['uniq_id'],
                    "loc": loc_str,
                    "tracks": tracks_on
                }
            }

        return view_specs

    def add_default_session(self, default_views):
        """
        Add some default session settings: set some assemblies/tracks on/off

        This allows to select a default view:
        - jb type (Linear, Circular, etc)
        - default location on an assembly
        - default tracks
        - ...

        Now using this method:
        https://github.com/GMOD/jbrowse-components/pull/4907

        Different methods that were tested/discussed earlier:
        - using a defaultSession item in config.json before PR 4970: this proved to be difficult:
          forced to write a full session block, including hard-coded/hard-to-guess items,
          no good way to let Jbrowse2 display a scaffold without knowing its size
        - using JBrowse2 as an embedded React component in a tool-generated html file:
          it works but it requires generating js code to actually do what we want = chosing default view, assembly, tracks, ...
        - writing a session-spec inside the config.json file: this is not yet supported as of 2.10.2 (see PR 4148 below)
          a session-spec is a kind of simplified defaultSession where you don't need to specify every aspect of the session
        - passing a session-spec through URL params by embedding the JBrowse2 index.html inside an iframe

        Xrefs to understand the choices:
        https://github.com/GMOD/jbrowse-components/issues/2708
        https://github.com/GMOD/jbrowse-components/discussions/3568
        https://github.com/GMOD/jbrowse-components/pull/4148
        """

        session_name = ', '.join(x['init']['assembly'] for x in default_views.values())

        session_spec = {
            "name": session_name,
            "views": list(default_views.values())
        }

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

        config_json["defaultSession"].update(session_spec)

        with open(config_path, "w") as config_file:
            json.dump(config_json, config_file, indent=2)

    def add_general_configuration(self, data):
        """
        Add some general configuration to the config.json file
        """

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

        config_data = {}

        config_data["disableAnalytics"] = data.get("analytics", "false") == "true"

        config_data["theme"] = {
            "palette": {
                "primary": {"main": data.get("primary_color", "#0D233F")},
                "secondary": {"main": data.get("secondary_color", "#721E63")},
                "tertiary": {"main": data.get("tertiary_color", "#135560")},
                "quaternary": {"main": data.get("quaternary_color", "#FFB11D")},
            },
            "typography": {"fontSize": int(data.get("font_size", 10))},
        }

        config_json["configuration"].update(config_data)

        with open(config_path, "w") as config_file:
            json.dump(config_json, config_file, indent=2)

    def clone_jbrowse(self, jbrowse_dir, destination):
        """
            Clone a JBrowse directory into a destination directory.

            Not using `jbrowse create` command to allow running on internet-less compute + to make sure code is frozen
        """

        copytree(jbrowse_dir, destination)
        try:
            shutil.rmtree(os.path.join(destination, "test_data"))
        except OSError as e:
            log.error(f"Error: {e.filename} - {e.strerror}.")

        if not os.path.exists(os.path.join(destination, "data")):
            # It can already exist if upgrading an instance
            os.makedirs(os.path.join(destination, "data"))
            log.info(f"makedir {os.path.join(destination, 'data')}")

        os.symlink("./data/config.json", os.path.join(destination, "config.json"))


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def parse_style_conf(item):
    if "type" in item.attrib and item.attrib["type"] in ["boolean", "integer"]:
        if item.attrib["type"] == "boolean":
            return item.text in ("yes", "true", "True")
        elif item.attrib["type"] == "integer":
            return int(item.text)
    else:
        return item.text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("xml", type=argparse.FileType("r"), help="Track Configuration")

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument("--update", help="Update an existing JBrowse2 instance", action="store_true")
    parser.add_argument("--outdir", help="Output directory", default="out")
    args = parser.parse_args()

    tree = ET.parse(args.xml.name)
    real_root = tree.getroot()

    # This should be done ASAP
    # Sometimes this comes as `localhost` without a protocol
    GALAXY_INFRASTRUCTURE_URL = real_root.find("metadata/galaxyUrl").text
    if not GALAXY_INFRASTRUCTURE_URL.startswith("http"):
        # so we'll prepend `http://` and hope for the best. Requests *should*
        # be GET and not POST so it should redirect OK
        GALAXY_INFRASTRUCTURE_URL = "http://" + GALAXY_INFRASTRUCTURE_URL

    jc = JbrowseConnector(
        jbrowse=args.jbrowse,
        outdir=args.outdir,
        update=args.update,
    )

    for assembly in real_root.findall("assembly"):
        genome_el = assembly.find('genome')
        genome = {
            "path": os.path.realpath(genome_el.attrib["path"]),
            "meta": metadata_from_node(genome_el.find("metadata")),
            "label": genome_el.attrib["label"],
        }

        log.debug("Processing genome", genome)
        genome["uniq_id"] = jc.add_assembly(genome["path"], genome["label"])

        # TODO add metadata to tracks
        track_num = 0
        for track in assembly.findall("tracks/track"):
            track_conf = {}
            track_conf["trackfiles"] = []
            track_conf["track_num"] = track_num

            default_tracks_on = []

            is_multi_bigwig = False
            try:
                if track.find("options/wiggle/multibigwig") and (
                    track.find("options/wiggle/multibigwig").text == "True"
                ):
                    is_multi_bigwig = True
                    multi_bigwig_paths = []
            except KeyError:
                pass

            trackfiles = track.findall("files/trackFile") or []
            if trackfiles:
                for x in track.findall("files/trackFile"):
                    if is_multi_bigwig:
                        multi_bigwig_paths.append(
                            (x.attrib["label"], os.path.realpath(x.attrib["path"]))
                        )
                    else:
                        metadata = metadata_from_node(x.find("metadata"))
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
            track_conf["style"] = {
                item.tag: parse_style_conf(item) for item in (track.find("options/style") or [])
            }

            track_conf["style"] = {
                item.tag: parse_style_conf(item) for item in (track.find("options/style") or [])
            }

            track_conf["style_labels"] = {
                item.tag: parse_style_conf(item)
                for item in (track.find("options/style_labels") or [])
            }
            track_conf["formatdetails"] = {
                item.tag: parse_style_conf(item) for item in (track.find("options/formatdetails") or [])
            }

            track_conf["conf"] = etree_to_dict(track.find("options"))

            track_labels = jc.process_annotations(track_conf, genome)

            if track.attrib["visibility"] == "default_on":
                for tlabel in track_labels:
                    default_tracks_on.append(tlabel)

            default_loc = assembly.find("defaultLocation").text
            jc.default_views[genome['uniq_id']] = jc.add_default_view(genome, default_loc, default_tracks_on)

            track_num += 1

    general_data = {
        "analytics": real_root.find("metadata/general/analytics").text,
        "primary_color": real_root.find("metadata/general/primary_color").text,
        "secondary_color": real_root.find("metadata/general/secondary_color").text,
        "tertiary_color": real_root.find("metadata/general/tertiary_color").text,
        "quaternary_color": real_root.find("metadata/general/quaternary_color").text,
        "font_size": real_root.find("metadata/general/font_size").text,
    }

    jc.add_default_session(jc.default_views)
    jc.add_general_configuration(general_data)
    jc.text_index()
