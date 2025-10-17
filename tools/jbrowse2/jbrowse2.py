#!/usr/bin/env python
import argparse
import csv
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

import requests


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
    "\n": "__cn__",
    "\r": "__cr__",
    "\t": "__tc__",
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


def unsanitize(input):

    for key, value in mapped_chars.items():
        input = input.replace(value, key)

    return input


def metadata_from_node(node):
    metadata = {}

    if len(node.findall("dataset")) == 1:

        for key, value in node.findall("dataset")[0].attrib.items():
            metadata[f"dataset_{key}"] = value

        for key, value in node.findall("history")[0].attrib.items():
            metadata[f"history_{key}"] = value

        for key, value in node.findall("metadata")[0].attrib.items():
            metadata[f"metadata_{key}"] = unsanitize(value)

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

    # Load additional metadata from a TSV file if any given by user
    bonus = node.findall("bonus")
    if bonus and "src" in bonus[0].attrib and bonus[0].attrib["src"]:
        with open(bonus[0].attrib["src"], "r") as bonus_tsv:
            bonus_content = csv.reader(bonus_tsv, delimiter="\t", quotechar='"')
            for row in bonus_content:
                if len(row) == 2:
                    if row[0] in metadata:
                        log.warning(f"Overwriting existing metadata {row[0]} with value from bonus file {row[1]}")
                    metadata[row[0]] = row[1]
                else:
                    log.warning(f"Skipping invalid bonus metadata line: {row}")

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

        self.plugins = []

        self.use_synteny_viewer = False

        self.synteny_tracks = []

        self.clone_jbrowse(self.jbrowse, self.outdir)

        # If upgrading, look at the existing data
        self.check_existing(self.outdir)

    def get_cwd(self, cwd):
        if cwd:
            return self.outdir
        else:
            return subprocess.check_output(['pwd']).decode('utf-8').strip()
            # return None

    def subprocess_check_call(self, command, output=None, cwd=True):
        if output:
            log.debug(f"cd {self.get_cwd(cwd)} && {' '.join(command)} > {output.name}")
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
        # Use to support symlinking in jbrowse1, in jbrowse2 prefer to use remote uri
        cmd = ["cp", src, dest]

        return self.subprocess_check_call(cmd)

    def _prepare_track_style(self, xml_conf):
        style_data = {
            "type": "LinearBasicDisplay",  # No ideal default, but should be overwritten anyway
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

        elif display_type == "MultiLinearWiggleDisplay":

            wig_renderer = xml_conf.get("renderer", "multirowxy")
            style_data["defaultRendering"] = wig_renderer

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
            feat_jexl = unsanitize(xml_conf["formatdetails"]["feature"])
            formatDetails["feature"] = feat_jexl

        if "subfeature" in xml_conf["formatdetails"]:
            sfeat_jexl = unsanitize(xml_conf["formatdetails"]["subfeature"])
            formatDetails["subfeatures"] = sfeat_jexl

        if "depth" in xml_conf["formatdetails"]:
            formatDetails["depth"] = int(xml_conf["formatdetails"]["depth"])

        return {"formatDetails": formatDetails}

    def _prepare_track_metadata(self, xml_conf):
        metadata = {
        }

        metadata = xml_conf["metadata"]

        return {"metadata": metadata}

    def check_existing(self, destination):
        existing = os.path.join(destination, "config.json")
        if os.path.exists(existing):
            with open(existing, "r") as existing_conf:
                conf = json.load(existing_conf)
                if "assemblies" in conf:
                    for assembly in conf["assemblies"]:
                        if "name" in assembly:

                            # Look for a default scaffold
                            default_seq = None
                            if 'defaultSession' in conf and 'views' in conf['defaultSession']:
                                for view in conf['defaultSession']['views']:
                                    if 'init' in view and 'assembly' in view['init'] and 'loc' in view['init']:
                                        if view['init']['assembly'] == assembly["name"]:
                                            default_seq = view['init']['loc'].split(":")[0]
                                    if "views" in view:
                                        subviews = view["views"]
                                        for subview in subviews:
                                            if 'init' in subview and 'assembly' in subview['init'] and 'loc' in subview['init']:
                                                if subview['init']['assembly'] == assembly["name"]:
                                                    default_seq = subview['init']['loc'].split(":")[0]

                            self.assembly_ids[assembly["name"]] = default_seq

    def _load_old_genome_views(self):

        views = {}

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

            # Find default synteny views existing from a previous jbrowse dataset
            if 'defaultSession' in config_json and 'views' in config_json['defaultSession']:
                for view in config_json['defaultSession']['views']:
                    if view['type'] != "LinearSyntenyView":
                        if 'init' in view and 'assembly' in view['init']:
                            views[view['init']['assembly']] = view

        return views

    def _load_old_synteny_views(self):

        views = []

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

            # Find default synteny views existing from a previous jbrowse dataset
            if 'defaultSession' in config_json and 'views' in config_json['defaultSession']:
                for view in config_json['defaultSession']['views']:
                    if view['type'] == "LinearSyntenyView":
                        views.append(view)

        return views

    def add_assembly(self, path, label, is_remote=False, cytobands=None, ref_name_aliases=None):

        if not is_remote:
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

                        # Find default views existing for this assembly
                        if 'defaultSession' in config_json and 'views' in config_json['defaultSession']:
                            for view in config_json['defaultSession']['views']:
                                if 'init' in view and 'assembly' in view['init']:
                                    if view['init']['assembly'] == label:

                                        log.info("Found existing assembly from existing JBrowse2 instance, preserving it")

                                        self.default_views[view['init']['assembly']] = view

                        return label

        # Copy ref alias file if any
        if ref_name_aliases:
            copied_ref_name_aliases = seq_path + ".aliases"
            shutil.copy(ref_name_aliases, copied_ref_name_aliases)
            copied_ref_name_aliases = rel_seq_path + ".aliases"

        # Copy cytobands file if any
        if cytobands:
            copied_cytobands = seq_path + ".cytobands"
            shutil.copy(cytobands, copied_cytobands)
            copied_cytobands = rel_seq_path + ".cytobands"

        # Find a non-existing label for the new genome
        # (to avoid colision when upgrading an existing instance)
        lab_try = 1
        uniq_label = label
        while uniq_label in self.assembly_ids:
            uniq_label = label + str(lab_try)
            lab_try += 1

        if is_remote:

            # Find a default scaffold to display
            with requests.get(path + ".fai", stream=True) as response:
                response.raise_for_status()
                first_seq = next(response.iter_lines())
                first_seq = first_seq.decode("utf-8").split('\t')[0]

            self.assembly_ids[uniq_label] = first_seq

            # We assume we just need to suffix url with .fai and .gzi for indexes.
            cmd_jb = [
                "jbrowse",
                "add-assembly",
                "--name",
                uniq_label,
                "--type",
                "bgzipFasta",
                "--out",
                self.outdir,
                "--skipCheck",
            ]

            if ref_name_aliases:
                cmd_jb.extend([
                    "--refNameAliases",
                    copied_ref_name_aliases,
                ])

            cmd_jb.append(path)  # Path is an url in remote mode

            self.subprocess_check_call(cmd_jb)
        else:
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
            ]

            if ref_name_aliases:
                cmd_jb.extend([
                    "--refNameAliases",
                    copied_ref_name_aliases,
                ])

            cmd_jb.append(rel_seq_path + ".fasta.gz")

            self.subprocess_check_call(cmd_jb)

        if cytobands:
            self.add_cytobands(uniq_label, copied_cytobands)

        return uniq_label

    def add_cytobands(self, assembly_name, cytobands_path):

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

        config_data = {}

        config_data["cytobands"] = {
            "adapter": {
                "type": "CytobandAdapter",
                "cytobandLocation": {
                    "uri": cytobands_path
                }
            }
        }

        filled_assemblies = []
        for assembly in config_json["assemblies"]:
            if assembly["name"] == assembly_name:
                assembly.update(config_data)
            filled_assemblies.append(assembly)
        config_json["assemblies"] = filled_assemblies

        with open(config_path, "w") as config_file:
            json.dump(config_json, config_file, indent=2)

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

    def add_gc_content(self, parent, trackData, **kwargs):

        adapter = {}
        existing = os.path.join(self.outdir, "config.json")
        if os.path.exists(existing):
            with open(existing, "r") as existing_conf:
                conf = json.load(existing_conf)
                if "assemblies" in conf:
                    for assembly in conf["assemblies"]:
                        if assembly.get('name', "") == parent['uniq_id']:
                            adapter = assembly.get('sequence', {}).get('adapter', {})

        json_track_data = {
            "type": "GCContentTrack",
            "trackId": trackData["label"],
            "name": trackData["key"],
            "adapter": adapter,
            "category": [trackData["category"]],
            "assemblyNames": [parent['uniq_id']],
        }

        style_json = self._prepare_track_style(trackData)

        json_track_data.update(style_json)

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json.dumps(json_track_data),
            ]
        )

    def add_bigwig(self, parent, data, trackData, wiggleOpts, **kwargs):

        if trackData['remote']:
            rel_dest = data
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".bw")
            dest = os.path.join(self.outdir, rel_dest)
            self.symlink_or_copy(os.path.realpath(data), dest)

        style_json = self._prepare_track_style(trackData)

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_bigwig_multi(self, parent, data_files, trackData, wiggleOpts, **kwargs):

        subadapters = []

        sub_num = 0
        for data in data_files:
            if trackData['remote']:
                rel_dest = data[1]
            else:
                rel_dest = os.path.join("data", f"{trackData['label']}_sub{sub_num}.bw")
                dest = os.path.join(self.outdir, rel_dest)
                self.symlink_or_copy(os.path.realpath(data[1]), dest)

            subadapters.append({
                "type": "BigWigAdapter",
                "name": data[0],
                "bigWigLocation": {
                    "uri": rel_dest,
                    "locationType": "UriLocation"
                }
            })
            sub_num += 1

        json_track_data = {
            "type": "MultiQuantitativeTrack",
            "trackId": trackData["label"],
            "name": trackData["key"],
            "adapter": {
                "type": "MultiWiggleAdapter",
                "subadapters": subadapters
            },
            "category": [trackData["category"]],
            "assemblyNames": [parent['uniq_id']],
        }

        style_json = self._prepare_track_style(trackData)

        json_track_data.update(style_json)

        track_metadata = self._prepare_track_metadata(trackData)

        json_track_data.update(track_metadata)

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json.dumps(json_track_data),
            ]
        )

    # Anything ending in "am" (Bam or Cram)
    def add_xam(self, parent, data, trackData, xamOpts, index=None, ext="bam", **kwargs):
        index_ext = "bai"
        if ext == "cram":
            index_ext = "crai"

        if trackData['remote']:
            rel_dest = data
            # Index will be set automatically as xam url + xai .suffix by add-track cmd
        else:
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

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_vcf(self, parent, data, trackData, vcfOpts={}, zipped=False, **kwargs):
        if trackData['remote']:
            rel_dest = data
        else:
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

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_gff(self, parent, data, format, trackData, gffOpts, **kwargs):
        if trackData['remote']:
            rel_dest = data
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".gff")
            dest = os.path.join(self.outdir, rel_dest)
            rel_dest = rel_dest + ".gz"

            self._sort_gff(data, dest)

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        if gffOpts.get('index', 'false') in ("yes", "true", "True"):
            if parent['uniq_id'] not in self.tracksToIndex:
                self.tracksToIndex[parent['uniq_id']] = []
            self.tracksToIndex[parent['uniq_id']].append(trackData["label"])

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_gtf(self, parent, data, format, trackData, gffOpts, **kwargs):
        # Not a super recommended format
        # https://github.com/GMOD/jbrowse-components/pull/2389
        # https://github.com/GMOD/jbrowse-components/issues/3876
        if trackData['remote']:
            rel_dest = data
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".gtf")
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), dest)

        json_track_data = {
            "type": "FeatureTrack",
            "trackId": trackData["label"],
            "name": trackData["key"],
            "adapter": {
                "type": "GtfAdapter",
                "gtfLocation": {
                    "uri": rel_dest,
                    "locationType": "UriLocation"
                },
            },
            "category": [trackData["category"]],
            "assemblyNames": [parent['uniq_id']],
        }

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        json_track_data.update(style_json)

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json.dumps(json_track_data),
            ]
        )

    def add_bed(self, parent, data, format, trackData, gffOpts, **kwargs):
        if trackData['remote']:
            rel_dest = data
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".bed")
            dest = os.path.join(self.outdir, rel_dest)
            rel_dest = rel_dest + ".gz"

            self._sort_bed(data, dest)

        style_json = self._prepare_track_style(trackData)

        formatdetails = self._prepare_format_details(trackData)

        style_json.update(formatdetails)

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        if gffOpts.get('index', 'false') in ("yes", "true", "True"):
            if parent['uniq_id'] not in self.tracksToIndex:
                self.tracksToIndex[parent['uniq_id']] = []
            self.tracksToIndex[parent['uniq_id']].append(trackData["label"])

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_paf(self, parent, data, trackData, pafOpts, **kwargs):

        if trackData['remote']:
            rel_dest = data

            if rel_dest.endswith('pif') or rel_dest.endswith('pif.gz'):
                adapter = "pif"
            else:
                adapter = "paf"
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".pif.gz")
            dest = os.path.join(self.outdir, rel_dest)

            cmd = ["jbrowse", "make-pif", "--out", dest, os.path.realpath(data)]
            self.subprocess_check_call(cmd)

            adapter = "pif"

        if trackData["style"]["display"] == "LinearBasicDisplay":
            # Normal style track

            json_track_data = {
                "type": "SyntenyTrack",
                "trackId": trackData["label"],
                "name": trackData["key"],
                "adapter": {
                    "type": "PairwiseIndexedPAFAdapter",
                    "pifGzLocation": {
                        "uri": rel_dest,
                    },
                    "index": {
                        "location": {
                            "uri": rel_dest + ".tbi",
                        }
                    },
                },
                "category": [trackData["category"]],
                "assemblyNames": [parent['uniq_id']],
            }
        else:
            # Synteny viewer

            json_track_data = {
                "type": "SyntenyTrack",
                "trackId": trackData["label"],
                "name": trackData["key"],
                "adapter": {
                    "assemblyNames": [
                        parent['uniq_id'],
                        "",  # Placeholder until we know the next genome id
                    ],
                },
                "category": [trackData["category"]],
                "assemblyNames": [
                    parent['uniq_id'],
                    "",  # Placeholder until we know the next genome id
                ]
            }

            if adapter == "pif":
                json_track_data["adapter"].update({
                    "type": "PairwiseIndexedPAFAdapter",
                    "pifGzLocation": {
                        "uri": rel_dest,
                    },
                    "index": {
                        "location": {
                            "uri": rel_dest + ".tbi",
                        }
                    },
                })
            else:
                json_track_data["adapter"].update({
                    "type": "PAFAdapter",
                    "pafLocation": {
                        "uri": rel_dest,
                    },
                })

        style_json = self._prepare_track_style(trackData)

        json_track_data.update(style_json)

        track_metadata = self._prepare_track_metadata(trackData)

        json_track_data.update(track_metadata)

        if trackData["style"]["display"] == "LinearBasicDisplay":
            self.subprocess_check_call(
                [
                    "jbrowse",
                    "add-track-json",
                    "--target",
                    self.outdir,
                    json.dumps(json_track_data),
                ]
            )
        else:
            self.synteny_tracks.append(json_track_data)

    def add_hic(self, parent, data, trackData, hicOpts, **kwargs):
        if trackData['remote']:
            rel_dest = data
        else:
            rel_dest = os.path.join("data", trackData["label"] + ".hic")
            dest = os.path.join(self.outdir, rel_dest)
            self.symlink_or_copy(os.path.realpath(data), dest)

        style_json = self._prepare_track_style(trackData)

        track_metadata = self._prepare_track_metadata(trackData)

        style_json.update(track_metadata)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            parent,
            config=style_json,
            remote=trackData['remote']
        )

    def add_maf(self, parent, data, trackData, mafOpts, **kwargs):

        # Add needed plugin
        plugin_def = {
            "name": "MafViewer",
            "url": "https://unpkg.com/jbrowse-plugin-mafviewer/dist/jbrowse-plugin-mafviewer.umd.production.min.js"
        }
        self.plugins.append(plugin_def)

        rel_dest = os.path.join("data", trackData["label"] + ".maf")
        dest = os.path.join(self.outdir, rel_dest)

        assembly_name = mafOpts.get("assembly_name", "")
        if not assembly_name:
            # Guess from assembly
            assembly_name = parent['uniq_id']

        self._convert_maf(data, dest, assembly_name)

        # Extract samples list
        mafs = open(data, "r").readlines()
        mafss = [x for x in mafs if (x.startswith("s\t") or x.startswith("s "))]
        samp = [x.split()[1] for x in mafss if len(x.split()) > 0]
        sampu = list(dict.fromkeys(samp))
        samples = [x.split(".")[0] for x in sampu]
        samples.sort()

        json_track_data = {
            "type": "MafTrack",
            "trackId": trackData["label"],
            "name": trackData["key"],
            "adapter": {
                "type": "MafTabixAdapter",
                "samples": samples,
                "bedGzLocation": {
                    "uri": rel_dest + ".gz",
                },
                "index": {
                    "location": {
                        "uri": rel_dest + ".gz.tbi",
                    },
                },
            },
            "category": [trackData["category"]],
            "assemblyNames": [parent['uniq_id']],
        }

        style_json = self._prepare_track_style(trackData)

        json_track_data.update(style_json)

        track_metadata = self._prepare_track_metadata(trackData)

        json_track_data.update(track_metadata)

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json.dumps(json_track_data),
            ]
        )

    def add_sparql(self, parent, url, query, query_refnames, trackData):
        json_track_data = {
            "type": "FeatureTrack",
            "trackId": trackData["label"],
            "name": trackData["key"],
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

        # TODO handle metadata somehow for sparql too

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                self.outdir,
                json.dumps(json_track_data),
            ]
        )

    def _add_track(self, track_id, label, category, path, assembly, config=None, trackType=None, load_action="inPlace", assemblies=None, remote=False):
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
            "--name",
            label,
            "--category",
            category,
            "--target",
            self.outdir,
            "--trackId",
            track_id,
            "--assemblyNames",
            assemblies if assemblies else assembly['uniq_id'],
        ]

        if not remote:
            cmd.append("--load")
            cmd.append(load_action)

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
            # Not using jbrowse sort-gff because it uses sort and has the problem exposed on https://github.com/tao-bioinfo/gff3sort
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

    def _convert_maf(self, data, dest, assembly_name):
        # Only convert if not already done
        if not os.path.exists(dest):

            dest_bed = dest + ".bed"
            cmd = ["python", os.path.join(SELF_LOCATION, "maf2bed.py"), assembly_name, data, dest_bed]
            self.subprocess_check_call(cmd, cwd=False)

            cmd = ["sort", "-k1,1", "-k2,2n", dest_bed]
            with open(dest, "w") as handle:
                self.subprocess_check_call(cmd, output=handle)

            self.subprocess_check_call(["bgzip", "-f", dest], cwd=False)
            self.subprocess_check_call(["tabix", "-f", "-p", "bed", dest + ".gz"], cwd=False)

    def process_annotations(self, track, parent):
        category = track["category"].replace("__pd__date__pd__", TODAY)

        track_labels = []

        for i, (
            dataset_path,
            dataset_ext,
            track_human_label,
            extra_metadata,
        ) in enumerate(track["trackfiles"]):
            # Unsanitize labels (element_identifiers are always sanitized by Galaxy)
            track_human_label = unsanitize(track_human_label)

            is_multi = type(dataset_path) is list

            log.info(
                f"-----> Processing track {category} / {track_human_label} ({dataset_ext}, {len(dataset_path) if is_multi else 1} files)"
            )

            outputTrackConfig = {
                "category": category,
            }

            outputTrackConfig["key"] = track_human_label
            # We add extra data to hash for the case of non-file tracks
            if (
                "conf" in track
                and "options" in track["conf"]
                and "url" in track["conf"]["options"]
            ):
                non_file_info = track["conf"]["options"]["url"]
            else:
                non_file_info = ""

            # I chose to use track['category'] instead of 'category' here. This
            # is intentional. This way re-running the tool on a different date
            # will not generate different hashes and make comparison of outputs
            # much simpler.
            hashData = [
                str(dataset_path),
                track_human_label,
                track["category"],
                non_file_info,
                parent["uniq_id"],
            ]
            hashData = "|".join(hashData).encode("utf-8")
            outputTrackConfig["label"] = hashlib.md5(hashData).hexdigest() + f"_{track['track_num']}_{i}"
            outputTrackConfig["metadata"] = extra_metadata

            outputTrackConfig["style"] = track["style"]

            outputTrackConfig["formatdetails"] = track["formatdetails"]

            outputTrackConfig["remote"] = track["remote"]

            # Guess extension for remote data
            if dataset_ext == "gff,gff3,bed":
                if dataset_path.endswith(".bed") or dataset_path.endswith(".bed.gz"):
                    dataset_ext = "bed"
                else:
                    dataset_ext = "gff"
            elif dataset_ext == "vcf,vcf_bgzip":
                if dataset_path.endswith(".vcf.gz"):
                    dataset_ext = "vcf_bgzip"
                else:
                    dataset_ext = "vcf"

            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    parent,
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                )
            elif dataset_ext in ("gtf"):
                self.add_gtf(
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
                if is_multi:
                    self.add_bigwig_multi(
                        parent,
                        dataset_path, outputTrackConfig, track["conf"]["options"]["wiggle"]
                    )
                else:
                    self.add_bigwig(
                        parent,
                        dataset_path, outputTrackConfig, track["conf"]["options"]["wiggle"]
                    )
            elif dataset_ext == "maf":
                self.add_maf(
                    parent,
                    dataset_path, outputTrackConfig, track["conf"]["options"]["maf"]
                )
            elif dataset_ext == "bam":

                if track["remote"]:
                    bam_index = dataset_path + '.bai'
                else:
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

                    bam_index = real_indexes[i]

                self.add_xam(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["pileup"],
                    index=bam_index,
                    ext="bam",
                )
            elif dataset_ext == "cram":

                if track["remote"]:
                    cram_index = dataset_path + '.crai'
                else:
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

                    cram_index = real_indexes[i]

                self.add_xam(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["cram"],
                    index=cram_index,
                    ext="cram",
                )
            elif dataset_ext == "vcf":
                self.add_vcf(
                    parent,
                    dataset_path,
                    outputTrackConfig
                )
            elif dataset_ext == "vcf_bgzip":
                self.add_vcf(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    zipped=True
                )
            elif dataset_ext == "paf":  # https://fr.wikipedia.org/wiki/Paf_le_chien
                self.add_paf(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["synteny"]
                )
            elif dataset_ext in ("hic"):
                self.add_hic(
                    parent,
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["hic"]
                )
            elif dataset_ext == "sparql":
                sparql_query = unsanitize(track["conf"]["options"]["sparql"]["query"])
                sparql_query_refnames = track["conf"]["options"]["sparql"].get("query_refnames", "")
                if sparql_query_refnames:
                    sparql_query_refnames = unsanitize(sparql_query_refnames)
                self.add_sparql(
                    parent,
                    track["conf"]["options"]["sparql"]["url"],
                    sparql_query,
                    sparql_query_refnames,
                    outputTrackConfig,
                )
            elif dataset_ext == "gc":
                self.add_gc_content(
                    parent,
                    outputTrackConfig,
                )
            else:
                raise RuntimeError(f"Do not know how to handle dataset of type '{dataset_ext}'")

            track_labels.append(outputTrackConfig["label"])

        # Return non-human label for use in other fields
        return track_labels

    def add_default_view_genome(self, genome, default_loc, tracks_on):

        refName = ""
        start = end = None
        if default_loc:
            loc_match = re.search(r"^(\w+):(\d+)\.+(\d+)$", default_loc)
            if loc_match:
                refName = loc_match.group(1)
                start = int(loc_match.group(2))
                end = int(loc_match.group(3))

        if not refName and self.assembly_ids[genome['uniq_id']]:
            refName = self.assembly_ids[genome['uniq_id']]

        if start and end:
            loc_str = f"{refName}:{start}-{end}"
        else:
            loc_str = refName

        # Updating an existing jbrowse instance, merge with pre-existing view
        view_specs = None
        if self.update:
            for existing in self.default_views.values():
                if len(existing) and existing["type"] == "LinearGenomeView":
                    if existing['init']['assembly'] == genome['uniq_id']:
                        view_specs = existing
                        if loc_str:
                            view_specs['init']['loc'] = loc_str
                        view_specs['init']['tracks'].extend(tracks_on)

        if view_specs is None:  # Not updating, or updating from synteny
            view_specs = {
                "type": "LinearGenomeView",
                "init": {
                    "assembly": genome['uniq_id'],
                    "loc": loc_str,
                    "tracks": tracks_on
                }
            }

        return view_specs

    def add_default_view_synteny(self, genome_views, synteny_tracks):

        # Add json for cached synteny tracks
        # We cache them because we need to know the target genome uniq_id
        for strack in synteny_tracks:

            # Target assembly is the next genome, find its uniq_id
            query_assembly = strack["assemblyNames"][0]
            ass_uniq_ids = list(self.assembly_ids.keys())
            query_index = ass_uniq_ids.index(query_assembly)
            target_assembly = ass_uniq_ids[query_index + 1]

            strack["assemblyNames"][1] = target_assembly
            strack["adapter"]["assemblyNames"][1] = target_assembly

            self.subprocess_check_call(
                [
                    "jbrowse",
                    "add-track-json",
                    "--target",
                    self.outdir,
                    json.dumps(strack),
                ]
            )

        # Configure the synteny view
        levels = []

        for strack in synteny_tracks:
            lev = {
                "type": "LinearSyntenyViewHelper",
                "tracks": [
                    {
                        "type": "SyntenyTrack",
                        "configuration": strack["trackId"],
                        "displays": [
                            {
                                "type": "LinearSyntenyDisplay",
                                "configuration": strack["trackId"] + "_LinearSyntenyDisplay"
                            }
                        ]
                    }
                ],
                "height": 100,
                "level": len(levels)
            }
            levels.append(lev)

        view_specs = {
            "type": "LinearSyntenyView",
            "views": genome_views,
            "levels": levels
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

        if self.use_synteny_viewer:
            session_name = "Synteny"
        else:
            session_name = ', '.join(x['init']['assembly'] for x in default_views)

        session_spec = {
            "name": session_name,
            "views": default_views
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

    def add_plugins(self, data):
        """
        Add plugins to the config.json file
        """

        config_path = os.path.join(self.outdir, "config.json")
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

        if "plugins" not in config_json:
            config_json["plugins"] = []

        config_json["plugins"].extend(data)

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


def validate_synteny(real_root):

    if len(real_root.findall('assembly/tracks/track[@format="synteny"]')) == 0:
        # No synteny data, all good
        return False

    assemblies = real_root.findall("assembly")

    if len(assemblies[-1].findall('tracks/track[@format="synteny"]')) > 0 and \
       assemblies[-1].find('tracks/track[@format="synteny"]/options/style/display').text == "LinearSyntenyDisplay":
        raise RuntimeError("You should not set a synteny track on the last genome.")

    for assembly in assemblies[1:0]:
        if len(assembly.findall('tracks/track[@format="synteny"]')) != 1 and \
           assembly.find('tracks/track[@format="synteny"]/options/style/display').text == "LinearSyntenyDisplay":
            raise RuntimeError("To use the synteny viewer, you should add a synteny track to each assembly, except the last one.")

    return True


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

    # Synteny options are special, check them first
    jc.use_synteny_viewer = validate_synteny(real_root)

    for assembly in real_root.findall("assembly"):
        genome_el = assembly.find('genome')

        is_remote = genome_el.attrib.get("remote", "false") == "true"

        genome = {
            "path": genome_el.attrib["path"] if is_remote else os.path.realpath(genome_el.attrib["path"]),
            "meta": metadata_from_node(genome_el.find("metadata")),
            "label": genome_el.attrib["label"],
        }

        cytobands = None
        cytobands_el = genome_el.find("cytobands")
        if cytobands_el is not None and "path" in cytobands_el.attrib:
            cytobands = cytobands_el.attrib["path"]

        ref_name_aliases = None
        ref_name_aliases_el = genome_el.find("ref_name_aliases")
        if ref_name_aliases_el is not None and "path" in ref_name_aliases_el.attrib:
            ref_name_aliases = ref_name_aliases_el.attrib["path"]

        log.debug("Processing genome", genome)
        genome["uniq_id"] = jc.add_assembly(genome["path"], genome["label"], is_remote, cytobands, ref_name_aliases)

        default_tracks_on = []

        track_num = 0
        for track in assembly.findall("tracks/track"):
            track_conf = {}
            track_conf["trackfiles"] = []
            track_conf["track_num"] = track_num

            trackfiles = track.findall("files/trackFile") or []

            is_multi = False
            multi_paths = []
            multi_type = None
            multi_metadata = {}
            try:
                multi_in_xml = track.find("options/multitrack")
                if multi_in_xml is not None and parse_style_conf(multi_in_xml):
                    is_multi = True
                    multi_paths = []
                    multi_type = trackfiles[0].attrib["ext"]
            except KeyError:
                pass

            is_remote = False
            if trackfiles:
                for x in trackfiles:
                    if is_multi:
                        is_remote = x.attrib.get("remote", "false") == "true"
                        multi_paths.append(
                            (x.attrib["label"], x.attrib["path"] if is_remote else os.path.realpath(x.attrib["path"]))
                        )
                        multi_metadata.update(metadata_from_node(x.find("metadata")))
                    else:
                        metadata = metadata_from_node(x.find("metadata"))
                        is_remote = x.attrib.get("remote", "false") == "true"
                        track_conf["trackfiles"].append(
                            (
                                x.attrib["path"] if is_remote else os.path.realpath(x.attrib["path"]),
                                x.attrib["ext"],
                                x.attrib["label"],
                                metadata,
                            )
                        )
            else:
                # For tracks without files (sparql, gc)
                track_conf["trackfiles"].append(
                    (
                        "",  # N/A, no path for sparql or gc
                        track.attrib["format"],
                        track.find("options/label").text,
                        {},
                    )
                )

            if is_multi:
                etal_tracks_nb = len(multi_paths[1:])
                multi_label = f"{multi_paths[0][0]} + {etal_tracks_nb} other track{'s' if etal_tracks_nb > 1 else ''}"

                track_conf["trackfiles"].append(
                    (
                        multi_paths,  # Passing an array of paths to represent as one track
                        multi_type,  # First file type
                        multi_label,  # First file label
                        multi_metadata,  # Mix of all metadata for multiple bigwig => only last file metadata coming from galaxy + custom oness
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

            track_conf["remote"] = is_remote

            track_labels = jc.process_annotations(track_conf, genome)

            if track.attrib["visibility"] == "default_on" and \
               (track_conf["format"] != "synteny" or track_conf["style"]["display"] != "LinearSyntenyDisplay"):
                for tlabel in track_labels:
                    default_tracks_on.append(tlabel)

            track_num += 1

        default_loc = assembly.find("defaultLocation").text

        jc.default_views[genome['uniq_id']] = jc.add_default_view_genome(genome, default_loc, default_tracks_on)

    if jc.use_synteny_viewer:
        synteny_view = jc.add_default_view_synteny(list(jc.default_views.values()), jc.synteny_tracks)

        views_for_session = jc._load_old_synteny_views()

        views_for_session.append(synteny_view)
    else:
        old_views = jc._load_old_genome_views()

        for old_view in old_views:
            if old_view not in jc.default_views:
                jc.default_views[old_view] = old_views[old_view]

        views_for_session = list(jc.default_views.values())

    general_data = {
        "analytics": real_root.find("metadata/general/analytics").text,
        "primary_color": real_root.find("metadata/general/primary_color").text,
        "secondary_color": real_root.find("metadata/general/secondary_color").text,
        "tertiary_color": real_root.find("metadata/general/tertiary_color").text,
        "quaternary_color": real_root.find("metadata/general/quaternary_color").text,
        "font_size": real_root.find("metadata/general/font_size").text,
    }

    jc.add_default_session(views_for_session)
    jc.add_general_configuration(general_data)
    jc.add_plugins(jc.plugins)
    jc.text_index()
