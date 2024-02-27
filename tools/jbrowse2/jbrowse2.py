#!/usr/bin/env python
import argparse
import binascii
import datetime
import hashlib
import json
import logging
import os
import re
import shutil
import struct
import subprocess
import tempfile
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


class ColorScaling(object):

    COLOR_FUNCTION_TEMPLATE = """
    function(feature, variableName, glyphObject, track) {{
        var score = {score};
        {opacity}
        return 'rgba({red}, {green}, {blue}, ' + opacity + ')';
    }}
    """

    COLOR_FUNCTION_TEMPLATE_QUAL = r"""
    function(feature, variableName, glyphObject, track) {{
        var search_up = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.parent() === undefined) {{
                return;
            }}else{{
                return self(sf.parent(), attr);
            }}
        }};

        var search_down = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.children() === undefined) {{
                return;
            }}else{{
                var kids = sf.children();
                for(var child_idx in kids){{
                    var x = self(kids[child_idx], attr);
                    if(x !== undefined){{
                        return x;
                    }}
                }}
                return;
            }}
        }};

        var color = ({user_spec_color} || search_up(feature, 'color') || search_down(feature, 'color') || {auto_gen_color});
        var score = (search_up(feature, 'score') || search_down(feature, 'score'));
        {opacity}
        if(score === undefined){{ opacity = 1; }}
        var result = /^#?([a-f\d]{{2}})([a-f\d]{{2}})([a-f\d]{{2}})$/i.exec(color);
        var red = parseInt(result[1], 16);
        var green = parseInt(result[2], 16);
        var blue = parseInt(result[3], 16);
        if(isNaN(opacity) || opacity < 0){{ opacity = 0; }}
        return 'rgba(' + red + ',' + green + ',' + blue + ',' + opacity + ')';
    }}
    """

    OPACITY_MATH = {
        "linear": """
            var opacity = (score - ({min})) / (({max}) - ({min}));
        """,
        "logarithmic": """
            var opacity = Math.log10(score - ({min})) / Math.log10(({max}) - ({min}));
        """,
        "blast": """
            var opacity = 0;
            if(score == 0.0) {{
                opacity = 1;
            }} else {{
                opacity = (20 - Math.log10(score)) / 180;
            }}
        """,
    }

    BREWER_COLOUR_IDX = 0
    BREWER_COLOUR_SCHEMES = [
        (166, 206, 227),
        (31, 120, 180),
        (178, 223, 138),
        (51, 160, 44),
        (251, 154, 153),
        (227, 26, 28),
        (253, 191, 111),
        (255, 127, 0),
        (202, 178, 214),
        (106, 61, 154),
        (255, 255, 153),
        (177, 89, 40),
        (228, 26, 28),
        (55, 126, 184),
        (77, 175, 74),
        (152, 78, 163),
        (255, 127, 0),
    ]

    BREWER_DIVERGING_PALLETES = {
        "BrBg": ("#543005", "#003c30"),
        "PiYg": ("#8e0152", "#276419"),
        "PRGn": ("#40004b", "#00441b"),
        "PuOr": ("#7f3b08", "#2d004b"),
        "RdBu": ("#67001f", "#053061"),
        "RdGy": ("#67001f", "#1a1a1a"),
        "RdYlBu": ("#a50026", "#313695"),
        "RdYlGn": ("#a50026", "#006837"),
        "Spectral": ("#9e0142", "#5e4fa2"),
    }

    def __init__(self):
        self.brewer_colour_idx = 0

    def rgb_from_hex(self, hexstr):
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back
        return struct.unpack("BBB", binascii.unhexlify(hexstr))

    def min_max_gff(self, gff_file):
        min_val = None
        max_val = None
        with open(gff_file, "r") as handle:
            for line in handle:
                try:
                    value = float(line.split("\t")[5])
                    min_val = min(value, (min_val or value))
                    max_val = max(value, (max_val or value))

                    if value < min_val:
                        min_val = value

                    if value > max_val:
                        max_val = value
                except Exception:
                    pass
        return min_val, max_val

    def hex_from_rgb(self, r, g, b):
        return "#%02x%02x%02x" % (r, g, b)

    def _get_colours(self):
        r, g, b = self.BREWER_COLOUR_SCHEMES[
            self.brewer_colour_idx % len(self.BREWER_COLOUR_SCHEMES)
        ]
        self.brewer_colour_idx += 1
        return r, g, b

    def parse_menus(self, track):
        trackConfig = {"menuTemplate": [{}, {}, {}, {}]}

        if "menu" in track["menus"]:
            menu_list = [track["menus"]["menu"]]
            if isinstance(track["menus"]["menu"], list):
                menu_list = track["menus"]["menu"]

            for m in menu_list:
                tpl = {
                    "action": m["action"],
                    "label": m.get("label", "{name}"),
                    "iconClass": m.get("iconClass", "dijitIconBookmark"),
                }
                if "url" in m:
                    tpl["url"] = m["url"]
                if "content" in m:
                    tpl["content"] = m["content"]
                if "title" in m:
                    tpl["title"] = m["title"]

                trackConfig["menuTemplate"].append(tpl)

        return trackConfig

    def parse_colours(self, track, trackFormat, gff3=None):
        # Wiggle tracks have a bicolor pallete
        trackConfig = {"style": {}}
        if trackFormat == "wiggle":

            trackConfig["style"]["pos_color"] = track["wiggle"]["color_pos"]
            trackConfig["style"]["neg_color"] = track["wiggle"]["color_neg"]

            if trackConfig["style"]["pos_color"] == "__auto__":
                trackConfig["style"]["neg_color"] = self.hex_from_rgb(
                    *self._get_colours()
                )
                trackConfig["style"]["pos_color"] = self.hex_from_rgb(
                    *self._get_colours()
                )

            # Wiggle tracks can change colour at a specified place
            bc_pivot = track["wiggle"]["bicolor_pivot"]
            if bc_pivot not in ("mean", "zero"):
                # The values are either one of those two strings
                # or a number
                bc_pivot = float(bc_pivot)
            trackConfig["bicolor_pivot"] = bc_pivot
        elif "scaling" in track:
            if track["scaling"]["method"] == "ignore":
                if track["scaling"]["scheme"]["color"] != "__auto__":
                    trackConfig["style"]["color"] = track["scaling"]["scheme"]["color"]
                else:
                    trackConfig["style"]["color"] = self.hex_from_rgb(
                        *self._get_colours()
                    )
            else:
                # Scored method
                algo = track["scaling"]["algo"]
                # linear, logarithmic, blast
                scales = track["scaling"]["scales"]
                # type __auto__, manual (min, max)
                scheme = track["scaling"]["scheme"]
                # scheme -> (type (opacity), color)
                # ==================================
                # GENE CALLS OR BLAST
                # ==================================
                if trackFormat == "blast":
                    red, green, blue = self._get_colours()
                    color_function = self.COLOR_FUNCTION_TEMPLATE.format(
                        **{
                            "score": "feature._parent.get('score')",
                            "opacity": self.OPACITY_MATH["blast"],
                            "red": red,
                            "green": green,
                            "blue": blue,
                        }
                    )
                    trackConfig["style"]["color"] = color_function.replace("\n", "")
                elif trackFormat == "gene_calls":
                    # Default values, based on GFF3 spec
                    min_val = 0
                    max_val = 1000
                    # Get min/max and build a scoring function since JBrowse doesn't
                    if scales["type"] == "automatic" or scales["type"] == "__auto__":
                        min_val, max_val = self.min_max_gff(gff3)
                    else:
                        min_val = scales.get("min", 0)
                        max_val = scales.get("max", 1000)

                    if scheme["color"] == "__auto__":
                        user_color = "undefined"
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())
                    elif scheme["color"].startswith("#"):
                        user_color = "'%s'" % self.hex_from_rgb(
                            *self.rgb_from_hex(scheme["color"][1:])
                        )
                        auto_color = "undefined"
                    else:
                        user_color = "undefined"
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())

                    color_function = self.COLOR_FUNCTION_TEMPLATE_QUAL.format(
                        **{
                            "opacity": self.OPACITY_MATH[algo].format(
                                **{"max": max_val, "min": min_val}
                            ),
                            "user_spec_color": user_color,
                            "auto_gen_color": auto_color,
                        }
                    )

                    trackConfig["style"]["color"] = color_function.replace("\n", "")
        return trackConfig


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
        metadata["dataset_%s" % key] = value

    for key, value in node.findall("history")[0].attrib.items():
        metadata["history_%s" % key] = value

    for key, value in node.findall("metadata")[0].attrib.items():
        metadata["metadata_%s" % key] = value

    for key, value in node.findall("tool")[0].attrib.items():
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
        self.cs = ColorScaling()
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.genome_paths = genomes
        self.tracksToIndex = []

        # This is the id of the current assembly
        self.assembly_ids = {}
        self.current_assembly_id = []

        # If upgrading, look at the existing data
        self.check_existing(self.outdir)

        self.clone_jbrowse(self.jbrowse, self.outdir)

    def get_cwd(self, cwd):
        if cwd:
            return self.outdir
        else:
            return subprocess.check_output(['pwd']).decode('utf-8').strip()
            # return None

    def subprocess_check_call(self, command, output=None, cwd=False):
        if output:
            log.debug("cd %s && %s >  %s", self.get_cwd(cwd), " ".join(command), output)
            subprocess.check_call(command, cwd=self.get_cwd(cwd), stdout=output)
        else:
            log.debug("cd %s && %s", self.get_cwd(cwd), " ".join(command))
            subprocess.check_call(command, cwd=self.get_cwd(cwd))

    def subprocess_popen(self, command, cwd=False):
        log.debug("cd %s && %s", self.get_cwd(cwd), command)
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
            log.error("cd %s && %s", self.get_cwd(cwd), command)
            log.error(output)
            log.error(err)
            raise RuntimeError("Command failed with exit code %s" % (retcode))

    def subprocess_check_output(self, command, cwd=False):
        log.debug("cd %s && %s", self.get_cwd(cwd), " ".join(command))
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
            "type": "LinearBasicDisplay",
        }

        if "display" in xml_conf["style"]:
            style_data["type"] = xml_conf["style"]["display"]
            del xml_conf["style"]["display"]

        style_data["displayId"] = "%s_%s" % (xml_conf["label"], style_data["type"])

        style_data.update(xml_conf["style"])

        return {"displays": [style_data]}

    def symlink_or_copy_load_action(self):
        if "GALAXY_JBROWSE_SYMLINKS" in os.environ and bool(
            os.environ["GALAXY_JBROWSE_SYMLINKS"]
        ):
            return "symlink"
        else:
            return "copy"

    def check_existing(self, destination):
        existing = os.path.join(destination, "config.json")
        if os.path.exists(existing):
            with open(existing, "r") as existing_conf:
                conf = json.load(existing_conf)
                if "assemblies" in conf:
                    for assembly in conf["assemblies"]:
                        if "name" in assembly:
                            self.assembly_ids[assembly["name"]] = None

    def process_genomes(self):
        for genome_node in self.genome_paths:
            # We only expect one input genome per run. This for loop is just
            # easier to write than the alternative / catches any possible
            # issues.
            log.debug("Processing genome", genome_node)
            self.add_assembly(genome_node["path"], genome_node["label"])

    def add_assembly(self, path, label, default=True):
        # Find a non-existing filename for the new genome
        # (to avoid colision when upgrading an existing instance)
        rel_seq_path = os.path.join(label)
        seq_path = os.path.join(self.outdir, rel_seq_path)
        fn_try = 1
        while (
            os.path.exists(seq_path + ".fasta")
            or os.path.exists(seq_path + ".fasta.gz")
            or os.path.exists(seq_path + ".fasta.gz.fai")
            or os.path.exists(seq_path + ".fasta.gz.gzi")
        ):
            rel_seq_path = os.path.join(f"{label}{fn_try}")
            seq_path = os.path.join(self.outdir, rel_seq_path)
            fn_try += 1

        # Find a non-existing label for the new genome
        # (to avoid colision when upgrading an existing instance)
        lab_try = 1
        uniq_label = label
        while uniq_label in self.assembly_ids:
            uniq_label = label + str(lab_try)
            lab_try += 1

        # Find a default scaffold to display
        # TODO this may not be necessary in the future, see https://github.com/GMOD/jbrowse-components/issues/2708
        with open(path, "r") as fa_handle:
            fa_header = fa_handle.readline()[1:].strip().split(" ")[0]

        self.assembly_ids[uniq_label] = fa_header
        if default:
            self.current_assembly_id = uniq_label

        copied_genome = rel_seq_path + ".fasta"
        shutil.copy(path, copied_genome)

        # Compress with bgzip
        cmd = ["bgzip", copied_genome]
        self.subprocess_check_call(cmd)

        # FAI Index
        cmd = ["samtools", "faidx", copied_genome + ".gz"]
        self.subprocess_check_call(cmd)

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-assembly",
                "--load",
                self.symlink_or_copy_load_action(),
                "--name",
                uniq_label,
                "--type",
                "bgzipFasta",
                "--out",
                self.outdir,
                # "--target",
                # os.path.join(self.outdir, 'config.json'),
                "--skipCheck",
                rel_seq_path + ".fasta.gz",
            ]
        )

        return uniq_label

    def text_index(self):
        # Index tracks
        args = [
            "jbrowse",
            "text-index",
            "--target",
            self.outdir,
            "--assemblies",
            self.current_assembly_id,
        ]

        tracks = ",".join(self.tracksToIndex)
        if tracks:
            args += ["--tracks", tracks]

            self.subprocess_check_call(args)

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
        log.debug("cd %s && %s > %s", self.outdir, " ".join(cmd), gff3_unrebased.name)
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
            log.debug("cd %s && %s > %s", self.outdir, " ".join(cmd), gff3_rebased.name)
            subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_rebased)
            gff3_rebased.close()

            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)

        rel_dest = os.path.join(trackData["label"] + ".gff")
        # dest = os.path.join(self.outdir, rel_dest)

        self._sort_gff(gff3, rel_dest)
        os.unlink(gff3)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest + ".gz",
            config=style_json,
        )

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):
        rel_dest = os.path.join(trackData["label"] + ".bw")
        dest = os.path.join(self.outdir, rel_dest)
        self.symlink_or_copy(os.path.realpath(data), dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            config=style_json,
        )

    # Anything ending in "am" (Bam or Cram)
    def add_xam(self, data, trackData, xamOpts, index=None, ext="bam", **kwargs):
        index_ext = "bai"
        if ext == "cram":
            index_ext = "crai"

        rel_dest = os.path.join(trackData["label"] + ".%s" % ext)
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        if index is not None and os.path.exists(os.path.realpath(index)):
            # xai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
            self.subprocess_check_call(
                ["cp", os.path.realpath(index), dest + ".%s" % index_ext]
            )
        else:
            # Can happen in exotic condition
            # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
            #      => no index generated by galaxy, but there might be one next to the symlink target
            #      this trick allows to skip the bam sorting made by galaxy if already done outside
            if os.path.exists(os.path.realpath(data) + ".%s" % index_ext):
                self.symlink_or_copy(
                    os.path.realpath(data) + ".%s" % index_ext, dest + ".%s" % index_ext
                )
            else:
                log.warn(
                    "Could not find a bam index (.%s file) for %s", (index_ext, data)
                )

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            config=style_json,
        )

    def add_vcf(self, data, trackData, vcfOpts={}, zipped=False, **kwargs):
        if zipped:
            rel_dest = os.path.join(trackData["label"] + ".vcf.gz")
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), rel_dest)
        else:
            rel_dest = os.path.join(trackData["label"] + ".vcf")
            # dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), rel_dest)

            cmd = ["bgzip", rel_dest]
            self.subprocess_check_call(cmd)
            cmd = ["tabix", rel_dest + ".gz"]
            self.subprocess_check_call(cmd)

            rel_dest = os.path.join(trackData["label"] + ".vcf.gz")

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            config=style_json,
        )

    def add_gff(self, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join(trackData["label"] + ".gff")
        # dest = os.path.join(self.outdir, rel_dest)

        self._sort_gff(data, rel_dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest + ".gz",
            config=style_json,
        )

    def add_bed(self, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join(trackData["label"] + ".bed")
        # dest = os.path.join(self.outdir, rel_dest)

        self._sort_bed(data, rel_dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest + ".gz",
            config=style_json,
        )

    def add_paf(self, data, trackData, pafOpts, parentgenome, **kwargs):
        # TODO: how to get both parent genomes.

        # print(trackData)
        rel_dest = os.path.join(trackData["label"] + ".paf")
        # dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), rel_dest)

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
            assemblies=[trackData['key'], parentgenome],
            config=style_json,
            trackType="SyntenyTrack",
        )

    def add_hic(self, data, trackData, hicOpts, **kwargs):
        rel_dest = os.path.join(trackData["label"] + ".hic")
        # dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), rel_dest)

        style_json = self._prepare_track_style(trackData)

        self._add_track(
            trackData["label"],
            trackData["key"],
            trackData["category"],
            rel_dest,
            config=style_json,
        )

    def add_sparql(self, url, query, query_refnames, trackData):
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
            "assemblyNames": [self.current_assembly_id],
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

        # Doesn't work as of 1.6.4, might work in the future
        # self.subprocess_check_call([
        #     'jbrowse', 'add-track',
        #     '--trackType', 'sparql',
        #     '--name', trackData['label'],
        #     '--category', trackData['category'],
        #     '--target', os.path.join(self.outdir, 'data'),
        #     '--trackId', id,
        #     '--config', '{"queryTemplate": "%s"}' % query,
        #     url])

    def _add_track(self, id, label, category, path, assemblies=[], config=None, trackType=None):
        assemblies_opt = self.current_assembly_id
        if assemblies:
            assemblies_opt = ",".join(assemblies)

        cmd = [
            "jbrowse",
            "add-track",
            "--load",
            self.symlink_or_copy_load_action(),
            "--name",
            label,
            "--category",
            category,
            "--target",
            self.outdir,
            "--trackId",
            id,
            "--assemblyNames",
            assemblies_opt,
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
            cmd = "gff3sort.pl --precise '%s' | grep -v \"^$\" > '%s'" % (data, dest)
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
        _parent_genome = parent.attrib['label']
        category = track["category"].replace("__pd__date__pd__", TODAY)
        outputTrackConfig = {
            "category": category,
        }

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
                "Processing track %s / %s (%s)",
                category,
                track_human_label,
                dataset_ext,
            )
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
                self.current_assembly_id,
            ]
            hashData = "|".join(hashData).encode("utf-8")
            outputTrackConfig["label"] = hashlib.md5(hashData).hexdigest() + "_%s" % i
            outputTrackConfig["metadata"] = extra_metadata

            outputTrackConfig["style"] = track["style"]

            if "menus" in track["conf"]["options"]:
                menus = self.cs.parse_menus(track["conf"]["options"])
                outputTrackConfig.update(menus)

            print(f"Adding track {dataset_ext}")
            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                )
            elif dataset_ext == "bed":
                self.add_bed(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                    track["conf"]["options"]["gff"],
                )
            elif dataset_ext == "bigwig":
                self.add_bigwig(
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
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["cram"],
                    index=real_indexes[i],
                    ext="cram",
                )
            elif dataset_ext == "blastxml":
                self.add_blastxml(
                    dataset_path, outputTrackConfig, track["conf"]["options"]["blast"]
                )
            elif dataset_ext == "vcf":
                self.add_vcf(dataset_path, outputTrackConfig)
            elif dataset_ext == "vcf_bgzip":
                self.add_vcf(dataset_path, outputTrackConfig, zipped=True)
            elif dataset_ext == "rest":
                self.add_rest(
                    track["conf"]["options"]["rest"]["url"], outputTrackConfig
                )
            elif dataset_ext == "paf":
                log.debug("===== PAF =====")
                self.add_paf(
                    dataset_path, outputTrackConfig, track["conf"]["options"]["synteny"],
                    _parent_genome
                )
            elif dataset_ext == "hic":
                self.add_hic(
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
                    track["conf"]["options"]["sparql"]["url"],
                    sparql_query,
                    sparql_query_refnames,
                    outputTrackConfig,
                )
            else:
                log.warn("Do not know how to handle %s", dataset_ext)

            # Return non-human label for use in other fields
            yield outputTrackConfig["label"]

    def add_default_session(self, data):
        """
        Add some default session settings: set some assemblies/tracks on/off
        """
        tracks_data = []

        # TODO using the default session for now, but check out session specs in the future https://github.com/GMOD/jbrowse-components/issues/2708

        # We need to know the track type from the config.json generated just before
        config_path = os.path.join(self.outdir, "config.json")
        print(config_path)
        track_types = {}
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)

        for track_conf in config_json["tracks"]:
            track_types[track_conf["trackId"]] = track_conf["type"]

        for on_track in data["visibility"]["default_on"]:
            # TODO several problems with this currently
            # - we are forced to copy the same kind of style config as the per track config from _prepare_track_style (not exactly the same though)
            # - we get an error when refreshing the page
            # - this could be solved by session specs, see https://github.com/GMOD/jbrowse-components/issues/2708
            style_data = {"type": "LinearBasicDisplay", "height": 100}

            if on_track in data["style"]:
                if "display" in data["style"][on_track]:
                    style_data["type"] = data["style"][on_track]["display"]
                    del data["style"][on_track]["display"]

                style_data.update(data["style"][on_track])

            if on_track in data["style_labels"]:
                # TODO fix this: it should probably go in a renderer block (SvgFeatureRenderer) but still does not work
                # TODO move this to per track displays?
                style_data["labels"] = data["style_labels"][on_track]

            tracks_data.append(
                {
                    "type": track_types[on_track],
                    "configuration": on_track,
                    "displays": [style_data],
                }
            )

        # The view for the assembly we're adding
        view_json = {"type": "LinearGenomeView", "tracks": tracks_data}

        refName = None
        if data.get("defaultLocation", ""):
            loc_match = re.search(r"^(\w+):(\d+)\.+(\d+)$", data["defaultLocation"])
            if loc_match:
                refName = loc_match.group(1)
                start = int(loc_match.group(2))
                end = int(loc_match.group(3))
        elif self.assembly_ids[self.current_assembly_id] is not None:
            refName = self.assembly_ids[self.current_assembly_id]
            start = 0
            end = 1000000  # Booh, hard coded! waiting for https://github.com/GMOD/jbrowse-components/issues/2708

        if refName is not None:
            # TODO displayedRegions is not just zooming to the region, it hides the rest of the chromosome
            view_json["displayedRegions"] = [
                {
                    "refName": refName,
                    "start": start,
                    "end": end,
                    "reversed": False,
                    "assemblyName": self.current_assembly_id,
                }
            ]

        session_name = data.get("session_name", "New session")
        if not session_name:
            session_name = "New session"

        # Merge with possibly existing defaultSession (if upgrading a jbrowse instance)
        session_json = {}
        if "defaultSession" in config_json:
            session_json = config_json["defaultSession"]

        session_json["name"] = session_name

        if "views" not in session_json:
            session_json["views"] = []

        session_json["views"].append(view_json)

        config_json["defaultSession"] = session_json

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
            log.error("Error: %s - %s." % (e.filename, e.strerror))

        if not os.path.exists(os.path.join(destination, "data")):
            # It can already exist if upgrading an instance
            os.makedirs(os.path.join(destination, "data"))
            log.info("makedir %s" % (os.path.join(destination, "data")))

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
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument("--version", "-V", action="version", version="%(prog)s 0.8.0")
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
        genomes=[
            {
                "path": os.path.realpath(x.attrib["path"]),
                "meta": metadata_from_node(x.find("metadata")),
                "label": x.attrib["label"],
            }
            for x in real_root.findall("assembly/genomes/genome")
        ],
    )

    # Process genomes
    jc.process_genomes()

    default_session_data = {
        "visibility": {
            "default_on": [],
            "default_off": [],
        },
        "style": {},
        "style_labels": {},
    }

    for assembly in real_root.findall("assembly"):
        genome = assembly.find('genomes/genome')

        # TODO add metadata to tracks
        for track in assembly.findall("tracks/track"):
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

            track_conf["conf"] = etree_to_dict(track.find("options"))

            keys = jc.process_annotations(track_conf, genome)
            keys = list(keys)

            # for key in keys:
            #     default_session_data["visibility"][
            #         track.attrib.get("visibility", "default_off")
            #     ].append(key)
            #
            # default_session_data["style"][key] = track_conf[
            #     "style"
            # ]  # TODO do we need this anymore?
            # default_session_data["style_labels"][key] = track_conf["style_labels"]

    default_session_data["defaultLocation"] = real_root.find(
        "metadata/general/defaultLocation"
    ).text
    default_session_data["session_name"] = real_root.find(
        "metadata/general/session_name"
    ).text

    general_data = {
        "analytics": real_root.find("metadata/general/analytics").text,
        "primary_color": real_root.find("metadata/general/primary_color").text,
        "secondary_color": real_root.find("metadata/general/secondary_color").text,
        "tertiary_color": real_root.find("metadata/general/tertiary_color").text,
        "quaternary_color": real_root.find("metadata/general/quaternary_color").text,
        "font_size": real_root.find("metadata/general/font_size").text,
    }

    jc.add_default_session(default_session_data)
    jc.add_general_configuration(general_data)
    jc.text_index()
