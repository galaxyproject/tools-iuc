#!/usr/bin/env python

import argparse
import binascii
import datetime
import json
import logging
import hashlib
import os
import re
import shutil
import ssl
import struct
import subprocess
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("jbrowse")

JB2VER = "v2.11.0"
# version pinned if cloning - but not used until now
logCommands = True
# useful for seeing what's being written but not for production setups
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


INDEX_TEMPLATE = """<!doctype html>
<html lang="en" style="height:100%">
<head>
<meta charset="utf-8"/>
<link rel="shortcut icon" href="./favicon.ico"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<meta name="theme-color" content="#000000"/>
<meta name="description" content="A fast and flexible genome browser"/>
<link rel="manifest" href="./manifest.json"/>
<title>JBrowse</title>
</script>
</head>
<body style="overscroll-behavior:none; height:100%; margin: 0;">
<iframe
  id="jbframe"
  title="JBrowse2"
  frameborder="0"
  width="100%"
  height="100%"
  src='index_noview.html?config=config.json__SESSION_SPEC__'>
</iframe>
</body>
</html>
"""


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

    for (key, value) in node.findall("dataset")[0].attrib.items():
        metadata["dataset_%s" % key] = value

    if node.findall("history"):
        for (key, value) in node.findall("history")[0].attrib.items():
            metadata["history_%s" % key] = value

    if node.findall("metadata"):
        for (key, value) in node.findall("metadata")[0].attrib.items():
            metadata["metadata_%s" % key] = value
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
            encoded_hist_id=metadata.get("history_id", "not available"),
            hist_name=metadata.get("history_display_name", "not available"),
        )
    if node.findall("tool"):
        for (key, value) in node.findall("tool")[0].attrib.items():
            metadata["tool_%s" % key] = value
        metadata[
            "tool_tool"
        ] = '<a target="_blank" href="{galaxy}/datasets/{encoded_id}/show_params">{tool_id}{tool_version}</a>'.format(
            galaxy=GALAXY_INFRASTRUCTURE_URL,
            encoded_id=metadata.get("dataset_id", ""),
            tool_id=metadata.get("tool_tool_id", ""),
            tool_version=metadata.get("tool_tool_version", ""),
        )
    return metadata


class JbrowseConnector(object):
    def __init__(self, outdir, jbrowse2path):
        self.trackCounter = 0  # to avoid name clashes
        self.assemblies = []  # these require more than a few line diff.
        self.assmeta = {}
        self.ass_first_contigs = (
            []
        )  # for default session - these are read as first line of the assembly .fai
        self.giURL = GALAXY_INFRASTRUCTURE_URL
        self.outdir = outdir
        self.jbrowse2path = jbrowse2path
        os.makedirs(self.outdir, exist_ok=True)
        self.genome_names = []
        self.trackIdlist = []
        self.tracksToAdd = {}
        self.config_json = {}
        self.config_json_file = os.path.join(outdir, "config.json")
        self.clone_jbrowse()

    def get_cwd(self, cwd):
        if cwd:
            return self.outdir
        else:
            return subprocess.check_output(["pwd"]).decode("utf-8").strip()

    def subprocess_check_call(self, command, output=None, cwd=True):
        if output:
            if logCommands:
                log.debug(
                    "cd %s && %s >  %s", self.get_cwd(cwd), " ".join(command), output
                )
            subprocess.check_call(command, cwd=self.get_cwd(cwd), stdout=output)
        else:
            if logCommands:
                log.debug("cd %s && %s", self.get_cwd(cwd), " ".join(command))
            subprocess.check_call(command, cwd=self.get_cwd(cwd))

    def subprocess_popen(self, command, cwd=True):
        if logCommands:
            log.debug(command)
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
            log.error(command)
            log.error(output)
            log.error(err)
            raise RuntimeError(f"Command ( {command} ) failed with exit code {retcode}")

    def subprocess_check_output(self, command):
        if logCommands:
            log.debug(" ".join(command))
        return subprocess.check_output(command, cwd=self.outdir)

    def symlink_or_copy(self, src, dest):
        if "GALAXY_JBROWSE_SYMLINKS" in os.environ and bool(
            os.environ["GALAXY_JBROWSE_SYMLINKS"]
        ):
            cmd = ["ln", "-s", src, dest]
        else:
            cmd = ["cp", src, dest]

        return self.subprocess_check_call(cmd)

    def _prepare_track_style(self, trackDict):

        style_data = {
            "type": "LinearBasicDisplay",
            "displayId": "%s-LinearBasicDisplay" % trackDict["trackId"],
        }

        if trackDict.get("displays", None):  # use first if multiple like bed
            style_data["type"] = trackDict["displays"][0]["type"]
            style_data["displayId"] = trackDict["displays"][0]["displayId"]
        return style_data

    def getNrow(self, url):
        useuri = url.startswith("https://") or url.startswith("http://")
        if not useuri:
            fl = open(url, "r").readlines()
            nrow = len(fl)
        else:
            try:
                scontext = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
                scontext.check_hostname = False
                scontext.verify_mode = ssl.VerifyMode.CERT_NONE
                with urllib.request.urlopen(url, context=scontext) as f:
                    fl = f.readlines()
                nrow = len(fl)
            except Exception:
                nrow = 0
        logging.debug("getNrow %s returning %d" % (url, nrow))
        return nrow

    def process_genomes(self, genomes):
        assembly = []
        assmeta = []
        useuri = False
        primaryGenome = None
        for i, genome_node in enumerate(genomes):
            this_genome = {}
            if genome_node["useuri"] == "yes":
                useuri = True
            genome_name = genome_node["label"].strip()
            if len(genome_name) == 0:
                genome_name = os.path.splitext(os.path.basename(genome_node["path"]))[0]
            if len(genome_name.split()) > 1:
                genome_name = genome_name.split()[0]
                # spaces and cruft break scripts when substituted
            if not primaryGenome:
                primaryGenome = genome_name
            if genome_name not in self.genome_names:
                self.genome_names.append(genome_name)
                fapath = genome_node["path"]
                if not useuri:
                    fapath = os.path.realpath(fapath)
                assem, first_contig = self.make_assembly(fapath, genome_name, useuri)
                assembly.append(assem)
                self.ass_first_contigs.append(first_contig)
                if genome_name == primaryGenome:  # first one
                    this_genome["genome_name"] = genome_name  # first one for all tracks
                    this_genome["genome_sequence_adapter"] = assem["sequence"][
                        "adapter"
                    ]
                    this_genome["genome_firstcontig"] = first_contig
                assmeta.append(this_genome)
        self.assemblies += assembly
        self.assmeta[primaryGenome] = assmeta
        self.tracksToAdd[primaryGenome] = []
        return primaryGenome

    def make_assembly(self, fapath, gname, useuri):
        if useuri:
            faname = fapath
            scontext = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
            scontext.check_hostname = False
            scontext.verify_mode = ssl.VerifyMode.CERT_NONE
            with urllib.request.urlopen(url=faname + ".fai", context=scontext) as f:
                fl = f.readline()
            contig = fl.decode("utf8").strip()
            # Merlin  172788  8       60      61
        else:
            faname = gname + ".fa.gz"
            fadest = os.path.realpath(os.path.join(self.outdir, faname))
            cmd = "bgzip -i -c %s -I %s.gzi > %s && samtools faidx %s" % (
                fapath,
                fadest,
                fadest,
                fadest,
            )
            self.subprocess_popen(cmd)
            contig = open(fadest + ".fai", "r").readline().strip()
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
        first_contig = contig.split()[:2]
        first_contig.insert(0, gname)
        trackDict = {
            "name": gname,
            "sequence": {
                "type": "ReferenceSequenceTrack",
                "trackId": gname,
                "adapter": adapter,
            },
            "displays": [
                {
                    "type": "LinearReferenceSequenceDisplay",
                    "displayId": "%s-LinearReferenceSequenceDisplay" % gname,
                },
                {
                    "type": "LinearGCContentDisplay",
                    "displayId": "%s-LinearGCContentDisplay" % gname,
                },
            ],
        }
        return (trackDict, first_contig)

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
        self.subprocess_check_call(cmd)

    def write_config(self):
        with open(self.config_json_file, "w") as fp:
            json.dump(self.config_json, fp, indent=2)

    def text_index(self):
        # Index tracks
        args = [
            "jbrowse",
            "text-index",
            "--target",
            self.outdir,
            "--assemblies",
            self.genome_name,
        ]

        tracks = ",".join(self.trackIdlist)
        if tracks:
            args += ["--tracks", tracks]

            self.subprocess_check_call(args)

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
        tId = trackData["label"]
        wasCool = trackData["wasCool"]
        # can be served - if public.
        # dsId = trackData["metadata"]["dataset_id"]
        # url = "%s/api/datasets/%s/display?to_ext=hic " % (self.giURL, dsId)
        useuri = trackData["useuri"].lower() == "yes"
        logging.debug("wasCool=%s, data=%s, tId=%s" % (wasCool, data, tId))
        if useuri:
            uri = data
        else:
            uri = tId + ".hic"
            if not wasCool:
                dest = os.path.join(self.outdir, uri)
                if not os.path.exists(dest):
                    cmd = ["cp", data, dest]
                    self.subprocess_check_call(cmd)
                else:
                    logging.error("not wasCool but %s exists" % dest)
        categ = trackData["category"]
        trackDict = {
            "type": "HicTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [trackData["assemblyNames"]],
            "category": [
                categ,
            ],
            "adapter": {"type": "HicAdapter", "hicLocation": {"uri": uri}},
        }
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_maf(self, data, trackData):
        """
        from https://github.com/cmdcolin/maf2bed
        Note: Both formats start with a MAF as input, and note that your MAF file should contain the species name and chromosome name
        e.g. hg38.chr1 in the sequence identifiers.
        need the reference id - eg hg18, for maf2bed.pl as the first parameter
        """
        tId = trackData["label"]
        mafPlugin = {
            "plugins": [
                {
                    "name": "MafViewer",
                    "url": "https://unpkg.com/jbrowse-plugin-mafviewer/dist/jbrowse-plugin-mafviewer.umd.production.min.js",
                }
            ]
        }
        categ = trackData["category"]
        fname = f"{tId}"
        dest = os.path.join(self.outdir, fname)
        gname = trackData["assemblyNames"]

        cmd = [
            "bash",
            os.path.join(INSTALLED_TO, "convertMAF.sh"),
            data,
            gname,
            INSTALLED_TO,
            dest,
        ]
        self.subprocess_check_call(cmd)
        mafs = open(data, "r").readlines()
        mafss = [x for x in mafs if (x.startswith("s\t") or x.startswith("s "))]
        samp = [x.split()[1] for x in mafss if len(x.split()) > 0]
        sampu = list(dict.fromkeys(samp))
        samples = [x.split(".")[0] for x in sampu]
        samples.sort()
        if logCommands:
            logging.debug(
                "$$$$ cmd=%s, mafss=%s samp=%s samples=%s"
                % (" ".join(cmd), mafss, samp, samples)
            )
        trackDict = {
            "type": "MafTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
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
            "assemblyNames": [trackData["assemblyNames"]],
            "displays": [
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearBasicDisplay" % tId,
                },
                {"type": "LinearArcDisplay", "displayId": "%s-LinearArcDisplay" % tId},
            ],
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[gname].append(trackDict)
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
        logging.debug("### blastxml to gff3 cmd = %s" % " ".join(cmd))
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
            logging.debug("### gff3rebase cmd = %s" % " ".join(cmd))
            gff3_rebased.close()
            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)
        self.add_gff(gff3, trackData, **kwargs)

    def add_bigwig(self, data, trackData):
        tId = trackData["label"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            url = tId
            # slashes in names cause path trouble
            dest = os.path.join(self.outdir, url)
            cmd = ["cp", data, dest]
            self.subprocess_check_call(cmd)
        bwloc = {"uri": url}
        categ = trackData["category"]
        trackDict = {
            "type": "QuantitativeTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
            "assemblyNames": [trackData["assemblyNames"]],
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
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_bam(self, data, trackData, bam_indexes=None, **kwargs):
        tId = trackData["label"]
        realFName = trackData["path"]
        useuri = trackData["useuri"].lower() == "yes"
        categ = trackData["category"]
        if useuri:
            url = data
        else:
            fname = tId
            dest = "%s/%s" % (self.outdir, fname)
            self.subprocess_check_call(["cp", data, dest])
            url = fname
            bindex = fname + ".bai"
            bi = bam_indexes.split(",")
            bam_index = [
                x.split(" ~ ")[1].strip()
                for x in bi
                if " ~ " in x and x.split(" ~ ")[0].strip() == realFName
            ]
            logging.debug(
                "===realFName=%s got %s as bam_indexes %s as bi, %s for bam_index"
                % (realFName, bam_indexes, bi, bam_index)
            )
            if len(bam_index) > 0 and os.path.exists(os.path.realpath(bam_index[0])):
                self.subprocess_check_call(["cp", bam_index[0], bindex])
            else:
                cmd = ["samtools", "index", "-b", "-o", bindex, data]
                self.subprocess_check_call(cmd)
        trackDict = {
            "type": "AlignmentsTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
            "assemblyNames": [trackData["assemblyNames"]],
            "adapter": {
                "type": "BamAdapter",
                "bamLocation": {"uri": url},
                "index": {
                    "location": {
                        "uri": bindex,
                    }
                },
            },
            "displays": [
                {
                    "type": "LinearAlignmentsDisplay",
                    "displayId": "%s-LinearAlignmentsDisplay" % tId,
                },
            ],
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_cram(self, data, trackData, cram_indexes=None, **kwargs):
        tId = trackData["label"]
        realFName = trackData["path"]
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        gsa = self.assmeta.get(trackData["assemblyNames"], None)
        if gsa:
            genseqad = gsa[0]["genome_sequence_adapter"]
        else:
            genseqad = "Not found"
            logging.warning("No adapter found for cram %s in gsa=%s" % (tId, gsa))
        if useuri:
            url = data
        else:
            fname = tId
            dest = os.path.join(self.outdir, fname)
            url = fname
            self.subprocess_check_call(["cp", data, dest])
            ci = cram_indexes.split(",")
            cram_index = [
                x.split(" ~ ")[1].strip()
                for x in ci
                if " ~ " in x and x.split(" ~ ")[0].strip() == realFName
            ]
            logging.debug(
                "===realFName=%s got %s as cram_indexes %s as ci, %s for cram_index"
                % (realFName, cram_indexes, ci, cram_index)
            )
            if len(cram_index) > 0 and os.path.exists(cram_index[0]):
                if not os.path.exists(dest + ".crai"):
                    # most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
                    self.subprocess_check_call(
                        ["cp", os.path.realpath(cram_index[0]), dest + ".crai"]
                    )
            else:
                cpath = os.path.realpath(dest) + ".crai"
                cmd = ["samtools", "index", "-c", "-o", cpath, os.path.realpath(dest)]
                self.subprocess_check_call(cmd)
        trackDict = {
            "type": "AlignmentsTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
            "assemblyNames": [trackData["assemblyNames"]],
            "adapter": {
                "type": "CramAdapter",
                "cramLocation": {"uri": url},
                "craiLocation": {
                    "uri": url + ".crai",
                },
                "sequenceAdapter": genseqad,
            },
            "displays": [
                {
                    "type": "LinearAlignmentsDisplay",
                    "displayId": "%s-LinearAlignmentsDisplay" % tId,
                },
            ],
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_vcf(self, data, trackData):
        tId = trackData["label"]
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            url = tId
            dest = os.path.join(self.outdir, url)
            cmd = "bgzip -c %s  > %s" % (data, dest)
            self.subprocess_popen(cmd)
            cmd = ["tabix", "-f", "-p", "vcf", dest]
            self.subprocess_check_call(cmd)
        trackDict = {
            "type": "VariantTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [trackData["assemblyNames"]],
            "category": [
                categ,
            ],
            "adapter": {
                "type": "VcfTabixAdapter",
                "vcfGzLocation": {"uri": url},
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
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "jbrowse sort-gff '%s' | bgzip -c > '%s'" % (
                data,
                dest,
            )
            self.subprocess_popen(cmd)
            self.subprocess_check_call(["tabix", "-f", "-p", "gff", dest])

    def _sort_bed(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "sort -k1,1 -k2,2n '%s' | bgzip -c > '%s'" % (data, dest)
            self.subprocess_popen(cmd)
            cmd = ["tabix", "-f", "-p", "bed", dest]
            self.subprocess_check_call(cmd)

    def add_gff(self, data, trackData):
        tId = trackData["label"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = trackData["path"]
        else:
            url = tId + ".gz"
            dest = os.path.join(self.outdir, url)
            self._sort_gff(data, dest)
        categ = trackData["category"]
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [trackData["assemblyNames"]],
            "category": [
                categ,
            ],
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
                {
                    "type": "LinearArcDisplay",
                    "displayId": "%s-LinearArcDisplay" % tId,
                },
            ],
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_bed(self, data, ext, trackData):
        tId = trackData["label"]
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            url = tId + ".gz"
            dest = os.path.join(self.outdir, url)
            self._sort_bed(data, dest)
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [trackData["assemblyNames"]],
            "adapter": {
                "category": [
                    categ,
                ],
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
                {
                    "type": "LinearPileupDisplay",
                    "displayId": "%s-LinearPileupDisplay" % tId,
                },
                {
                    "type": "LinearArcDisplay",
                    "displayId": "%s-LinearArcDisplay" % tId,
                },
            ],
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def add_paf(self, data, trackData, pafOpts, **kwargs):
        tname = trackData["name"]
        tId = trackData["label"]
        url = tId
        useuri = data.startswith("http://") or data.startswith("https://")
        if not useuri:
            dest = os.path.join(self.outdir, url)
            self.symlink_or_copy(os.path.realpath(data), dest)
            nrow = self.getNrow(dest)
        else:
            url = data
            nrow = self.getNrow(url)
        categ = trackData["category"]
        pg = pafOpts["genome"].split(",")
        pgc = [x.strip() for x in pg if x.strip() > ""]
        gnomes = [x.split(" ~ ") for x in pgc]
        logging.debug("pg=%s, gnomes=%s" % (pg, gnomes))
        passnames = [trackData["assemblyNames"]]  # always first
        for i, (gpath, gname) in enumerate(gnomes):
            # may have been forgotten by user for uri
            if len(gname) == 0:
                gn = os.path.basename(gpath)
                gname = os.path.splitext(gn)[0]
            # trouble from spacey names in command lines avoidance
            if len(gname.split()) > 1:
                gname = gname.split()[0]
            if gname not in passnames:
                passnames.append(gname)
            useuri = pafOpts["useuri"] == "true"
            if gname not in self.genome_names:
                # ignore if already there - eg for duplicates among pafs.
                asstrack, first_contig = self.make_assembly(gpath, gname, useuri)
                self.genome_names.append(gname)
                self.tracksToAdd[gname] = []
                self.assemblies.append(asstrack)
                self.ass_first_contigs.append(first_contig)
        trackDict = {
            "type": "SyntenyTrack",
            "trackId": tId,
            "assemblyNames": passnames,
            "category": [
                categ,
            ],
            "name": tname,
            "adapter": {
                "type": "PAFAdapter",
                "pafLocation": {"uri": url},
                "assemblyNames": passnames,
            },
            "displays": [
                {
                    "type": "LGVSyntenyDisplay",
                    "displayId": "%s-LGVSyntenyDisplay" % tId,
                },
                {
                    "type": "DotplotDisplay",
                    "displayId": "%s-DotplotDisplay" % tId,
                },
                {
                    "type": "LinearComparativeDisplay",
                    "displayId": "%s-LinearComparativeDisplay" % tId,
                },
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearSyntenyDisplay" % tId,
                },
            ],
        }
        if nrow > 10000:
            style_json = {
                "type": "LGVSyntenyDisplay",
                "displayId": "%s-LGVSyntenyDisplay" % tId,
            }
        else:
            style_json = {
                "type": "LinearBasicDisplay",
                "displayId": "%s-LinearBasicDisplay" % tId,
            }
        trackDict["style"] = style_json
        self.tracksToAdd[trackData["assemblyNames"]].append(trackDict)
        self.trackIdlist.append(tId)

    def process_annotations(self, track):
        category = track["category"].replace("__pd__date__pd__", TODAY)
        for trackIndex, (
            dataset_path,
            dataset_ext,
            useuri,
            track_human_label,
            extra_metadata,
        ) in enumerate(track["trackfiles"]):
            if not dataset_path.strip().startswith("http"):
                # Unsanitize labels (element_identifiers are always sanitized by Galaxy)
                for key, value in mapped_chars.items():
                    track_human_label = track_human_label.replace(value, key)
                track_human_label = track_human_label.replace(" ", "_")
            outputTrackConfig = {
                "category": category,
                "style": {},
            }

            hashData = [
                str(dataset_path),
                track_human_label,
                track["category"],
            ]
            hashData = "|".join(hashData).encode("utf-8")
            hash_string = hashlib.md5(hashData).hexdigest()

            outputTrackConfig["assemblyNames"] = track["assemblyNames"]
            outputTrackConfig["key"] = track_human_label
            outputTrackConfig["useuri"] = useuri
            outputTrackConfig["path"] = dataset_path
            outputTrackConfig["ext"] = dataset_ext
            outputTrackConfig["trackset"] = track.get("trackset", {})
            outputTrackConfig["label"] = track["label"]
            #outputTrackConfig["label"] = "%s_%i_%s_%s" % (
            #    dataset_ext,
            #    trackIndex,
            #    track_human_label,
            #    hash_string,
            #)

            outputTrackConfig["metadata"] = extra_metadata
            outputTrackConfig["name"] = track_human_label
            if track["label"] in self.trackIdlist:
                logging.error(
                    "### not adding %s already in %s"
                    % (track["label"], self.trackIdlist)
                )
                yield None
            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext in ("hic", "juicebox_hic"):
                outputTrackConfig["wasCool"] = False
                self.add_hic(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext in ("cool", "mcool", "scool"):
                hic_url = outputTrackConfig["label"]
                hic_path = os.path.join(self.outdir, hic_url) + ".hic"
                outputTrackConfig["wasCool"] = True
                self.subprocess_check_call(
                    [
                        "hictk",
                        "convert",
                        "-f",
                        "--output-fmt",
                        "hic",
                        dataset_path,
                        hic_path,
                    ]
                )
                self.add_hic(
                    hic_path,
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
                real_indexes = track["conf"]["options"]["bam"]["bam_index"]
                self.add_bam(
                    dataset_path,
                    outputTrackConfig,
                    bam_indexes=real_indexes,
                )
            elif dataset_ext == "cram":
                real_indexes = track["conf"]["options"]["cram"]["cram_index"]
                self.add_cram(
                    dataset_path,
                    outputTrackConfig,
                    cram_indexes=real_indexes,
                )
            elif dataset_ext == "blastxml":
                self.add_blastxml(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["blast"],
                )
            elif dataset_ext == "vcf":
                self.add_vcf(dataset_path, outputTrackConfig)
            elif dataset_ext == "paf":
                self.add_paf(
                    dataset_path,
                    outputTrackConfig,
                    track["conf"]["options"]["paf"],
                )
            else:
                logging.warning("Do not know how to handle %s", dataset_ext)
            # Return non-human label for use in other fields
            yield outputTrackConfig["label"]

    def add_default_session(self, default_data):
        """
        default session settings are hard and fragile.
        .add_default_view() and other configuration code adapted from
         https://github.com/abretaud/tools-iuc/blob/jbrowse2/tools/jbrowse2/jbrowse2.py
        """
        # TODO using the default session for now, but check out session specs in the future https://github.com/GMOD/jbrowse-components/issues/2708
        track_types = {}
        with open(self.config_json_file, "r") as config_file:
            config_json = json.load(config_file)
        if self.config_json:
            config_json.update(self.config_json)
        if "defaultSession" in config_json:
            session_json = config_json["defaultSession"]
            session_views = []
        else:
            session_json = {}
            session_views = []
        for gnome in self.assmeta.keys():  # assemblies have their own tracks
            tracks_data = []
            for track_conf in self.tracksToAdd[gnome]:
                tId = track_conf["trackId"]
                if tId in default_data[gnome]["visibility"]["default_on"]:
                    track_types[tId] = track_conf["type"]
                    style_data = default_data[gnome]["style"].get(tId, None)
                    if not style_data:
                        logging.debug(
                            "### No style data for %s in available default data %s"
                            % (tId, default_data)
                        )
                        style_data = {"type": "LinearBasicDisplay"}
                        if "displays" in track_conf:
                            disp = track_conf["displays"][0]["type"]
                            style_data["type"] = disp
                    if track_conf.get("style_labels", None):
                        # TODO fix this: it should probably go in a renderer block (SvgFeatureRenderer) but still does not work
                        # TODO move this to per track displays?
                        style_data["labels"] = track_conf["style_labels"]
                    tracks_data.append(
                        {
                            "type": track_types[tId],
                            "configuration": tId,
                            "displays": [style_data],
                        }
                    )
            view_json = {
                "type": "LinearGenomeView",
                "offsetPx": 0,
                "minimized": False,
                "tracks": tracks_data,
            }
            first = [x for x in self.ass_first_contigs if x[0] == gnome]
            if len(first) > 0:
                [gnome, refName, end] = first[0]
                start = 0
                end = int(end)
                drdict = {
                    "refName": refName,
                    "start": start,
                    "end": end,
                    "reversed": False,
                    "assemblyName": gnome,
                }
            else:
                ddl = default_data.get("defaultLocation", None)
                if ddl:
                    loc_match = re.search(r"^([^:]+):([\d,]*)\.*([\d,]*)$", ddl)
                    # allow commas like 100,000 but ignore as integer
                    if loc_match:
                        refName = loc_match.group(1)
                        drdict["refName"] = refName
                        if loc_match.group(2) > "":
                            drdict["start"] = int(loc_match.group(2).replace(",", ""))
                        if loc_match.group(3) > "":
                            drdict["end"] = int(loc_match.group(3).replace(",", ""))
                    else:
                        logging.info(
                            "@@@ regexp could not match contig:start..end in the supplied location %s - please fix"
                            % ddl
                        )
            if drdict.get("refName", None):
                # TODO displayedRegions is not just zooming to the region, it hides the rest of the chromosome
                view_json["displayedRegions"] = [
                    drdict,
                ]
                logging.info("@@@ defaultlocation %s for default session" % drdict)
            else:
                logging.info(
                    "@@@ no track location for default session - please add one!"
                )
            session_views.append(view_json)
        session_name = default_data.get("session_name", "New session")
        for key, value in mapped_chars.items():
            session_name = session_name.replace(value, key)
        session_json["name"] = session_name

        if "views" not in session_json:
            session_json["views"] = session_views
        else:
            session_json["views"] += session_views

        pp = json.dumps(session_views, indent=2)
        config_json["defaultSession"] = session_json
        self.config_json.update(config_json)
        logging.debug("defaultSession=%s" % (pp))
        with open(self.config_json_file, "w") as config_file:
            json.dump(self.config_json, config_file, indent=2)

    def add_defsess_to_index(self, data):
        """
        ----------------------------------------------------------
        Add some default session settings: set some assemblies/tracks on/off

        This allows to select a default view:
        - jb type (Linear, Circular, etc)
        - default location on an assembly
        - default tracks
        - ...

        Different methods to do that were tested/discussed:
        - using a defaultSession item in config.json: this proved to be difficult:
          forced to write a full session block, including hard-coded/hard-to-guess items,
          no good way to let Jbrowse2 display a scaffold without knowing its size
        - using JBrowse2 as an embedded React component in a tool-generated html file:
          it works but it requires generating js code to actually do what we want = chosing default view, assembly, tracks, ...
        - writing a session-spec inside the config.json file: this is not yet supported as of 2.10.2 (see PR 4148 below)
          a session-spec is a kind of simplified defaultSession where you don't need to specify every aspect of the session
        - passing a session-spec through URL params by embedding the JBrowse2 index.html inside an iframe
          we selected this option

        Xrefs to understand the choices:
        https://github.com/GMOD/jbrowse-components/issues/2708
        https://github.com/GMOD/jbrowse-components/discussions/3568
        https://github.com/GMOD/jbrowse-components/pull/4148
        """
        new_index = "Nothing written"
        session_spec = {"views": []}
        logging.debug("def ass_first=%s\ndata=%s" % (self.ass_first_contigs, data))
        for first_contig in self.ass_first_contigs:
            logging.debug("first contig=%s" % self.ass_first_contigs)
            [gnome, refName, end] = first_contig
            start = 0
            aview = {
                "assembly": gnome,
                "loc": "{}:{}..{}".format(refName, start, end),
                "type": "LinearGenomeView",
                "tracks": data[gnome]["tracks"],
            }
            session_spec["views"].append(aview)
        sess = json.dumps(session_spec, sort_keys=True, indent=2)
        new_index = INDEX_TEMPLATE.replace(
            "__SESSION_SPEC__", "&session=spec-{}".format(sess)
        )

        os.rename(
            os.path.join(self.outdir, "index.html"),
            os.path.join(self.outdir, "index_noview.html"),
        )

        with open(os.path.join(self.outdir, "index.html"), "w") as nind:
            nind.write(new_index)
        logging.debug(
            "#### add_defsession gnome=%s refname=%s\nsession_spec=%s\nnew_index=%s"
            % (gnome, refName, sess, new_index)
        )

    def add_general_configuration(self, data):
        """
        Add some general configuration to the config.json file
        """

        config_path = self.config_json_file
        if os.path.exists(config_path):
            with open(config_path, "r") as config_file:
                config_json = json.load(config_file)
        else:
            config_json = {}
        if self.config_json:
            config_json.update(self.config_json)
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
        if not config_json.get("configuration", None):
            config_json["configuration"] = {}
        config_json["configuration"].update(config_data)
        self.config_json.update(config_json)
        with open(config_path, "w") as config_file:
            json.dump(self.config_json, config_file, indent=2)

    def clone_jbrowse(self, realclone=False):
        """
            Clone a JBrowse directory into a destination directory.

            `realclone=true` will use the `jbrowse create` command.
            To allow running on internet-less compute and for reproducibility
            use frozen code with `realclone=false

        """
        dest = self.outdir
        if realclone:
            self.subprocess_check_call(
                ["jbrowse", "create", dest, "-f", "--tag", f"{JB2VER}"]
            )
        else:
            shutil.copytree(self.jbrowse2path, dest, dirs_exist_ok=True)
        for fn in [
            "asset-manifest.json",
            "favicon.ico",
            "robots.txt",
            "umd_plugin.js",
            "version.txt",
            "test_data",
        ]:
            try:
                path = os.path.join(dest, fn)
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
            except OSError as e:
                log.error("Error: %s - %s." % (e.filename, e.strerror))
        shutil.copyfile(os.path.join(INSTALLED_TO, "jb2_webserver.py"), os.path.join(dest, "jb2_webserver.py"))


def parse_style_conf(item):
    if item.text.lower() in ["false", "true", "yes", "no"]:
        return item.text.lower in ("yes", "true")
    else:
        return item.text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("--xml", help="Track Configuration")
    parser.add_argument(
        "--jbrowse2path", help="Path to JBrowse2 directory in BioContainer or Conda"
    )
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument("--version", "-V", action="version", version=JB2VER)
    args = parser.parse_args()
    tree = ET.parse(args.xml)
    root = tree.getroot()

    # This should be done ASAP
    GALAXY_INFRASTRUCTURE_URL = root.find("metadata/galaxyUrl").text
    # Sometimes this comes as `localhost` without a protocol
    if not GALAXY_INFRASTRUCTURE_URL.startswith("http"):
        # so we'll prepend `http://` and hope for the best. Requests *should*
        # be GET and not POST so it should redirect OK
        GALAXY_INFRASTRUCTURE_URL = "http://" + GALAXY_INFRASTRUCTURE_URL

    jc = JbrowseConnector(outdir=args.outdir, jbrowse2path=args.jbrowse2path)

    default_session_data = {}
    trackI = 0
    for ass in root.findall("assembly"):
        genomes = [
            {
                "path": x.attrib["path"],
                "label": x.attrib["label"].split(" ")[0].replace(",", ""),
                "useuri": x.attrib["useuri"],
                "meta": metadata_from_node(x.find("metadata")),
            }
            for x in ass.findall("metadata/genomes/genome")
        ]
        primaryGenome = jc.process_genomes(genomes)
        if not default_session_data.get(primaryGenome, None):
            default_session_data[primaryGenome] = {
                "tracks": [],
                "style": {},
                "style_labels": {},
                "visibility": {
                    "default_on": [],
                    "default_off": [],
                },
            }
        for track in ass.find("tracks"):
            track_conf = {}
            track_conf["trackfiles"] = []
            track_conf["assemblyNames"] = primaryGenome
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
                for x in trackfiles:
                    track_conf["label"] = "%s_%d" % (
                        x.attrib["label"].replace(" ", "_").replace(",", ""),
                        trackI,
                    )
                    trackI += 1
                    track_conf["useuri"] = x.attrib["useuri"]
                    if is_multi_bigwig:
                        multi_bigwig_paths.append(
                            (
                                track_conf["label"],
                                track_conf["useuri"],
                                os.path.realpath(x.attrib["path"]),
                            )
                        )
                    else:
                        if trackfiles:
                            metadata = metadata_from_node(x.find("metadata"))
                            track_conf["dataset_id"] = metadata.get(
                                "dataset_id", "None"
                            )
                            if x.attrib["useuri"].lower() == "yes":
                                tfa = (
                                    x.attrib["path"],
                                    x.attrib["ext"],
                                    x.attrib["useuri"],
                                    track_conf["label"],
                                    metadata,
                                )
                            else:
                                tfa = (
                                    os.path.realpath(x.attrib["path"]),
                                    x.attrib["ext"],
                                    x.attrib["useuri"],
                                    track_conf["label"],
                                    metadata,
                                )
                            track_conf["trackfiles"].append(tfa)

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
            track_conf["conf"] = etree_to_dict(track.find("options"))
            keys = jc.process_annotations(track_conf)
            if keys:
                for key in keys:
                    vis = track.attrib.get("visibility", "default_off")
                    if not vis:
                        vis = "default_off"
                    default_session_data[primaryGenome]["visibility"][vis].append(key)
                    trakdat = jc.tracksToAdd[primaryGenome]
                    stile = {}
                    for trak in trakdat:
                        if trak["trackId"] == key:
                            stile = trak.get("style", {})
                    if track.find("options/style"):
                        supdate = {
                            item.tag: parse_style_conf(item)
                            for item in track.find("options/style")
                        }
                        stile.update(supdate)
                    default_session_data[primaryGenome]["style"][key] = stile
                    if track.find("options/style_labels"):
                        default_session_data[primaryGenome]["style_labels"][key] = {
                            item.tag: parse_style_conf(item)
                            for item in track.find("options/style_labels")
                        }
                    default_session_data[primaryGenome]["tracks"].append(key)
    default_session_data["defaultLocation"] = root.find(
        "metadata/general/defaultLocation"
    ).text
    default_session_data["session_name"] = root.find(
        "metadata/general/session_name"
    ).text
    logging.debug("default_session=%s" % (json.dumps(default_session_data, indent=2)))
    jc.zipOut = root.find("metadata/general/zipOut").text == "true"
    general_data = {
        "analytics": root.find("metadata/general/analytics").text,
        "primary_color": root.find("metadata/general/primary_color").text,
        "secondary_color": root.find("metadata/general/secondary_color").text,
        "tertiary_color": root.find("metadata/general/tertiary_color").text,
        "quaternary_color": root.find("metadata/general/quaternary_color").text,
        "font_size": root.find("metadata/general/font_size").text,
    }
    jc.add_general_configuration(general_data)
    trackconf = jc.config_json.get("tracks", [])
    for gnome in jc.genome_names:
        gtracks = jc.tracksToAdd[gnome]
        if len(gtracks) > 0:
            logging.debug(
                "for genome %s adding gtracks %s"
                % (gnome, json.dumps(gtracks, indent=2))
            )
            trackconf += gtracks
    jc.config_json["tracks"] = trackconf
    assconf = jc.config_json.get("assemblies", [])
    assconf += jc.assemblies
    jc.config_json["assemblies"] = assconf
    logging.debug(
        "assmeta=%s, first_contigs=%s, assemblies=%s, gnames=%s, trackidlist=%s, tracks=%s"
        % (
            jc.assmeta,
            jc.ass_first_contigs,
            json.dumps(assconf, indent=2),
            jc.genome_names,
            jc.trackIdlist,
            json.dumps(trackconf, indent=2),
        )
    )
    jc.write_config()
    jc.add_default_session(default_session_data)
    # note that this can be left in the config.json but has NO EFFECT if add_defsess_to_index is called.
    # jc.add_defsess_to_index(default_session_data)
    # jc.text_index() not sure what broke here.
