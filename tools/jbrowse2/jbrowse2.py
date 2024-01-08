#!/usr/bin/env python
# change to accumulating all configuration for config.json based on the default from the clone
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
    def __init__(self, outdir, genomes):
        self.debug = False
        self.usejson = True
        self.giURL = GALAXY_INFRASTRUCTURE_URL
        self.outdir = outdir
        os.makedirs(self.outdir, exist_ok=True)
        self.genome_paths = genomes
        self.genome_name = None
        self.genome_names = []
        self.trackIdlist = []
        self.tracksToAdd = []
        self.config_json = {}
        self.config_json_file = os.path.join(outdir, "config.json")
        self.clone_jbrowse()

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

        if trackDict.get("displays", None):
            style_data["type"] = trackDict["displays"]["type"]
            style_data["displayId"] = trackDict["displays"]["displayId"]
        return {"displays": [style_data]}

    def process_genomes(self):
        assemblies = []
        for i, genome_node in enumerate(self.genome_paths):
            if self.debug:
                log.info("genome_node=%s" % str(genome_node))
            genome_name = genome_node["meta"]["dataset_dname"].strip()
            if len(genome_name.split()) > 1:
                genome_name = genome_name.split()[0]
                # spaces and cruft break scripts when substituted
            fapath = genome_node["path"]
            assem = self.make_assembly(fapath, genome_name)
            assemblies.append(assem)
            self.genome_names.append(genome_name)
            if self.genome_name is None:
                self.genome_name = (
                    genome_name  # first one for all tracks - other than paf
                )
            if self.config_json.get("assemblies", None):
                self.config_json["assemblies"] += assemblies
            else:
                self.config_json["assemblies"] = assemblies

    def make_assembly(self, fapath, gname):
        faname = gname + ".fa.gz"
        fadest = os.path.join(self.outdir, faname)
        # fadest = os.path.realpath(os.path.join(self.outdir, faname))
        cmd = "bgzip -i -c %s -I %s.gzi > %s && samtools faidx %s" % (
            fapath,
            fadest,
            fadest,
            fadest,
        )
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
            "name": gname,
            "sequence": {
                "type": "ReferenceSequenceTrack",
                "trackId": gname,
                "adapter": adapter,
            },
            "rendering": {"type": "DivSequenceRenderer"},
        }
        return trackDict

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

    def text_index(self):
        # Index tracks
        args = [
            "jbrowse",
            "text-index",
            "--target",
            os.path.join(self.outdir, "data"),
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
        log.info("#### trackData=%s" % trackData)
        tId = trackData["label"]
        # can be served - if public.
        # dsId = trackData["metadata"]["dataset_id"]
        # url = "%s/api/datasets/%s/display?to_ext=hic " % (self.giURL, dsId)
        hname = trackData["name"]
        dest = os.path.join(self.outdir, hname)
        cmd = ["cp", data, dest]
        # these can be very big.
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
            "displays": [
                {
                    "type": "LinearHicDisplay",
                    "displayId": "%s-LinearHicDisplay" % tId,
                },
            ],
        }
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
            "displays": [
                {
                    "type": "LinearAlignmentsDisplay",
                    "displayId": "%s-LinearAlignmentsDisplay" % tId,
                },
            ],
        }
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
        cmd = ["tabix", "-f", "-p", "vcf", dest]
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
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

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
        # style_json = self._prepare_track_style(trackDict)
        # trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_paf(self, data, trackData, pafOpts, **kwargs):
        tname = trackData["name"]
        tId = trackData["label"]
        pgname = pafOpts["genome_label"]
        if len(pgname.split() > 1):
            pgname = pgname.split()[
                0
            ]  # trouble from spacey names in command lines avoidance
        asstrack, gname = self.make_assembly(pafOpts["genome"], pgname)
        self.genome_names.append(pgname)
        if self.config_json.get("assemblies", None):
            self.config_json["assemblies"].append(asstrack)
        else:
            self.config_json["assemblies"] = [
                asstrack,
            ]

        style_json = self._prepare_track_style(trackData)
        url = "%s.paf" % (trackData["label"])
        dest = "%s/%s" % (self.outdir, url)
        self.symlink_or_copy(os.path.realpath(data), dest)

        if self.usejson:
            trackDict = {
                "type": "SyntenyTrack",
                "trackId": tId,
                "assemblyNames": [self.genome_name, pgname],
                "name": tname,
                "adapter": {
                    "type": "PAFAdapter",
                    "pafLocation": {"uri": url},
                    "assemblyNames": [self.genome_name, pgname],
                },
                "config": style_json,
            }
            self.tracksToAdd.append(trackDict)
            self.trackIdlist.append(tId)
        else:
            self._add_track(
                trackData["label"],
                trackData["key"],
                trackData["category"],
                dest,
                assemblies=[self.genome_name, pgname],
                config=style_json,
            )

    def add_hicab(self, data, trackData, hicOpts, **kwargs):
        rel_dest = os.path.join("data", trackData["label"] + ".hic")
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
            "assemblyNames": [self.genome_name],
        }

        if query_refnames:
            json_track_data["adapter"]["refNamesQueryTemplate"]: query_refnames

        self.subprocess_check_call(
            [
                "jbrowse",
                "add-track-json",
                "--target",
                os.path.join(self.outdir, "data"),
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
                "style": {},
            }

            outputTrackConfig["key"] = track_human_label
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
            # Return non-human label for use in other fields
            yield outputTrackConfig["label"]

    def add_default_session(self, data):
        """
        Add some default session settings: set some assemblies/tracks on/off
        """
        tracks_data = []

        # TODO using the default session for now, but check out session specs in the future https://github.com/GMOD/jbrowse-components/issues/2708

        # We need to know the track type from the config.json generated just before
        track_types = {}
        logging.info("### add default session has data = %s\n" % str(data))
        with open(self.config_json_file, "r") as config_file:
            config_json = json.load(config_file)
        logging.info("### config.json read \n%s\n" % (config_json))

        for track_conf in self.tracksToAdd:  # config_json["tracks"]:
            track_types[track_conf["trackId"]] = track_conf["type"]
        logging.info(
            "### self.tracksToAdd = %s; track_types = %s"
            % (str(self.tracksToAdd), str(track_types))
        )

        for on_track in data["visibility"]["default_on"]:
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
        elif self.genome_name is not None:
            refName = self.genome_name
            start = 0
            end = 1000  # Booh, hard coded! waiting for https://github.com/GMOD/jbrowse-components/issues/2708

        if refName is not None:
            # TODO displayedRegions is not just zooming to the region, it hides the rest of the chromosome
            view_json["displayedRegions"] = [
                {
                    "refName": refName,
                    "start": start,
                    "end": end,
                    "reversed": False,
                    "assemblyName": self.genome_name,
                }
            ]

        session_name = data.get("session_name", "New session")

        # Merge with possibly existing defaultSession (if upgrading a jbrowse instance)
        session_json = {}
        if "defaultSession" in config_json:
            session_json = config_json["defaultSession"]

        session_json["name"] = session_name

        if "views" not in session_json:
            session_json["views"] = []

        session_json["views"].append(view_json)

        config_json["defaultSession"] = session_json

        with open(self.config_json_file, "w") as config_file:
            json.dump(config_json, config_file, indent=2)

    def add_general_configuration(self, data):
        """
        Add some general configuration to the config.json file
        """

        config_path = self.config_json_file
        config_json = {}
        if os.path.exists(config_path):
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
        if not config_json.get("configuration", None):
            config_json["configuration"] = {}
        config_json["configuration"].update(config_data)

        with open(config_path, "w") as config_file:
            json.dump(config_json, config_file, indent=2)

    def clone_jbrowse(self):
        """Clone a JBrowse directory into a destination directory."""
        dest = self.outdir
        cmd = ["jbrowse", "create", "-f", dest]
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
    parser.add_argument("--xml", help="Track Configuration")
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument("--version", "-V", action="version", version="%(prog)s 2.0.1")
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
    jc = JbrowseConnector(
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

    # .add_default_view() replace from https://github.com/abretaud/tools-iuc/blob/jbrowse2/tools/jbrowse2/jbrowse2.py
    default_session_data = {
        "visibility": {
            "default_on": [],
            "default_off": [],
        },
        "style": {},
        "style_labels": {},
    }

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
        track_conf["style"] = {
            item.tag: parse_style_conf(item) for item in track.find("options/style")
        }

        track_conf["style"] = {
            item.tag: parse_style_conf(item) for item in track.find("options/style")
        }

        track_conf["style_labels"] = {
            item.tag: parse_style_conf(item)
            for item in track.find("options/style_labels")
        }

        track_conf["conf"] = etree_to_dict(track.find("options"))
        keys = jc.process_annotations(track_conf)

        if keys:
            for key in keys:
                default_session_data["visibility"][
                    track.attrib.get("visibility", "default_off")
                ].append(key)
                default_session_data["style"][key] = track_conf[
                    "style"
                ]  # TODO do we need this anymore?
                default_session_data["style_labels"][key] = track_conf["style_labels"]

        default_session_data["defaultLocation"] = root.find(
            "metadata/general/defaultLocation"
        ).text
        default_session_data["session_name"] = root.find(
            "metadata/general/session_name"
        ).text

        general_data = {
            "analytics": root.find("metadata/general/analytics").text,
            "primary_color": root.find("metadata/general/primary_color").text,
            "secondary_color": root.find("metadata/general/secondary_color").text,
            "tertiary_color": root.find("metadata/general/tertiary_color").text,
            "quaternary_color": root.find("metadata/general/quaternary_color").text,
            "font_size": root.find("metadata/general/font_size").text,
        }
        track_conf["category"] = track.attrib["cat"]
        track_conf["format"] = track.attrib["format"]
        try:
            # Only pertains to gff3 + blastxml. TODO?
            track_conf["style"] = {t.tag: t.text for t in track.find("options/style")}
        except TypeError:
            track_conf["style"] = {}
            pass
        track_conf["conf"] = etree_to_dict(track.find("options"))
        jc.add_general_configuration(general_data)
        print("## processed", str(track_conf), "trackIdlist", jc.trackIdlist)
    x = open(args.xml, "r").read()
    log.info(
        "###done processing xml=%s; trackIdlist=%s, config=%s"
        % (x, jc.trackIdlist, str(jc.config_json))
    )
    jc.config_json["tracks"] = jc.tracksToAdd
    if jc.usejson:
        jc.write_config()
    # jc.add_default_view()
    jc.add_default_session(default_session_data)

    # jc.text_index() not sure what broke here.
