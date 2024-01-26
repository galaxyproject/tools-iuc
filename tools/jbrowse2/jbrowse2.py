 #!/usr/bin/env python

import argparse
import binascii
import datetime
import json
import logging
import os
import re
import shutil
import struct
import subprocess
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("jbrowse")

JB2VER = "v2.10.3"
# version pinned for cloning

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
                encoded_hist_id=metadata["history_id"],
                hist_name=metadata["history_display_name"],
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
    def __init__(self, outdir, jbrowse2path, genomes):
        self.giURL = GALAXY_INFRASTRUCTURE_URL
        self.outdir = outdir
        self.jbrowse2path = jbrowse2path
        os.makedirs(self.outdir, exist_ok=True)
        self.genome_paths = genomes
        self.genome_name = None
        self.genome_names = []
        self.trackIdlist = []
        self.tracksToAdd = []
        self.config_json = {}
        self.config_json_file = os.path.join(outdir, "config.json")
        self.clone_jbrowse()

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

    def subprocess_popen(self, command):
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
            raise RuntimeError("Command failed with exit code %s" % (retcode))

    def subprocess_check_output(self, command):
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
        wstyle = {
            "displays": [
                style_data,
            ]
        }
        return wstyle

    def process_genomes(self):
        assemblies = []
        useuri = False
        for i, genome_node in enumerate(self.genome_paths):
            if genome_node["useuri"].strip().lower() == "yes":
                useuri = True
            genome_name = genome_node["meta"]["dataset_dname"].strip()
            if len(genome_name.split()) > 1:
                genome_name = genome_name.split()[0]
                # spaces and cruft break scripts when substituted
            if genome_name not in self.genome_names:
                # pafs with shared references
                fapath = genome_node["path"]
                if not useuri:
                    fapath = os.path.realpath(fapath)
                assem = self.make_assembly(fapath, genome_name, useuri)
                assemblies.append(assem)
                self.genome_names.append(genome_name)
                if self.genome_name is None:
                    self.genome_name = (
                        genome_name  # first one for all tracks
                    )
                    self.genome_sequence_adapter = assem["sequence"]["adapter"]
                    self.genome_firstcontig = None
                    if not useuri:
                        fl = open(fapath, "r").readline()
                        fls = fl.strip().split(">")
                        if len(fls) > 1:
                            fl = fls[1]
                            if len(fl.split()) > 1:
                                self.genome_firstcontig = fl.split()[0].strip()
                            else:
                                self.genome_firstcontig = fl
                    else:
                        try:
                            fl = urllib.request.urlopen(fapath + ".fai").readline()
                        except:
                            fl = None
                        if fl:  # is first row of the text fai so the first contig name
                            self.genome_firstcontig = (
                                fl.decode("utf8").strip().split()[0]
                            )
                        else:
                            self.genome_firstcontig = None
        if self.config_json.get("assemblies", None):
            self.config_json["assemblies"] += assemblies
        else:
            self.config_json["assemblies"] = assemblies

    def make_assembly(self, fapath, gname, useuri):
        if useuri:
            faname = fapath
            adapter = {
                "type": "BgzipFastaAdapter",
                "fastaLocation": {
                    "uri": faname,
                    "locationType": "UriLocation"
                },
                "faiLocation": {
                    "uri": faname + ".fai",
                    "locationType": "UriLocation"
                },
                "gziLocation": {
                    "uri": faname + ".gzi",
                    "locationType": "UriLocation"
                }
            }
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
                }
            }

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
            ]
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
        tId = trackData["label"]
        # can be served - if public.
        # dsId = trackData["metadata"]["dataset_id"]
        # url = "%s/api/datasets/%s/display?to_ext=hic " % (self.giURL, dsId)
        hic_path = trackData.get("hic_path", None)
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            uri = data
        else:
            uri = "%s.hic" % trackData["label"]
            # slashes in names cause path trouble
            dest = os.path.join(self.outdir, uri)
            cmd = ["cp", data, dest]
            self.subprocess_check_call(cmd)
        categ = trackData["category"]
        trackDict = {
            "type": "HicTrack",
            "trackId": tId,
            "name":  trackData["name"],
            "assemblyNames": [self.genome_name],
            "category": [
                categ,
            ],
            "adapter": {
                "type": "HicAdapter",
                "hicLocation": { "uri": uri }
            }
        }
        self.tracksToAdd.append(trackDict)
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
        fname = "%s" % tId
        dest = "%s/%s" % (self.outdir, fname)
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
        mafs = open(data,'r').readlines()
        mafss = [x for x in mafs if (x.startswith('s\t') or x.startswith('s '))]
        samp = [x.split()[1] for x in mafss if len(x.split()) > 0]
        sampu = list(dict.fromkeys(samp))
        samples = [x.split('.')[0] for x in sampu]
        samples.sort()
        logging.warn("$$$$ cmd=%s, mafss=%s samp=%s samples=%s" % (' '.join(cmd), mafss, samp, samples))
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
            "assemblyNames": [self.genome_name],
            "displays": [
                {
                    "type": "LinearBasicDisplay",
                    "displayId": "%s-LinearBasicDisplay" % tId
                },
                {
                    "type": "LinearArcDisplay",
                    "displayId": "%s-LinearArcDisplay" % tId
                },
            ]
        }
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
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
        logging.warn("### blastxml to gff3 cmd = %s" % ' '.join(cmd))
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
            logging.warn("### gff3rebase cmd = %s" % ' '.join(cmd))
            gff3_rebased.close()
            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)
        url = "%s.gff3.gz" % trackData["label"]
        dest = "%s/%s" % (self.outdir, url)
        self._sort_gff(gff3, dest)
        tId = trackData["label"]
        categ = trackData["category"]
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)
        os.unlink(gff3)

    def add_bigwig(self, data, trackData):
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            url = "%s.bigwig" % trackData["label"]
            # slashes in names cause path trouble
            dest = os.path.join(self.outdir, url)
            cmd = ["cp", data, dest]
            self.subprocess_check_call(cmd)
        bwloc = {"uri": url}
        tId = trackData["label"]
        categ = trackData["category"]
        trackDict = {
            "type": "QuantitativeTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
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
        style_json = self._prepare_track_style(trackDict)
        trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_bam(self, data, trackData, bam_index=None, **kwargs):
        tId = trackData["label"]
        useuri = trackData["useuri"].lower() == "yes"
        bindex = bam_index
        categ = trackData["category"]
        if useuri:
            url = data
        else:
            fname = "%s.bam" % trackData["label"]
            dest = "%s/%s" % (self.outdir, fname)
            url = fname
            bindex = fname + ".bai"
            self.subprocess_check_call(["cp", data, dest])
            if bam_index is not None and os.path.exists(bam_index):
                if not os.path.exists(bindex):
                    # bai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
                    self.subprocess_check_call(["cp", bam_index, bindex])
                else:
                    # Can happen in exotic condition
                    # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
                    #      => no index generated by galaxy, but there might be one next to the symlink target
                    #      this trick allows to skip the bam sorting made by galaxy if already done outside
                    if os.path.exists(os.path.realpath(data) + ".bai"):
                        self.symlink_or_copy(os.path.realpath(data) + ".bai", bindex)
                    else:
                        log.warn("Could not find a bam index (.bai file) for %s", data)
        trackDict = {
            "type": "AlignmentsTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
            "assemblyNames": [self.genome_name],
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_cram(self, data, trackData, cram_index=None, **kwargs):
        tId = trackData["label"]
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            fname = "%s.cram" % trackData["label"]
            dest = "%s/%s" % (self.outdir, fname)
            url = fname
            self.subprocess_check_call(["cp", data, dest])
            if cram_index is not None and os.path.exists(cram_index):
                if not os.path.exists(dest + ".crai"):
                    # most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
                    self.subprocess_check_call(
                        ["cp", os.path.realpath(cram_index), dest + ".crai"]
                    )
            else:
                cpath = os.path.realpath(dest) + ".crai"
                cmd = ["samtools", "index", "-c", "-o", cpath, os.path.realpath(dest)]
                logging.debug("executing cmd %s" % " ".join(cmd))
                self.subprocess_check_call(cmd)
        trackDict = {
            "type": "AlignmentsTrack",
            "trackId": tId,
            "name": trackData["name"],
            "category": [
                categ,
            ],
            "assemblyNames": [self.genome_name],
            "adapter": {
                "type": "CramAdapter",
                "cramLocation": {"uri": url},
                "craiLocation": {
                    "uri": url + ".crai",
                },
                "sequenceAdapter": self.genome_sequence_adapter,
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_vcf(self, data, trackData):
        tId = trackData["label"]
        # url = "%s/api/datasets/%s/display" % (
        # self.giURL,
        # trackData["metadata"]["dataset_id"],
        # )
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "jbrowse sort-gff '%s' | bgzip -c > '%s'" % (
                data,
                dest,
            )  # "gff3sort.pl --precise '%s' | grep -v \"^$\" > '%s'"
            self.subprocess_popen(cmd)
            self.subprocess_check_call(["tabix", "-f", "-p", "gff", dest])

    def _sort_bed(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "sort -k1,1 -k2,2n '%s' | bgzip -c > '%s'" % (data, dest)
            self.subprocess_popen(cmd)
            cmd = ["tabix", "-f", "-p", "bed", dest]
            self.subprocess_check_call(cmd)

    def add_gff(self, data, ext, trackData):
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = trackData["path"]
        else:
            url = "%s.%s.gz" % (trackData["label"], ext)
            dest = "%s/%s" % (self.outdir, url)
            self._sort_gff(data, dest)
        tId = trackData["label"]
        categ = trackData["category"]
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_bed(self, data, ext, trackData):
        tId = trackData["label"]
        categ = trackData["category"]
        useuri = trackData["useuri"].lower() == "yes"
        if useuri:
            url = data
        else:
            url = "%s.%s.gz" % (trackData["label"], ext)
            dest = "%s/%s" % (self.outdir, url)
            self._sort_bed(data, dest)
        trackDict = {
            "type": "FeatureTrack",
            "trackId": tId,
            "name": trackData["name"],
            "assemblyNames": [self.genome_name],
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
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_paf(self, data, trackData, pafOpts, **kwargs):
        tname = trackData["name"]
        tId = trackData["label"]
        categ = trackData["category"]
        pgnames = [x.strip() for x in pafOpts["genome_label"].split(",") if len(x.strip()) > 0]
        pgpaths = [x.strip() for x in pafOpts["genome"].split(",") if len(x.strip()) > 0]
        passnames = [self.genome_name]  # always first
        logging.debug("### add_paf got pafOpts=%s, pgnames=%s, pgpaths=%s for %s" % (pafOpts, pgnames, pgpaths, tId))
        for i, gname in enumerate(pgnames):
            if len(gname.split()) > 1:
                gname = gname.split()[0]
            passnames.append(gname)
            # trouble from spacey names in command lines avoidance
            if gname not in self.genome_names:
                # ignore if already there - eg for duplicates among pafs.
                useuri = pgpaths[i].startswith("http://") or pgpaths[i].startswith(
                    "https://"
                )
                asstrack = self.make_assembly(pgpaths[i], gname, useuri)
                self.genome_names.append(gname)
                if self.config_json.get("assemblies", None):
                    self.config_json["assemblies"].append(asstrack)
                else:
                    self.config_json["assemblies"] = [
                        asstrack,
                    ]
        url = "%s.paf" % (trackData["label"])
        dest = "%s/%s" % (self.outdir, url)
        self.symlink_or_copy(os.path.realpath(data), dest)
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
            }
        }
        style_json = {
        "displays": [
            { "type": "LinearBasicDisplay",
            "displayId": "%s-LinearBasicyDisplay" % trackDict["trackId"]
            }
            ]
        }
        trackDict["style"] = style_json
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def process_annotations(self, track):
        category = track["category"].replace("__pd__date__pd__", TODAY)
        for i, (
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

            outputTrackConfig["key"] = track_human_label
            outputTrackConfig["useuri"] = useuri
            outputTrackConfig["path"] = dataset_path
            outputTrackConfig["ext"] = dataset_ext

            outputTrackConfig["trackset"] = track.get("trackset", {})
            outputTrackConfig["label"] = "%s_%i_%s" % (
                dataset_ext,
                i,
                track_human_label,
            )
            outputTrackConfig["metadata"] = extra_metadata
            outputTrackConfig["name"] = track_human_label

            if dataset_ext in ("gff", "gff3"):
                self.add_gff(
                    dataset_path,
                    dataset_ext,
                    outputTrackConfig,
                )
            elif dataset_ext in ("hic", "juicebox_hic"):
                self.add_hic(
                    dataset_path,
                    outputTrackConfig,
                )
            elif dataset_ext in ("cool", "mcool", "scool"):
                hic_url = "%s_%d.hic" % (track_human_label, i)
                hic_path = os.path.join(self.outdir, hic_url)
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
                    hic_url,
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
                    bam_index=real_indexes,
                )
            elif dataset_ext == "cram":
                real_indexes = track["conf"]["options"]["cram"]["cram_index"]
                self.add_cram(
                    dataset_path,
                    outputTrackConfig,
                    cram_index=real_indexes,
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
                logging.warn("Do not know how to handle %s", dataset_ext)
            # Return non-human label for use in other fields
            yield outputTrackConfig["label"]

    def add_default_session(self, default_data):
        """
        default session settings are hard and fragile.
        .add_default_view() and other configuration code adapted from
         https://github.com/abretaud/tools-iuc/blob/jbrowse2/tools/jbrowse2/jbrowse2.py
        """
        tracks_data = []
        # TODO using the default session for now, but check out session specs in the future https://github.com/GMOD/jbrowse-components/issues/2708
        track_types = {}
        with open(self.config_json_file, "r") as config_file:
            config_json = json.load(config_file)
        if self.config_json:
            config_json.update(self.config_json)
        for track_conf in self.tracksToAdd:
            tId = track_conf["trackId"]
            track_types[tId] = track_conf["type"]
            style_data = default_data["style"].get(tId,  None)
            if not style_data:
                logging.warn("### No style data in default data %s for %s" % (default_data, tId))
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
        # The view for the assembly we're adding
        view_json = {"type": "LinearGenomeView", "tracks": tracks_data}
        refName = None
        drdict = {
            "reversed": False,
            "assemblyName": self.genome_name,
            "start": 1,
            "end": 100000,
            "refName": "x",
        }

        if default_data.get("defaultLocation", ""):
            ddl = default_data["defaultLocation"]
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
        else:
            drdict["refName"] = self.genome_firstcontig
        if drdict.get("refName", None):
            # TODO displayedRegions is not just zooming to the region, it hides the rest of the chromosome
            view_json["displayedRegions"] = [
                drdict,
            ]

            logging.info("@@@ defaultlocation %s for default session" % drdict)
        else:
            logging.info(
                "@@@ no contig name found for default session - please add one!"
            )
        session_name = default_data.get("session_name", "New session")
        for key, value in mapped_chars.items():
            session_name = session_name.replace(value, key)
        # Merge with possibly existing defaultSession (if upgrading a jbrowse instance)
        session_json = {}
        if "defaultSession" in config_json:
            session_json = config_json["defaultSession"]

        session_json["name"] = session_name

        if "views" not in session_json:
            session_json["views"] = []

        session_json["views"].append(view_json)

        config_json["defaultSession"] = session_json
        self.config_json.update(config_json)

        with open(self.config_json_file, "w") as config_file:
            json.dump(self.config_json, config_file, indent=2)

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

    def clone_jbrowse(self, realclone=True):
        """Clone a JBrowse directory into a destination directory. This also works in Biocontainer testing now
        Leave as True between version updates on temporary tools - requires manual conda trigger :(
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
            cmd = ["rm", "-rf", os.path.join(dest, fn)]
            self.subprocess_check_call(cmd)
        cmd = ["cp", os.path.join(INSTALLED_TO, "jb2_webserver.py"), dest]
        self.subprocess_check_call(cmd)


def parse_style_conf(item):
    if item.text.lower() in ['false','true','yes','no']:
            return item.text.lower in ("yes", "true")
    else:
        return item.text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("--xml", help="Track Configuration")
    parser.add_argument(
        "--jbrowse2path", help="Path to JBrowse2 directory in biocontainer or Conda"
    )
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
        jbrowse2path=args.jbrowse2path,
        genomes=[
            {
                "path": x.attrib["path"],
                "label": x.attrib["label"],
                "useuri": x.attrib["useuri"],
                "meta": metadata_from_node(x.find("metadata")),
            }
            for x in root.findall("metadata/genomes/genome")
        ],
    )
    jc.process_genomes()

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
                track_conf["label"] = x.attrib["label"]
                trackkey = track_conf["label"]
                track_conf["useuri"] = x.attrib["useuri"]
                if is_multi_bigwig:
                    multi_bigwig_paths.append(
                        (
                            x.attrib["label"],
                            x.attrib["useuri"],
                            os.path.realpath(x.attrib["path"]),
                        )
                    )
                else:
                    if trackfiles:
                        metadata = metadata_from_node(x.find("metadata"))
                        track_conf["dataset_id"] = metadata["dataset_id"]
                        if x.attrib["useuri"].lower() == "yes":
                            tfa = (
                                x.attrib["path"],
                                x.attrib["ext"],
                                x.attrib["useuri"],
                                x.attrib["label"],
                                metadata,
                            )
                        else:
                            tfa = (
                                os.path.realpath(x.attrib["path"]),
                                x.attrib["ext"],
                                x.attrib["useuri"],
                                x.attrib["label"],
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
            )
        track_conf["category"] = track.attrib["cat"]
        track_conf["format"] = track.attrib["format"]
        track_conf["conf"] = etree_to_dict(track.find("options"))
        track_conf["category"] = track.attrib["cat"]
        track_conf["format"] = track.attrib["format"]
        keys = jc.process_annotations(track_conf)

        if keys:
            for key in keys:
                default_session_data["visibility"][
                    track.attrib.get("visibility", "default_off")
                ].append(key)
            if track.find("options/style"):
                default_session_data["style"][key] = {
                    item.tag: parse_style_conf(item) for item in track.find("options/style")
                }
            else:
                default_session_data["style"][key] = {}
                logging.warn("@@@@ no options/style found for %s" % (key))

            if track.find("options/style_labels"):
                default_session_data["style_labels"][key] = {
                    item.tag: parse_style_conf(item)
                    for item in track.find("options/style_labels")
                }
    default_session_data["defaultLocation"] = root.find(
        "metadata/general/defaultLocation"
    ).text
    default_session_data["session_name"] = root.find(
        "metadata/general/session_name"
    ).text
    logging.debug("default_session=%s" % (default_session_data))
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
    trackconf = jc.config_json.get("tracks", None)
    if trackconf:
        jc.config_json["tracks"].update(jc.tracksToAdd)
    else:
        jc.config_json["tracks"] = jc.tracksToAdd
    jc.write_config()
    jc.add_default_session(default_session_data)
    # jc.text_index() not sure what broke here.
