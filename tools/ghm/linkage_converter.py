#!/usr/bin/env python2
from bisect import bisect_left
from getopt import getopt
from sys import stderr, argv


class Converter:
    """
    Superclass for converting formats
    """
    haplo_individual_buffer = "%-8d %8d %8d %8d %2d %2d  "
    lod_line_buffer = "%8s %12s %8s %8s  %-s"

    def __init__(self, pedin, mapin):
        self.out_lod = "linkage.allegro_lod"
        self.out_haplo = "linkage.allegro_haplo"
        self.out_descent = "linkage.allegro_descent"

        self.haplo_map = {}  # fam_id -> indiv_id -> (all1,all2)
        self.descent_map = {}
        self.pos_marker = {}  # genpos -> marker
        self.marker_order = []
        self.lod_array = []
        self.pedigree = {}

        self.__populatePedigree(pedin)
        self.__populateMarkerMap(mapin)

    def __populatePedigree(self, input_ped):
        with open(input_ped, "r") as pio:
            for line in pio:
                tokens = map(int, line.split())
                f_id, p_id, father, mother, gender, affect = tokens[:6]

                if f_id not in self.pedigree:
                    self.pedigree[f_id] = {}
                if p_id not in self.pedigree:
                    self.pedigree[f_id][p_id] = (
                        father, mother, gender, affect
                    )
                else:
                    print >> stderr, "Duplicate individual:", f_id, p_id
                    exit(-1)

    def __populateMarkerMap(self, mapin):
        with open(mapin, "r") as mio:
            mio.readline()  # chomp header

            for line in mio:
                chrom, gpos, marker, ppos, nr = Converter.tokenizer(line)
                marker = marker.strip()
                self.pos_marker[float(gpos)] = marker
                self.marker_order.append(marker)

    def __annotateClosestMarker(self, pos_lod):

        # Produce sorted list of gpos
        pos_rsid_array = sorted(self.pos_marker.keys())
        pos_lod_array = sorted(pos_lod.keys())

        keys_to_annotate = {}

        # bisect and find closest marker
        for marker_pos in pos_rsid_array:
            marker = self.pos_marker[marker_pos]
            lm_ind = bisect_left(pos_lod_array, marker_pos)
            lm_pos = pos_lod_array[lm_ind - 1]

            diff1 = marker_pos - lm_pos
            diff1 = diff1 if diff1 > 0 else -diff1

            if diff1 > 0.1:
                continue

            if lm_pos not in keys_to_annotate:
                keys_to_annotate[lm_pos] = marker
            else:
                diff2 = marker_pos - lm_pos
                diff2 = diff2 if diff2 > 0 else -diff2

                if diff1 < diff2:
                    # update with closer marker
                    keys_to_annotate[lm_pos] = marker

        return keys_to_annotate, pos_lod_array

    def makeLODArray(self, pos_lod):
        (
            postns_with_closest_markers, sorted_poslod
        ) = self.__annotateClosestMarker(pos_lod)

        # update map with marker info
        for pos in sorted_poslod:
            marker = ""
            if pos in postns_with_closest_markers:
                marker = postns_with_closest_markers[pos]

            lod, alpha, hlod = pos_lod[pos]
            self.lod_array.append((pos, lod, alpha, hlod, marker))

    def __generateHeaders(self, npad_left=10):
        marker_order = self.marker_order
        max_len = -1
        markerpadd = []
        for marker in marker_order:
            nmark = len(marker)
            if nmark > max_len:
                max_len = nmark

        # paddleft
        for marker in marker_order:
            markerpadd.append(("%%-%ds" % max_len) % marker)

        # transpose
        buffer_left = ("%%%ds" % npad_left) % " "
        return '\n'.join(
            [
                buffer_left + "  ".join(x)
                for x in zip(*markerpadd)
            ][::-1]
        )

    @staticmethod
    def tokenizer(line):
        return line.splitlines()[0].split()

    # Override
    def extractLOD(self, file1):
        return ""

    # Override
    def extractHaplo(self, file1, extrafile=None):
        return ""

    # Override
    def extractDescent(self, file1):
        return ""

    def writeLOD(self):

        out_lod = open(self.out_lod, "w")

        header = Converter.lod_line_buffer % (
            "location", "LOD", "alpha", "HLOD", "marker"
        )
        print >> out_lod, header

        for pos, lod, alpha, hlod, marker in self.lod_array:
            out_line = Converter.lod_line_buffer % (
                "%.4f" % pos,
                "%.4f" % lod,
                "%.4f" % alpha,
                "%.4f" % hlod,
                marker
            )
            print >> out_lod, out_line
        out_lod.close()
        print >> stderr, "Wrote: ", out_lod.name

    def writeDescent(self):
        self.writeHaplo(True)

    def writeHaplo(self, descent=False):
        out_file = open(self.out_haplo
                        if not descent
                        else self.out_descent, "w")

        map_map = self.haplo_map if not descent else self.descent_map

        dummy_l = Converter.haplo_individual_buffer % (
            1, 1, 1, 1, 1, 1
        )  # *range(6) (unpack in py3 only...)
        headers = self.__generateHeaders(len(dummy_l)) + '\n'
        print >> out_file, headers

        for fam_id in map_map:
            for indiv_id in map_map[fam_id]:
                alleles = map_map[fam_id][indiv_id]
                ped_data = self.pedigree[fam_id][indiv_id]
                father, mother, gender, affect = ped_data

                indiv_data = Converter.haplo_individual_buffer % (
                    fam_id, indiv_id, father, mother, gender, affect
                )

                print >> out_file, indiv_data + " ".join(map(
                    lambda x: "%-2d" % x, alleles[0]
                ))
                print >> out_file, indiv_data + " ".join(map(
                    lambda x: "%-2d" % x, alleles[1]
                ))
        out_file.close()
        print >> stderr, "Wrote: ", out_file.name


class Swiftlink(Converter):

    def extractLOD(self, lodfile):
        pos_lod = {}
        with open(lodfile, 'r') as lio:
            lio.readline()

            for line in lio:
                if line[0] == "-":
                    pos, lod = map(float, Converter.tokenizer(line)[1:])
                    pos_lod[pos] = (lod, 1, 0)

        self.makeLODArray(pos_lod)


class Merlin(Converter):

    def extractDescent(self, file1):
        self.extractHaplo(file1, True)

    def extractHaplo(self, file1, use_flow=False):
        ped_map = {}

        tmp = [
            None,  # fam id
            [],    # [individuals]
            []     # [assosciated allele pairs]
        ]

        def flushTmpData(tmp):
            # finish populating alleles
            if len(tmp[1]) != len(tmp[2]):
                print >> stderr, "length mismatch"
                exit(-1)

            for tpa in range(len(tmp[1])):

                perc_alleles = tmp[2][tpa]
                perc_id = tmp[1][tpa]
                fam_id = tmp[0]

                if fam_id not in ped_map:
                    ped_map[fam_id] = {}

                if perc_id not in ped_map[fam_id]:
                    ped_map[fam_id][perc_id] = (
                        perc_alleles[0], perc_alleles[1]
                    )

            # clear
            tmp[1] = []
            tmp[2] = []

        with open(file1, 'r') as fio:

            for line in fio:

                if line.startswith("FAMILY"):
                    flushTmpData(tmp)

                    tmp[0] = int(line.split()[1])
                    continue

                if len(tmp[1]) == 0:  # hunt names after a flush
                    if line.find("(") != -1 and line.find(")") != -1:

                        people = [
                            x.strip() for x in line.splitlines()[0].split("  ")
                            if x.strip() != ""
                        ]

                        for p in range(len(people)):
                            perc = people[p].split(" ")
                            perc_id = int(perc[0])

                            # Add new perc to tmp array with blank alleles
                            tmp[1].append(perc_id)
                            tmp[2].append([[], []])
                    #
                    continue

                # Allele pairs
                trimmed = line.strip()
                if len(trimmed) == 0:
                    flushTmpData(tmp)
                    continue

                multiple_alleles = [
                    x.strip() for x in trimmed.split("   ")
                    if x.strip() != ""
                ]

                if len(multiple_alleles) != len(tmp[1]):
                    print >> stderr, "Num alleles and num percs mismatch"
                    exit(-1)

                for a in range(len(multiple_alleles)):
                    alleles = multiple_alleles[a]
                    left_b_right = alleles.split()

                    if not use_flow:
                        # pick first phasing
                        left_b_right[0] = int(
                            left_b_right[0]
                            .split(",")[0]
                            .replace("A", "").replace("?", "0")
                        )
                        left_b_right[2] = int(
                            left_b_right[2]
                            .split(",")[0]
                            .replace("A", "").replace("?", "0")
                        )
                    else:
                        # convert letters to numbers (A - Z)
                        left_b_right[0] = ord(left_b_right[0]) - 64
                        left_b_right[2] = ord(left_b_right[2]) - 64

                    tmp[2][a][0].append(left_b_right[0])
                    tmp[2][a][1].append(left_b_right[2])

            flushTmpData(tmp)

        if use_flow:
            self.descent_map = ped_map
        else:
            self.haplo_map = ped_map

    def extractLOD(self, file1):
        pos_lod = {}

        with open(file1, 'r') as fio:
            line = ""
            while not line.startswith(
                    "       POSITION        LOD      ALPHA       HLOD"
            ):
                line = fio.readline()

            tokens = Converter.tokenizer(fio.readline())

            while len(tokens) == 4:
                gpos = float(tokens[0])

                lod = tokens[1]
                if lod == "-INFINITY":
                    lod = -10000
                lod = float(lod)

                alph = float(tokens[2])
                hlod = float(tokens[3])

                pos_lod[gpos] = (lod, alph, hlod)
                tokens = Converter.tokenizer(fio.readline())

        self.makeLODArray(pos_lod)


class Simwalk(Converter):

    def extractHaplo(self, hef):
        self.extractHaploAndDescent(hef)

    def extractDescent(self, hef):
        self.extractHaploAndDescent(hef)

    def extractHaploAndDescent(self, hef):

        # clear marker map
        self.pos_marker = {}
        self.marker_order = []

        hio = open(hef, 'r')

        # extract markerloci
        line = ""
        while not line.startswith(" Name    Female   Male   Alleles"):
            line = hio.readline()
        line = hio.readline()

        while not line.startswith("Pedigree Name,       ID   Father"):
            line = hio.readline()
            if not line.startswith("  "):
                tokens = Converter.tokenizer(line)

                if len(tokens) != 4:
                    break

                gpos = float(tokens[2])
                marker = tokens[0].strip()
                self.pos_marker[gpos] = marker
                self.marker_order.append(marker)

        tmp = {
            "_fam": None,
            "_perc": None,
            "_allpat": [],  # alleles
            "_allmat": [],  #
            "_decpat": [],  # descent
            "_decmat": []
        }

        def insertDat():
            if len(tmp["_allpat"]) > 0:
                fid = tmp["_fam"]
                pid = tmp["_perc"]

                all1 = tmp["_allpat"]
                all2 = tmp["_allmat"]
                dec1 = tmp["_decpat"]
                dec2 = tmp["_decmat"]

                if fid not in self.haplo_map:
                    self.haplo_map[fid] = {}
                    self.descent_map[fid] = {}

                self.haplo_map[fid][pid] = (all1, all2)
                self.descent_map[fid][pid] = (dec1, dec2)

                tmp["_perc"] = None
                tmp["_allmat"] = []
                tmp["_allpat"] = []
                tmp["_decmat"] = []
                tmp["_decpat"] = []

        dashedlines_found = False

        for line in hio:
            if line.startswith("____"):
                if tmp["_perc"] is not None:
                    insertDat()

                dashedlines_found = True
                continue

            if dashedlines_found and not line.startswith(" "):
                fam = line.split("(")[0].strip()
                tmp["_fam"] = int(fam)
                dashedlines_found = False
                continue

            tokens = line.splitlines()[0].split()

            if tmp["_fam"] is not None and len(tokens) == 5:
                insertDat()
                tmp["_perc"] = int(tokens[0])

            # Allele data
            if (
                    tmp["_fam"] is not None and
                    tmp["_perc"] is not None and
                    len(tokens) == 6
            ):
                tokens = map(int, tokens)

                tmp["_allpat"].append(tokens[0])
                tmp["_allmat"].append(tokens[1])
                tmp["_decpat"].append(tokens[2])
                tmp["_decmat"].append(tokens[3])

        hio.close()
        insertDat()

    def extractLOD(self, scorefile):
        # This code is terrible, but simwalk is a horrible format

        pos_marker_lod = {}
        sio = open(scorefile, 'r')

        line = sio.readline()
        while not line.startswith("  NAME  , Haldane cM ,   alpha=1.00   ,"):
            line = sio.readline()

        sio.readline()  # skip empty

        recomb_fracts = False
        map_headers = False

        current_marker = None

        # Now add offsets and calc absolute cM
        pos_lod = {}

        for line in sio:
            if not recomb_fracts:
                if line[:2] == "__":
                    recomb_fracts = True
                    continue

                tokens = [x for x in line.split() if x != ","]
                if len(tokens) == 1 and not line.startswith(" "):
                    current_marker = tokens[0].splitlines()[0].strip()
                    if current_marker in pos_marker_lod:
                        print >> stderr, "Duplicate", current_marker
                        exit(-1)

                    pos_marker_lod[current_marker] = {}
                    continue

                if len(tokens) != 4:
                    continue

                gpos, lod, hlod, alpha = map(float, tokens)
                if gpos < 0 or gpos > 5:
                    continue

                pos_marker_lod[current_marker][gpos] = (lod, 1, hlod)

            else:
                if line.startswith(
                        "Haldane cM   NAME     FRACTION    OBSERVED & EXPECTED"
                ):
                    map_headers = True
                    continue

                if not map_headers:
                    continue

                if line[:2] == "__":
                    break

                tokens = Converter.tokenizer(line)
                if len(tokens) == 2:
                    gpos, marker = tokens
                    gpos = float(gpos)
                    self.pos_marker[gpos] = marker

                    offset_map = pos_marker_lod[marker]

                    for relative_gpos in offset_map:
                        abs_gpos = gpos + relative_gpos
                        if abs_gpos in pos_lod:
                            print >> stderr, "Duplicate abs cM", abs_gpos
                            exit(-1)

                        pos_lod[abs_gpos] = offset_map[relative_gpos]

        sio.close()
        self.makeLODArray(pos_lod)


class Genehunter(Converter):

    def extractHaplo(self, hapfile):
        with open(hapfile, 'r') as hio:
            line = "   "
            fam_id = -1

            try:
                while True:
                    line = hio.readline()
                    fam_id = int(line.split()[1])

                    if fam_id not in self.haplo_map:
                        self.haplo_map[fam_id] = {}

                    line = ""

                    while not line.startswith("*****"):
                        line = hio.readline()
                        indiv_all1 = map(int, line.split())
                        line = hio.readline()
                        indiv_all2 = map(int, line.split())

                        if indiv_all1 == [] and indiv_all2 == []:
                            return True

                        indiv_id = int(indiv_all1[0])

                        if len(indiv_all1) != len(indiv_all2) + 4:
                            print >> stderr, (
                                "Allele mismatch for indiv",
                                fam_id, indiv_id)
                            exit(-1)

                        self.haplo_map[fam_id][indiv_id] = (
                            indiv_all1[4:],
                            indiv_all2
                        )
            except IOError:
                hio.close()
                return True

    def extractLOD(self, file):

        pos_lod = {}  # pos -> lod

        with open(file, 'r') as fio:
            # Jump to marker positions
            line = ""
            while not line.startswith("Current map ("):
                line = fio.readline()
                # iterate

            # # Extract gts
            # while not line.startswith("npl"):
            #     line = fio.readline()
            #     gt_line += line.splitlines()[0] + " "

            # Get num pedigrees
            while not line.startswith("Totalling pedigrees"):
                line = fio.readline()
                # iterate

            # Jump to LOD scoring
            while not line.startswith("position LOD_score  (alpha, HLOD)"):
                line = fio.readline()

            # Extract lod
            line = fio.readline()

            while len(line) != 1:
                # print len(line), line
                (gpos, lod, alpha, hlod,
                 npl, pval, info) = Converter.tokenizer(line)
                gpos = float(gpos)

                if lod == "-INFINITY":
                    lod = -10000
                lod = float(lod)

                alpha = float(alpha.split("(")[-1].split(",")[0])
                hlod = float(hlod.split(")")[0])

                # insert
                if gpos not in pos_lod:
                    pos_lod[gpos] = [lod, alpha, hlod]
                else:
                    # update if new lod is larger
                    old_lod = pos_lod[gpos][0]
                    if lod > old_lod:
                        pos_lod[gpos][0] = lod

                line = fio.readline()
            fio.close()

        self.makeLODArray(pos_lod)


if __name__ == "__main__":

    def help():
        print >> stderr, '''%s <pedin> <mapin> <PROGRAM> [OPTIONS]

PROGRAM:  genehunter, merlin, simwalk, swiftlink
OPTIONS:
    -l lodfile
    -h haplofile
    -d descentfile
''' % argv[0]
        exit(-1)

    opts, rem = getopt(argv[4:], "-l:-h:-d:")

    if len(argv) < 3:
        help()

    pedin, mapin, program = argv[1:4]
    opts = dict(opts)

    obj = None
    if program == "genehunter":
        obj = Genehunter(pedin, mapin)
    elif program == "merlin":
        obj = Merlin(pedin, mapin)
    elif program == "simwalk":
        obj = Simwalk(pedin, mapin)
    elif program == "swiftlink":
        obj = Swiftlink(pedin, mapin)
    else:
        print >> stderr, "Nuts! No such prog:", program
        help()
        exit(-1)

    if '-l' in opts:
        obj.extractLOD(opts['-l'])
        obj.writeLOD()

    if '-h' in opts:
        obj.extractHaplo(opts['-h'])
        obj.writeHaplo()

    if '-d' in opts:
        obj.extractDescent(opts['-d'])
        obj.writeDescent()
