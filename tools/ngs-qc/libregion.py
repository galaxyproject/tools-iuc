#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 26/03/15

@author: Matthias Blum
"""

import libngs
import os
import sqlite3
import sys
import tempfile

from subprocess import Popen, PIPE


def median(lst):
    """
    Return a list's median
    :param lst:     List (sorted or not) of int
    :return:        int
    """
    lst.sort()
    length = len(lst)

    if length == 0:
        return None
    elif length % 2 == 0:
        # Even number of elements
        return (lst[int((length+1)/2 - 1)] + lst[int((length+1)/2)]) / 2.0
    else:
        return lst[int((length+1)/2 - 1)]


def make_region_around_gene(genome, gene, resolution, centerontss=False):
    """
    Look for a gene in a genome and return a genomic region of <resolution> bp aroung the gene if found
    :param genome:          Genome assembly
    :param gene:            Gene name
    :param resolution:      Resolution
    :param centerontss:     If True, center the region on the gene's TSS,
    otherwise center the region on the gene's center
    :return:
    """
    region = {}

    if genome.startswith('TAIR'):
        genome = 'TAIR10'

    con = sqlite3.connect(os.path.join(os.path.dirname(__file__), 'db', 'genes.db'))
    cur = con.cursor()
    cur.execute('SELECT chrm, start, end, strand, name '
                'FROM gene '
                'WHERE assembly=? COLLATE NOCASE AND name=? COLLATE NOCASE', (genome, gene))
    row = cur.fetchone()

    if row:
        if centerontss:
            center = row[1] if row[3] == '+' else row[2]
        else:
            center = (row[1] + row[2]) / 2

        r_start = center - resolution / 2
        r_end = center + resolution / 2
        gene = row[4]  # Gene's corrected name

        # Sorted list, but we want <gene> to be the first element of the list (more important)
        genes = find_genes(genome, row[0], r_start, r_end)
        idx = [i[3] for i in genes].index(gene)
        # Swap <gene> with the first gene of the list
        genes[0], genes[idx] = genes[idx], genes[0]

        region = {
            'chrm': row[0],
            'start': center - resolution / 2,
            'end': center + resolution / 2,
            'genes': genes
        }

    cur.close()
    con.close()

    return region


def find_genes(genome, chrm, start_pos, end_pos, minsize=1):
    """
    Look for genes present in a given genomic region
    :param genome:      Genome assembly
    :param chrm:        Chromosome
    :param start_pos:   Minimal position
    :param end_pos:     Maximal position
    :param minsize:     Minimum size
    :return:            List of gene (a gene: [start, end, strand, name])
    """
    if genome.startswith('TAIR'):
        # Same DB for TAIR8, TAIR9, and TAIR10
        genome = 'TAIR'

    genes = {}
    con = sqlite3.connect(os.path.join(os.path.dirname(__file__), 'db', 'genes.db'))
    cur = con.cursor()

    if not chrm.lower().startswith('chr'):
        chrm = 'chr{}'.format(chrm)

    query = ('SELECT start, end, strand, name '
             'FROM gene '
             'WHERE assembly=:assembly COLLATE NOCASE AND chrm=:chrm COLLATE NOCASE '
             'AND ((start BETWEEN :min AND :max) OR (end BETWEEN :min AND :max) OR (:min BETWEEN start AND end))')
    cur.execute(query, {'assembly': genome, 'chrm': chrm, 'min': start_pos, 'max': end_pos})
    for row in cur.fetchall():
        g_start = row[0]
        g_end = row[1]

        if g_end - g_start >= minsize:
            # Gene is long enough to be kept

            # Useless to have genes that go outside of the region
            if g_start < start_pos:
                g_start = start_pos
            if g_end > end_pos:
                g_end = end_pos

            name = row[3]
            if name in genes:
                # Transcript: several genes with the same name but different positions
                # Merging them if overlapping
                if g_start < genes[name][0]:
                    genes[name][0] = g_start

                if g_end > genes[name][1]:
                    genes[name][1] = g_end
            else:
                genes[name] = [g_start, g_end, row[2]]

    cur.close()
    con.close()

    # Return genes sorted by start position
    return sorted([v + [k] for k, v in genes.items()], key=lambda x: x[0])


def filter_by_gene(regions, genome, resolution=0, minsize=0, mingenes=0):
    """
    Filter regions with at least <mingenes> genes in the centered sub-region of <resolution> bp
    Let's say we have a region of 500kb (0-500kb on the genome) and a given <resolution> of 50kb.
    The region will be selected if we found genes in 225kb-275kb (center of the region +/- 50kb/2)
    :param regions:     Dictionary of regions (key: chrm, value: list of regions [start, end, reads, bin_start, sd])
    :param genome:      Genome assembly
    :param resolution:  Sub-region resolution
    :param mingenes:    Minimum number of genes to keep the region
    :return:            List of genes (gene: start, end, strand, name)
    """
    filtered_regions = {}
    for chrm in regions:
        filtered_regions[chrm] = []

    for chrm, reg_lst in regions.items():
        for r in reg_lst:
            r_start = r[0]
            r_end = r[1]
            r_length = r_end - r_start

            # If <resolution> is longer than the region,
            # it means we're going to find genes that are outside the region, which has no sense.
            if resolution and resolution < r_length:

                region_center = r_start + r_length * 0.5  # Start + region length / 2
                min_pos = region_center - resolution * 0.5
                max_pos = region_center + resolution * 0.5
            else:
                min_pos = r_start
                max_pos = r_end

            genes = find_genes(genome, chrm, start_pos=min_pos, end_pos=max_pos, minsize=minsize)

            if len(genes) >= mingenes:
                r.append(genes)
                filtered_regions[chrm].append(r)

    return filtered_regions


def get_disp(localqc_table, regions, bg=1):
    """
    Find localQC bins for each region
    :param localqc_table:   Table of localQC regions
    :param regions:         List dictionary of regions
    :param bg:              Background threshold
    :return:
    """
    # Group regions by chromosome
    regions_chrm = {}
    for i, r in enumerate(regions):
        r['id'] = i
        r['disps'] = []
        if r['chrm'] not in regions_chrm:
            regions_chrm[r['chrm']] = [r]
        else:
            regions_chrm[r['chrm']].append(r)

    with open(localqc_table) as fh:
        for line in fh:
            cols = line.rstrip().split('\t')

            if cols[0] in regions_chrm and int(cols[3]) >= bg:
                bin_start = int(cols[1])

                for r in regions_chrm[cols[0]]:

                    if r['start'] <= bin_start <= r['end']:
                        r['disps'].append([bin_start, abs(float(cols[9]))])

    regions = [None] * len(regions)
    for chrm, reg_lst in regions_chrm.items():
        for r in reg_lst:
            idx = r.pop('id')
            regions[idx] = r

    return regions


def enqueue_regions(regions, limit):
    """
    Keep the best regions of each chromosome
    :param regions:     Dictionary of regions (key: chrm, value: [start, end, reads, bin_start, sd, genes])
    :param limit:       Number of regions to keep
    :return:
    """
    queue = []
    queue2 = []
    empty_chrms = {}  # chromosomes without any region

    while len(queue) < limit:
        for chrm, reg_lst in regions.items():
            if reg_lst:
                # Remove the first (best) region of the list and add it to the queue
                r = reg_lst.pop(0)
                queue.append({
                    'chrm': chrm,
                    'start': r[0],
                    'end': r[1],
                    'binstart': r[3],
                    'genes': r[5]
                })
            elif chrm not in empty_chrms:
                empty_chrms[chrm] = 1

        if len(empty_chrms) == len(regions):
            break

    # Place all remaining regions in the second queue
    for chrm, reg_lst in regions.items():
        # Add the chrm to each region, now it's: [chrm, start, end, reads, bin_start, sd, genes]
        queue2 += [[chrm] + r for r in reg_lst]

    # Sort the supplementary regions by standard deviation
    queue2.sort(key=lambda x: x[5])

    for r in queue2:
        queue.append({
            'chrm': r[0],
            'start': r[1],
            'end': r[2],
            'binstart': r[4],
            'genes': r[6]
        })

    return queue


def filter_by_overlap(regions):
    """
    Iterate all regions and remove those that overlap more than 25% with a better region
    :param regions: Dictionary of regions (key: chrm, value: list of regions [start, end, reads, bin_start sd])
    :return:
    """
    # For each chromosome, we look if some regions are overlapping more than 25%
    for chrm, reg_lst in regions.items():
        regions2 = []
        if len(reg_lst) > 1:
            # List of regions to skip because they are overlapping a better region
            to_skip = []

            for i, r1 in enumerate(reg_lst):
                if i in to_skip:
                    continue

                # Add the region
                regions2.append(r1)
                resolution = r1[1] - r1[0]
                start_limit = r1[0] + 0.25 * resolution
                end_limit = r1[1] - 0.25 * resolution

                # Compare it with region not as good (since sorted)
                for j, r2 in enumerate(reg_lst):
                    if j > i:
                        if start_limit < r2[0] < end_limit or start_limit < r2[1] < end_limit:
                            # r2 is overlapping with r1, so we're going to discard it
                            to_skip.append(j)

            regions[chrm] = regions2

    return regions


def sortregions(regions, total_reads):
    """
    Sort regions by standard deviation, to avoid having too big peaks
    :param regions:         Dictionary of regions (key: chrm, value: list of regions [start, end, reads, bin_start])
    :param total_reads:     Sum of all reads in regions
    :return:                Dictionary of regions (key is the chromosome)
    """
    if len(total_reads) == 0:
        return None

    med = median(total_reads)
    mean = float(sum(total_reads)) / len(total_reads)

    # We expect the mean is greater than the median
    v = med + abs(mean - med)

    for chrm, reg_lst in regions.items():
        for r in reg_lst:
            # Standard deviation: absolute difference between the number of reads in bin and the mean of reads
            # Now: [start, end, reads, bin_start, sd]
            r.append(abs(r[2] - v))

        # Sort by the standard deviation
        reg_lst.sort(key=lambda x: x[4])

    return regions


def select_regions(localqc_table, genomeinfo, resolution, elongation=150, maxdisp=2.5):
    """
    Select regions to display based on their dispersion and their read count intensity
    :param localqc_table:   LocalQC regions file
    :param genomeinfo:      Dictionary of chromosomes' size
    :param resolution:      Total size of the desired regions (bp)
    :param elongation:      Length (bp) to stretch the region in both directions
    :param maxdisp:         Dispersion filtering (the lower the better)
    :return:                A list of regions and the total number of reads in all selected regions
    """
    regions = {}
    total_reads = []

    with open(localqc_table) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            chrm = cols[0]

            # Dispersion 10% for 50% sampling
            # Use the new localqc table format (10 columns instead of 13)!
            disp = abs(float(cols[9]))

            if disp <= maxdisp and not chrm.lower().startswith(('chrx', 'chry', 'chrm')):
                total_reads.append(int(cols[3]))
                bin_start = int(cols[1])
                start = bin_start - resolution / 2 - elongation
                end = int(cols[2]) + resolution / 2 + elongation

                # Shift region start/end if it's outside of the chromosome
                if 0 <= start < end <= genomeinfo[chrm]:
                    # Region entirely in the chromosome

                    if chrm not in regions:
                        regions[chrm] = []

                    # Store chrm, start, end, reads, bin_start
                    regions[chrm].append([start, end, int(cols[3]), bin_start])

    return regions, total_reads


def get_wigs(wigit, regions, bed):
    """
    Parse the BED file and attribute wiggles to the given regions
    :param regions:     List of regions
    :param bed:    BED file
    :return:
    """
    # Group regions by chromosome
    regions_chrm = {}
    for i, r in enumerate(regions):
        r['id'] = i
        r['wigs'] = []
        if r['chrm'] not in regions_chrm:
            regions_chrm[r['chrm']] = [r]
        else:
            regions_chrm[r['chrm']].append(r)

    params = ['{}:{}:{}'.format(r['chrm'], r['start'], r['end']) for r in regions]
    cmd = [wigit, bed, '--sorted'] + params

    pop = Popen(cmd, stdout=PIPE)
    parse= False
    for i, line in enumerate(pop.stdout):
        if line.decode('utf-8').startswith('variableStep'):
            # Line like: variableStep chrom=chr21 span=20
            chrm = line.decode('utf-8').split()[1][6:]

            if chrm in regions_chrm:
                parse = True
            else:
                parse = False
        elif parse:
            try:
                wig, reads = line.decode('utf-8').rstrip().split('\t')
                wig = int(wig)
                reads = int(reads)
            except:
                continue
            else:
                for r in regions_chrm[chrm]:
                    if r['start'] <= wig <= r['end']:
                        r['wigs'].append([wig, reads])
    pop.wait()

    regions = [None] * len(regions)

    for chrm, lst_reg in regions_chrm.items():
        for r in lst_reg:
            idx = r.pop('id')
            regions[idx] = r

    return regions


def colorbar(dest, gnuplot='gnuplot'):
    """
    Draw a colorbar
    :param dest:        Output image
    :param gnuplot:     Gnuplot binary path
    :return:
    """
    fh, path = tempfile.mkstemp()
    os.close(fh)

    with open(path, 'w') as fo:
        # Define terminal and output
        fo.write('set terminal pngcairo nocrop enhanced font "sans,10" size 100,400\n')
        fo.write('set output "{0}"\n'.format(dest))

        # Define the colorbox's params
        fo.write('set colorbox vertical user origin 0.1, 0.02 size 0.2,0.96 front\n')

        # Define the palette
        fo.write('set palette defined (0 "#000000", 1 "#331900", 2 "#663300", 3 "#994c00", 4 "#cc6600", 5 "#ff7f00", 6 "#ff9800", 7 "#ffb100", 8 "#ffcb00", 9 "#fff400")\n')

        # Colorbox's title
        fo.write('set cblabel "Enrichment quality" font "sans,15"\n')

        # Colorbox's tics
        fo.write('set cbtics ("Low" 0,"High" 10) font "sans,15"\n')

        # Define the plot (we don't care since we won't plot anything but it's required)
        fo.write('set xrange [0:1]\n')
        fo.write('set yrange [0:1]\n')

        # Define colorbox range (dispersion so 0-10)
        fo.write('set cbrange [0:10]\n')

        # Hide everything related to the plot (tics, title, borders)
        fo.write('unset xtics\n')
        fo.write('unset ytics\n')
        fo.write('unset key\n')
        fo.write('unset border\n')

        # Plot nothing :)
        fo.write('plot -1 with boxes palette\n')

    with open(os.devnull, 'w') as dn:
        pop = Popen([gnuplot, path], stdout=dn, stderr=dn)
        pop.wait()
    os.unlink(path)
    return pop.returncode


def place_labels(genes, region_start, region_end, num_letters=196):
    """
    Handle the collision between genes' label by placing overlapping labels on different levels
    :param genes:           List of dictionary of genes (gene: [start, end, strand, name])
    :param region_start:    Region's start position
    :param region_end:      Region's end position
    :param num_letters:     226 for font size 8, 195 for font size 9 (text used: 'HELLO WORLD ')
    :return:                List of lists (of genes)
    """
    # Size (bp) of a letter: region's size / number of letters that can be displayed in the region
    letter_bp = (region_end - region_start) / num_letters

    levels = [[]]

    for gene in genes:
        # Elements (5): name, start, end, strand, emphasis

        on_this_level = True

        strand_sign = '> ' if gene[2] == '+' else '< '
        label_start = gene[0]
        # Label's end is equal to the label's start plus the length of the full label (gene's name + strand sign (2))
        label_end = label_start + (len(gene[3]) + 2) * letter_bp

        if label_end > region_end:
            # Align right so we shift label's positions
            label_start -= label_end - region_end + letter_bp
            label_end = region_end - letter_bp
            align = 'right'
        else:
            align = 'left'

        for i, lvl in enumerate(levels):
            on_this_level = True

            if not lvl:
                # Empty level: allocating here
                levels[i].append([label_start, label_end, gene[3], strand_sign, align])
            else:
                for g in lvl:
                    if g[0] <= label_start <= g[1] or label_start <= g[0] <= label_end:
                        # Already a gene allocated here: we'll see on the next level
                        on_this_level = False
                        break

                if on_this_level:
                    # There is enough space for this label!
                    levels[i].append([label_start, label_end, gene[3], strand_sign, align])
                    break

        if not on_this_level:
            # Couldn't find a nice small place for this label: need to add another level
            levels.append([[label_start, label_end, gene[3], strand_sign, align]])

    return levels


def place_genes(genes, limit=0):
    """
    Handle the collision between genes by placing overlapping genes on different levels
    :param genes:   Sorted list of genes (gene: [start, end, strand, name])
    :return:        List of lists (of genes)
    """
    levels = [[]]

    # For each gene within the same region
    for gene in genes:
        on_this_level = True

        # For each level
        for i, lvl in enumerate(levels):
            on_this_level = True

            if not lvl:
                # Level empty: allocate here
                levels[i].append(gene)
            else:
                # For each gene already allocated on the current level
                for g in lvl:
                    if g[0] <= gene[0] <= g[1] + 500:
                        # Already a gene allocated where the current gene starts: need to go deeper (next level)
                        on_this_level = False
                        break

                if on_this_level:
                    # Allocate on this level
                    levels[i].append(gene)
                    break

        if not on_this_level:
            if limit and len(levels) == limit:
                return levels

            # Add a new level and allocate the current gene within in
            levels.append([gene])

    return levels


def shrink_region(region, resolution):
    """
    Reduce a region's size and get rid of the dispersions/wigs/genes outside the new region
    :param region:              Dictionary of the region to shrink
    :param new_resolution:      Size in bp of the new region
    :return:
    """
    midpoint = (region['start'] + region['end']) / 2
    new_start = midpoint - resolution / 2
    new_end = midpoint + resolution / 2

    new_disps = [disp for disp in region['disps'] if new_start <= disp[0] < new_end]
    new_wigs = [wig for wig in region['wigs'] if new_start <= wig[0] < new_end]
    new_genes = []

    for gene in region['genes']:
        if new_start <= gene[0] < new_end or new_start < gene[1] <= new_end or gene[0] <= new_start < new_end <= gene[1]:
            if gene[0] < new_start:
                gene[0] = new_start
            if gene[1] > new_end:
                gene[1] = new_end
            new_genes.append(gene)

    return {
        'chrm': region['chrm'],
        'start': new_start,
        'end': new_end,
        'disps': new_disps,
        'genes': new_genes,
        'wigs': new_wigs,
        'binstart': region.get('binstart')
    }


def plot(region, dest, **kwargs):
    """
    Plot the wiggles, the dispersion, and the gene of a given region
    :param region:  Dictionary of the given region (key are: chrm, start, end, genes, wigs, disps)
    :param dest:    Output image
    :param kwargs:
    :return:
    """
    details = kwargs.get('details', True)           # If True, draw each gene instead of grouping them
    gnuplot = kwargs.get('gnuplot', 'gnuplot')      # Gnuplot binary
    bg = kwargs.get('bg', 0)                        # Background threshold. If given, draw a horizontal line.
    binsize = kwargs.get('binsize', 500)            # Size of each bin

    region_length = region['end'] - region['start']
    xtics = []

    for i in range(5):
        x = region['start'] + region_length / 6. * (i + 1)
        xtics.append(int(x) / 100 * 100)

    xtics = ','.join('"{}" {}'.format(libngs.num_format(x), x) for x in xtics)

    fd, wigs_file = tempfile.mkstemp()
    os.close(fd)

    if region['wigs']:
        if region.get('binstart') is not None:
            idx = region.get('binstart') / 500
            try:
                # Find the highest peak in the bin we selected
                highest_peak_bin = max([wig[1] for wig in region['wigs'] if wig[0]/500 == idx])
            except ValueError as e:
                # max() arg is an empty sequence: there is no wig in the selected bin (weird)
                try:
                    highest_peak_bin = max([wig[1] for wig in region['wigs']])
                except ValueError as e:
                    highest_peak_bin = 0

            try:
                highest_peak_region = max([wig[1] for wig in region['wigs']])
            except ValueError as e:
                # wtf? That means there is no bloody wig in the entire region?
                highest_peak_region = 0

            if highest_peak_region > highest_peak_bin * 5:
                # Intensity of the region is 75% of the image height
                ymax = int(highest_peak_bin / 0.75)
            else:
                ymax = (highest_peak_region / 10 + 1) * 10
        else:
            ymax = (max([wig[1] for wig in region['wigs']]) / 10 + 1) * 10
    else:
        ymax = 10

    # Write wigs file
    with open(wigs_file, 'w') as fo:
        for wig in region['wigs']:
            fo.write('{}\t{}\n'.format(wig[0], wig[1]))

    fd, disps_file = tempfile.mkstemp()
    os.close(fd)

    # Write disps file
    with open(disps_file, 'w') as fo:
        for disp in region['disps']:
            # Add <binsize>/2 to the local QC position in order to align the left of the box with the position
            # Otherwise, the box is centered on the position
            fo.write('{}\t1\t{}\n'.format(disp[0]+binsize/2, disp[1]))

    region['genes'] = region['genes'][-4:]
    gene_levels = place_genes(region['genes'], limit=3)

    # Keep only selected genes and sort them by start position
    region['genes'] = sorted([g for lvl in gene_levels for g in lvl], key=lambda x: x[0])
    if details:
        label_levels = place_labels(region['genes'], region['start'], region['end'])
    else:
        label_levels = []

    if details:
        cnt_levels = len(gene_levels) + len(label_levels)
        subplots_ratio = [0] * 3

        # For 8 levels, the genes plot needs a ratio of 0.4
        subplots_ratio[2] = cnt_levels * 0.4 / 8
        # Minimal height
        if subplots_ratio[2] < 0.25:
            subplots_ratio[2] = 0.25
        # The dispersion plot has a static ratio
        subplots_ratio[1] = 0.1
        # We adapt the last plot (wigs)
        subplots_ratio[0] = 1 - subplots_ratio[1] - subplots_ratio[2]
    else:
        # Had to play with the ratio, since the last plot (gene) has some of its height used for the bottom margin
        subplots_ratio = [0.8, 0.05, 0.15]

    # We need to adapt the levels' height/margin
    # For a ratio of 0.4 (8 levels), height: 0.05 and margin: 0.01
    # If the ratio decrease (less levels), we increase the height/margin
    level_height = 0.4 / subplots_ratio[2] * 0.05
    level_margin = 0.4 / subplots_ratio[2] * 0.01

    fd, script = tempfile.mkstemp()
    os.close(fd)

    with open(script, 'w') as fo:
        # Output to PNG, with Sans font
        fo.write('set terminal pngcairo nocrop font "DejaVuSans,9" size 1600,400\n')
        fo.write('set output "{}"\n'.format(dest))

        fo.write('set tmargin 2\n')
        fo.write('set bmargin 0\n')  # No space between plots

        # Plots correctly aligned
        fo.write('set lmargin 7\n')
        fo.write('set rmargin 7\n')

        fo.write('set multiplot\n')

        ###########################
        # Plot A
        ###########################
        fo.write('set origin 0, {}\n'.format(1 - subplots_ratio[0]))
        fo.write('set size 1, {}\n'.format(subplots_ratio[0]))

        # Title (format: <chrm>:<start>-<end>)
        if not region['chrm'].lower().startswith('chr'):
            region['chrm'] = 'chr{}'.format(region['chrm'])

        fo.write('set title "{}: {} - {}" font "DejaVuSans,14"\n'.format(region['chrm'],
                                                                         libngs.num_format(region['start']),
                                                                         libngs.num_format(region['end'])))

        # Hide legend
        fo.write('set nokey\n')

        # If we ever want a grid. But:
        # 1) Hinrich does not like them
        # 2) We have not ytics so there wouldn't have any horizontal lines
        # 3) We would unset xtics to have remove vertical lines
        # fo.write('set grid\n')

        # Filled bars
        fo.write('set style fill solid 1.00 \n')

        # Delete auto-labelling for xtics, hide tic marks (scale) and hide mirror tics
        fo.write('set xtics format " " scale 0 nomirror\n')

        # Define the Y-axis range and display only the top tic
        fo.write('set yrange [0:{}]\n'.format(ymax))

        if bg:
            fo.write('set ytics ("{0}" {0}, "{1}" {1})\n'.format(bg, ymax))
        else:
            fo.write('set ytics ("{0}" {0})\n'.format(ymax))

        # Region start/end
        fo.write('set xrange[{}:{}]\n'.format(region['start'], region['end']))

        # Wig-sized bars
        fo.write('set boxwidth 20\n')

        if bg:
            fo.write('set object 1 rectangle from {},0 '
                     'to {},{} fc rgb "#BFBFBF" fs solid noborder\n'.format(region['start'], region['end'], bg))

        if region['wigs']:
            fo.write('plot "{}" using 1:2 with boxes lt rgb "#000000"\n'.format(wigs_file))
        else:
            fo.write('plot 0\n')

        if bg:
            fo.write('unset object 1\n')

        ###########################
        # Plot B
        ###########################
        fo.write('set origin 0, {}\n'.format(1 - subplots_ratio[0] - subplots_ratio[1]))
        fo.write('set size 1, {}\n'.format(subplots_ratio[1]))

        # No space between plots
        fo.write('set tmargin 0\n')

        # Title only for top plot
        fo.write('unset title\n')

        # Dispersion is always between 0 and 1
        fo.write('set yrange [0:1]\n')

        # Turn off ytics
        fo.write('unset ytics\n')

        # Colorbox range (dispersion from 0 to 10)
        fo.write('set cbrange [0:10]\n')

        # Hide the colorbox
        fo.write('unset colorbox\n')

        # Bin-sized bars
        fo.write('set boxwidth {}\n'.format(binsize))

        # Define colours for colorbox
        fo.write('set palette negative defined (0 "#000000", 1 "#331900", 2 "#663300", 3 "#994c00", '
                 '4 "#cc6600", 5 "#ff7f00", 6 "#ff9800", 7 "#ffb100", 8 "#ffcb00", 9 "#fff400")\n')

        if region['disps']:
            fo.write('plot "{}" using 1:2:3 with boxes palette\n'.format(disps_file))
        else:
            fo.write('plot 0\n')

        ###########################
        # Plot C
        ###########################
        fo.write('set origin 0, 0\n')
        fo.write('set size 1, {}\n'.format(subplots_ratio[2]))

        # Margin for xtics
        fo.write('set bmargin 2\n')
        fo.write('set xtics ({}) font "DejaVuSans,10" scale 0.75 out\n'.format(xtics))

        # Default style for rectangle
        # For transparency, add "fs transparent solid 0.5" after the fillcolor definition
        if details:
            fo.write('set style rectangle fc rgb "#0088cc" fs solid noborder\n')
        else:
            fo.write('set style rectangle fc rgb "#0088cc" fs transparent solid 0.5 noborder\n')

        y = 1  # Top

        # Draw genes
        if details:
            for lvl in gene_levels:
                y -= level_height + level_margin
                for gene in lvl:
                    fo.write('set object rect from {},{} to {},{}\n'.format(gene[0], y, gene[1], y+level_height))
        else:
            for lvl in gene_levels:
                for gene in lvl:
                    fo.write('set object rect from {},0.1 to {},0.9\n'.format(gene[0], gene[1]))

        # Draw labels
        for lvl in label_levels:
            y -= level_height * 2 + level_margin

            for label in lvl:
                if label[4] == 'left':
                    label_cmd = "set label '{}' at {},{} left front".format(label[3] + label[2], label[0], y)
                else:
                    label_cmd = "set label '{}' at {},{} right front".format(label[3] + label[2], label[1], y)

                """if label[4]:
                    # emphasis on gene
                    label_cmd += ' font "DejaVuSans-Bold"\n'
                else:
                    label_cmd += '\n'"""
                label_cmd += '\n'
                fo.write(label_cmd)

        # We do not have anything to plot here, genes are drawn
        fo.write('plot 0 \n')

        fo.write('unset multiplot\n')

    with open(os.devnull, 'w') as dn:
        pop = Popen([gnuplot, script], stdout=dn, stderr=dn)
        pop.wait()

    os.unlink(wigs_file)
    os.unlink(disps_file)
    os.unlink(script)

    if pop.returncode != 0:
        if os.path.isfile(dest):
            os.unlink(dest)
        #sys.stderr.write('Gnuplot error: local QC displays might not be produced.\n')

    return pop.returncode
