#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 31/01/15

@author: Matthias Blum
"""

import math
import os
import sqlite3
import sys

from matplotlib import pyplot as plt


def init_db(dbtxt):
    """
    Create an SQlite database containing the globalQC/target for each entry in the database
    :param dbtxt:   File (tab-separated columns). Use: SELECT SECID, target_molecule_common, organism, TMRS, UMRS, `stampsQC_2.5`, stampsQC_5, stampsQC_10  FROM `ngs_sample`
    :return:
    """
    dbfile = os.path.join(os.path.dirname(__file__), 'db', 'targets.db')

    if os.path.isfile(dbfile):
        os.unlink(dbfile)

    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute("CREATE TABLE public_data (id INTEGER PRIMARY KEY AUTOINCREMENT, "
                "secid TEXT NOT NULL, target TEXT NOT NULL, "
                "organism TEXT NOT NULL, reads INTEGER NOT NULL, "
                "unique_reads INTEGER NOT NULL, "
                "stampsqc25 REAL NOT NULL, stampsqc5 REAL NOT NULL, stampsqc10 REAL NOT NULL)")

    commit = True
    with open(dbtxt) as f:
        for line in f:
            cols = line.rstrip().split('\t')

            try:
                cur.execute("INSERT INTO public_data (secid, target, organism, reads, "
                            "unique_reads, stampsqc25, stampsqc5, stampsqc10) "
                            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                            (cols[0], cols[1], cols[2], int(cols[3]), int(cols[4]),
                             float(cols[5]), float(cols[6]), float(cols[7])))
            except Exception as e:
                sys.stderr.write('Error (libtargets.init_db): {} {}\n'.format(type(e), e))
                commit = False
                break

    if commit:
        con.commit()

    cur.close()
    con.close()


def get_data(target, organism):
    """
    Return the public data results associated to a given target molecule
    :param target:      Target molecule
    :param organism:    Organism (not genome!)
    :return:            The public data, and the target as stored in the database (h3k4me3 -> H3K4me3)
    """
    dbfile = os.path.join(os.path.dirname(__file__), 'db', 'targets.db')
    data = []

    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    try:
        cur.execute("SELECT secid, target, organism, reads, unique_reads, stampsqc25, stampsqc5, stampsqc10 "
                    "FROM public_data WHERE organism=? AND target=? COLLATE NOCASE", (organism, target,))
        res = cur.fetchall()
    except Exception as e:
        sys.stderr.write('Error (libtargets.find_target): {} {}\n'.format(type(e), e))
    else:
        for row in res:
            target = row[1]
            reads = row[3]
            qcvals = {
                2.5: row[5],
                5: row[6],
                10: row[7],
            }
            data.append([reads, qcvals])

    cur.close()
    con.close()
    return data, target


def plot(data, disp, reads, qcval, quartiles, dest, public_label='', label=''):
    """
    Plot the number of reads against the QC values for the public data and the submitted profile
    :param data:            Public data (list: [reads, qcvals]
    :param disp:            Dispersion (2.5, 5, or 10)
    :param reads:           Input dataset's reads
    :param qcval:           Input dataset's qcval for a given disp
    :param quartiles:       Current DB quartiles for a given dispersion
    :param dest:            Outout path for the image
    :param public_label:    Label given to public data
    :return:
    """
    plt.figure()
    fig, ax = plt.subplots()

    public_reads = []
    public_qcvals = []

    for ds_reads, ds_qcvals in data:
        # We want the number of reads in millions
        public_reads.append(ds_reads / 1000000)
        # Public QCvals are not yet in log2 (but input dataset's are)
        public_qcvals.append(math.log(ds_qcvals, 2))

    # Plot public data
    ax.scatter(public_reads, public_qcvals, marker='o', color='#749BBD', label=public_label, s=25, alpha=0.5)

    # Plot user data
    ax.scatter(reads/1000000, qcval, marker='o', color='#D91E18', label=label, s=50, alpha=1)

    ymin = min([min(public_qcvals), qcval])
    ymax = max([max(public_qcvals), qcval])

    # Force X/Y limits
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(0, xmax)
    ax.set_ylim(ymin-0.5, ymax+0.5)  # +/- 0.5 to no have the min/max values directly one the plot's bottom/top

    # Draw the quartiles/stamp letters
    prev = ymin
    if ymin < quartiles[25] < ymax:
        ax.axhline(quartiles[25], color='black', lw=0.7, ls='dashed')

        # Write the stamp letter bellow the just drawn quartile if there is enough space
        if quartiles[25] - prev >= 0.5:
            text_y = prev + (quartiles[25] - prev) / 2
            ax.text(xmax+1, text_y, 'D', horizontalalignment='left', verticalalignment='center',
                    size='large', color='#e57c3c', weight='bold')

        prev = quartiles[25]

    if ymin < quartiles[50] < ymax:
        ax.axhline(quartiles[50], color='black', lw=0.7, ls='dashed')

        if quartiles[50] - prev >= 0.5:
            text_y = prev + (quartiles[50] - prev) / 2
            ax.text(xmax+1, text_y, 'C', horizontalalignment='left', verticalalignment='center',
                    size='large', color='#e57c3c', weight='bold')

        prev = quartiles[50]

    if ymin < quartiles[75] < ymax:
        ax.axhline(quartiles[75], color='black', lw=0.7, ls='dashed')

        if quartiles[75] - prev >= 0.5:
            text_y = prev + (quartiles[75] - prev) / 2
            ax.text(xmax+1, text_y, 'B', horizontalalignment='left', verticalalignment='center',
                    size='large', color='#e57c3c', weight='bold')

        prev = quartiles[75]

    # Write the stamp letter above the last drawn quartile if there is enough space
    if ymax - prev >= 0.5:
        text_y = prev + (ymax - prev) / 2

        if prev == ymin:
            c = 'D'
        elif prev == quartiles[25]:
            c = 'C'
        elif prev == quartiles[50]:
            c = 'B'
        else:
            c = 'A'

        ax.text(xmax+1, text_y, c, horizontalalignment='left', verticalalignment='center',
                size='large', color='#e57c3c', weight='bold')

    # Hide ticks on the top and on the right of the plot (they are useless)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if public_label:
        # Add the legend
        ax.legend(loc='lower right', shadow=False, scatterpoints=1, fontsize='small')

    plt.xlabel('Mapped reads (millions)', fontsize='large')
    plt.ylabel('log2(qc-stamp) dCRCI < {}'.format(disp), fontsize='large')
    plt.title('QC-stamps ({}%)'.format(disp), fontsize='x-large')

    plt.tight_layout()
    plt.savefig(dest)


def list_targets():
    """
    # Print the XML-formated list of target molecules for Galaxy
    :return:
    """
    targets = []

    dbfile = os.path.join(os.path.dirname(__file__), 'db', 'targets.db')
    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute('SELECT DISTINCT target COLLATE NOCASE from public_data ORDER BY target')
    for row in cur.fetchall():
        targets.append(row[0])
    cur.close()
    con.close()

    for target in sorted(targets, key=lambda s: s.lower()):
        print('<option value="{0}">{0}</option>'.format(target))


if __name__ == '__main__':
    init_db(sys.argv[1])
    list_targets()