#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 29/01/15

@author: Matthias Blum
"""


import libngs
import os

from datetime import datetime
from PIL import Image
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Frame, Table
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import cm


def frame(x, y, width, height, top=0, right=0, bottom=0, left=0):
    return Frame(x, y, width, height, topPadding=top, rightPadding=right, bottomPadding=bottom, leftPadding=left)


def parse_summary(summary):
    """
    Return a dictionary of the information found in a NGS-QC summary file
    :param summary:     Run's summary
    :return:
    """
    run = {
        'bg': None,             # Background threshold
        'bins': None,           # Number of bins
        'bins2': None,          # Bins without with a RCI greater than the background threshold
        'binsize': None,        # Size of each bin in bp
        'chrms': None,          # List of chromosomes
        'ci': None,             # Confidence interval
        'date': None,           # Run date
        'denqc_50_10': None,    # Density QCis
        'denqc_50_2.5': None,
        'denqc_50_5': None,
        'denqc_90_10': None,
        'denqc_90_2.5': None,
        'denqc_90_5': None,
        'genome':None,          # Genome assembly
        'infile': None,         # Input file name
        'mean_len': None,       # Average size of reads
        'nobgs': None,          # True if background subtraction was enabled
        'nodup': None,          # True if PCR duplicate reads ware removed
        'nreads': None,         # Number of mapped reads
        'nreplicates': None,    # Number of replicates
        'nureads': None,        # Number of unique reads
        'rep': None,            # Current replicate
        'sampling': None,       # Sampling performed (not in %)
        'simqc_10': None,       # Similarity QCis
        'simqc_2.5': None,
        'simqc_5': None,
        'target': None,         # Target molecule
        'version': None,        # NGS-QC version
    }

    with open(summary) as f:
        for i, line in enumerate(f):
            if i == 0:
                cols = line.rstrip().split()
                run['version'] = cols[5]
                run['date'] = cols[6] + ' ' + cols[7]
            else:
                field, value = line.rstrip().split('\t', 1)

                if field in ('Threshold', 'Threshold (100%)'):
                    run['bg'] = int(value)
                elif field == 'Bins':
                    run['bins'] = int(value)
                elif field == 'Bins w/o background':
                    run['bins2'] = int(value)
                elif field == 'Windows size':
                    run['binsize'] = int(value)
                elif field == 'Chromosomes':
                    run['chrms'] = sorted(value.split(','))
                elif field == 'Conf. interval':
                    run['ci'] = float(value)
                elif field == 'DenQC_s90_2.5pc(%)':
                    run['denqc_90_2.5'] = float(value)
                elif field == 'DenQC_s50_2.5pc(%)':
                    run['denqc_50_2.5'] = float(value)
                elif field == 'DenQC_s90_5pc(%)':
                    run['denqc_90_5'] = float(value)
                elif field == 'DenQC_s50_5pc(%)':
                    run['denqc_50_5'] = float(value)
                elif field == 'DenQC_s90_10pc(%)':
                    run['denqc_90_10'] = float(value)
                elif field == 'DenQC_s50_10pc(%)':
                    run['denqc_50_10'] = float(value)
                elif field == 'Genome assembly':
                    run['genome'] = value
                elif field == 'Input file':
                    run['infile'] = value
                elif field == 'Reads size mean':
                    run['mean_len'] = float(value)
                elif field in ('Background subtraction', 'Background subtracted'):
                    run['nobgs'] = value == 'Off'
                if field in ('Clonal reads removal', 'Clonal reads removed'):
                    run['nodup'] = value == 'On'
                elif field == 'Reads':
                    run['nreads'] = int(value)
                elif field == 'Replicate':
                    run['rep'] = int(value.split('/')[0])
                    run['nreplicates'] = int(value.split('/')[-1])
                elif field == 'Unique reads':
                    run['nureads'] = int(value)
                elif field == 'Sampling':
                    run['sampling'] = sorted([float(i.strip()) for i in value.split(',')])
                elif field == 'SimQC_2.5pc (DenQC_s90/s50)':
                    value = float(value)

                    if value == float('inf'):
                        continue
                    else:
                        run['simqc_2.5'] = value
                elif field == 'SimQC_5pc (DenQC_s90/s50)':
                    value = float(value)

                    if value == float('inf'):
                        continue
                    else:
                        run['simqc_5'] = value
                elif field == 'SimQC_10pc (DenQC_s90/s50)':
                    value = float(value)

                    if value == float('inf'):
                        continue
                    else:
                        run['simqc_10'] = value
                elif field == 'Target molecule':
                    if len(value) == 0 or value == 'None':
                        run['target'] = None
                    else:
                        run['target'] = value

        for k, v in run.items():
            if v is None and k not in ('mean_len', 'target'):
                return None

        return run


def get_stamp(qcval, quartiles):
    """
    Return the stamp letter for the given QC value
    :param qcval:       QC val for a dispersion
    :param quartiles:   Quartile dictionary
    :return:            str (one letter)
    """
    stamp = 'N'

    if qcval is not None:
        if qcval > quartiles[75]:
            stamp = 'A'
        elif quartiles[75] >= qcval > quartiles[50]:
            stamp = 'B'
        elif quartiles[50] >= qcval > quartiles[25]:
            stamp = 'C'
        elif qcval <= quartiles[25]:
            stamp = 'D'

    return stamp


def insert_header(canvas, page_height, page_width):
    canvas.saveState()
    canvas.setFillColor('#396A92')

    y = page_height - 2.5*cm
    canvas.rect(0, y, page_width, y + 2.5*cm, fill=1, stroke=0)

    p = Paragraph('NGS-QC Generator',
                  ParagraphStyle('header', fontName='Enriqueta-Bold', fontSize=25, textColor='#FFFFFF'))

    p.wrapOn(canvas, page_width, 2.5*cm)
    p.drawOn(canvas, 1.5*cm, y + 2.5*cm / 2)
    canvas.restoreState()


def insert_footer(canvas, width, db_version=None, page=None):
    canvas.setLineWidth(0.25)
    canvas.line(1.5*cm, 1.5*cm, 1.5*cm+width, 1.5*cm)

    style = ParagraphStyle('footer', fontSize=8, alignment=TA_LEFT)
    p = Paragraph('Copyright &copy; {0} '
                  '<a href="http://www.ngs-qc.org">ngs-qc.org</a>'.format(datetime.now().year), style)
    frame(1.5*cm, 0, width, 1.5*cm, top=10, left=10).addFromList([p], canvas)

    style.alignment = TA_CENTER
    if db_version:
        p = Paragraph('Database v{0}'.format(db_version), style)
        frame(1.5*cm, 0, width, 1.5*cm, top=10).addFromList([p], canvas)
    elif page:
        p = Paragraph('Page {0}'.format(page), style)
        frame(1.5*cm, 0, width, 1.5*cm, top=10).addFromList([p], canvas)

    style.alignment = TA_RIGHT
    p = Paragraph('<a href="mailto:contact@ngs-qc.org">contact@ngs-qc.org</a>', style)
    frame(1.5*cm, 0, width, 1.5*cm, top=10, right=10).addFromList([p], canvas)


def ngsqc_report(run, dest, rundir, quartiles=None, dbversion=None, uniqbar=True):
    """
    Generate a PDF report of a NGS-QC run
    :param run:         Dictionary of the run info
    :param dest:        PDF output path
    :param rundir:      Directory containing all images
    :param quartiles:   Quartile dictionary
    :param dbversion:   NGS-QC version
    :param uniqbar      If True, display a color bar of the ratio of unique/total reads next to the stamp
    :return:
    """
    # Margin: 1.5 cm left/right

    organism = libngs.getorganism(run['genome'])
    if organism:
        genus, species = organism.split()
        organism = '{}. {}'.format(genus[0].upper(), species.lower())

    qcvals = {
        2.5: libngs.getqcval(run['denqc_50_2.5'], run['simqc_2.5']),
        5: libngs.getqcval(run['denqc_50_5'], run['simqc_5']),
        10: libngs.getqcval(run['denqc_50_10'], run['simqc_10'])
    }

    if quartiles is not None:
        stamps = {}
        for disp in (2.5, 5, 10):
            stamps[disp] = get_stamp(qcvals[disp], quartiles[disp])
        stamp = stamps[2.5] + stamps[5] + stamps[10]
    else:
        stamp = None

    pdfmetrics.registerFont(TTFont('Enriqueta-Bold',
                                   os.path.join(os.path.dirname(__file__), 'fonts', 'Enriqueta-Bold.ttf')))
    pdfmetrics.registerFont(TTFont('Arial-Bold',
                                   os.path.join(os.path.dirname(__file__), 'fonts', 'arial-bold.ttf')))
    pdfmetrics.registerFont(TTFont('Arial',
                                   os.path.join(os.path.dirname(__file__), 'fonts', 'arial.ttf')))
    pdfmetrics.registerFont(TTFont('Capture',
                                   os.path.join(os.path.dirname(__file__), 'fonts', 'capture_it.ttf')))

    cell_style = ParagraphStyle('cell',
                                fontSize=10, textColor='#000000')
    title_style = ParagraphStyle('title',
                                 fontSize=16, leading=18, spaceBefore=4, spaceAfter=4)
    subtitle_style = ParagraphStyle('subtitle',
                                    fontSize=10, leading=12, spaceBefore=4, spaceAfter=4, alignment=TA_RIGHT)
    normal_style = ParagraphStyle('normal',
                                  fontSize=10, leading=12, spaceBefore=4, spaceAfter=4, alignment=TA_JUSTIFY)
    legend_style = ParagraphStyle('legend',
                                  fontSize=8, leading=12, spaceBefore=4, spaceAfter=4, alignment=TA_JUSTIFY)

    canvas = Canvas(dest, pagesize=A4)
    (page_width, page_height) = defaultPageSize
    p_width = page_width - 1.5 * 2 * cm

    canvas.setTitle("NGS-QC Generator Report")
    canvas.setAuthor("Hinrich Gronemeyer group")

    canvas.setLineWidth(0.25)

    ######################
    #   HEADER
    ######################
    insert_header(canvas, page_height, page_width)

    ######################
    #   DISP PLOT
    ######################
    y = page_height - 420
    path = os.path.join(rundir, 'pc_s50_s70_s90_replicate_{0}.png'.format(run['rep']))

    # Images generated by Gnuplots are in "P" mode. We want them in RGBA to avoid a Pillow warning
    im = Image.open(path).convert('RGBA')
    im.save(path)
    canvas.drawImage(path, page_width/2, y, width=p_width/2, height=p_width/2)

    # Legend
    y -= 100
    p = Paragraph('<b>Effect of random sampling on the profile</b>. '
                  'This figure illustrates the influence '
                  'of the random sampling subsets (90%: black; 70%: blue; 50%: red) '
                  'on the recovered read count Intensity (recRCI) per bin. '
                  'The dark-green vertical line '
                  'represents the background threshold ({} RCI).'.format(run['bg']), legend_style)
    frame(page_width/2, y, p_width/2, 100, left=5).addFromList([p], canvas)

    ######################
    # TITLE / SUBTITLE
    ######################
    y = page_height - 2.5 * cm - 20
    p = Paragraph('<b>Data Quality Report</b>', title_style)
    p.wrapOn(canvas, p_width, 20)
    p.drawOn(canvas, 1.5*cm, y)

    p = Paragraph(run['date'], subtitle_style)
    p.wrapOn(canvas, p_width, 20)
    p.drawOn(canvas, 1.5*cm, y)

    ######################
    # LINE
    ######################
    y -= 5
    canvas.line(1.5*cm, y, page_width - 1.5 * cm, y)

    ######################
    # FILE NAME
    ######################
    y -= 35

    p = Paragraph('File name: ' + run['infile'], normal_style)
    frame(1.5*cm, y, p_width, 30, left=5).addFromList([p], canvas)

    ######################
    # CHIP-SEQ STAMP
    ######################
    if stamp is not None:
        y = page_height - 125
        canvas.saveState()
        if run['nobgs']:
            canvas.setFillColor(colors.fidred)
            canvas.setStrokeColor(colors.fidred)
        else:
            canvas.setFillColor(colors.green)
            canvas.setStrokeColor(colors.green)

        canvas.rotate(-7)
        canvas.setFont('Arial-Bold', 14)
        canvas.drawCentredString(page_width - 1.5*cm - 185, y + 45, 'Global QC certification')
        canvas.setFont('Capture', 45)
        canvas.drawCentredString(page_width - 230, y + 5, stamp)
        canvas.setLineWidth(4)
        canvas.rect(page_width - 320, y, 180, 60, stroke=True, fill=False)

        if uniqbar:
            ######################
            # UNIQUE READS BAR
            ######################
            unique_ratio = float(run['nureads']) / run['nreads']
            # canvas.saveState()
            canvas.setLineWidth(1)
            canvas.setStrokeColor('#000000')
            canvas.setFillColor('#CE68F5')
            p = canvas.beginPath()

            xbar = page_width / 2 - 5
            ybar = page_height - 120
            p.rect(xbar, ybar, 7, 30)
            canvas.drawPath(p, fill=0, stroke=1)

            # Green part
            canvas.setFillColor(colors.green)
            p = canvas.beginPath()
            p.rect(xbar + 0.5, ybar+0.5, 6, 30*unique_ratio)
            canvas.drawPath(p, fill=1, stroke=0)

            # Red part
            canvas.setFillColor(colors.red)
            p = canvas.beginPath()
            p.rect(xbar + 0.5, ybar+30*unique_ratio, 6, 30*(1-unique_ratio))
            canvas.drawPath(p, fill=1, stroke=0)

            canvas.rotate(90)
            canvas.setFillColor(colors.black)
            canvas.setFont('Arial', 7)
            canvas.drawCentredString(737, -310, 'URs (%)')
        canvas.restoreState()

    ######################
    # DATA SET INFO
    ######################
    organism_str = ' ({})'.format(organism)
    unique_ratio = float(run['nureads']) / run['nreads']

    if unique_ratio >= 0.5:
        p = Paragraph('{0} ({1:.2f}%)'.format(libngs.num_format(run['nureads']), unique_ratio*100), cell_style)
    else:
        p = Paragraph("%s (<font color=fidred>%.2f%%</font>)" % (libngs.num_format(run['nureads']),
                                                                 unique_ratio*100), cell_style)

    data = [['Dataset informations', None],
            ['Total reads', libngs.num_format(run['nreads'])],
            ['Unique reads (URs)', p],
            #["Reads's size mean (bp)", run['mean_len']],
            ['Genome assembly', run['genome'] + organism_str],
            ['Target molecule', run['target']]]

    t = Table(data, style=[('GRID', (0, 0), (-1, -1), 0.5, colors.black),
                           ('FONTSIZE', (0, 0), (-1, -1), 10),
                           ('BACKGROUND', (0, 0), (1, 0), colors.grey),
                           ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
                           #('TEXTCOLOR', (1, 3), (1, 3), colors.orange),
                           ('SPAN', (0, 0), (1, 0)),
                           ('SPAN', (0, 0), (1, 0))], colWidths=[p_width/4, p_width/4])

    y = page_height - 225
    t.wrapOn(canvas, p_width / 2, 10)
    t.drawOn(canvas, 1.5*cm, y)

    ######################
    # QC PARAMETERS
    ######################
    y -= 135
    data = [['QC parameters', ''],
            ['Sampling percentages', ', '.join(str(int(e*100)) for e in run['sampling'])],
            ['Windows size (bp)', run['binsize']],
            ['Replicate number', '{0}/{1}'.format(run['rep'], run['nreplicates'])],
            #['Sampled strand(s)', run['strandmode']],
            ['Referenced chromosomes', len(run['chrms']) if run['chrms'] else None],
            ['Background subtraction', 'Off' if run['nobgs'] else 'On'],
            ['Clonal reads removal', 'On' if run['nodup'] else 'Off']]

    t = Table(data, style=[('GRID', (0, 0), (-1, -1), 0.5, colors.black),
                           ('BACKGROUND', (0, 0), (1, 0), colors.darkgrey),
                           ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
                           ('SPAN', (0, 0), (1, 0))], colWidths=[p_width/4, p_width/4])

    t.wrapOn(canvas, p_width / 2, 50)
    t.drawOn(canvas, 1.5*cm, y)

    ######################
    # INDICATORS
    ######################
    y -= 100

    considered_reads = run['nureads'] if run['nodup'] else run['nreads']
    data = [['Results', None, None],
            ['Considered reads*', None, libngs.num_format(considered_reads)],
            ['QC values\ndenQC (50%) / simQC', '2.5%', '{:.3f} / {:.3f}'.format(run['denqc_50_2.5'],
                                                                                run['simqc_2.5'])],
            [None, '5%', '{:.3f} / {:.3f}'.format(run['denqc_50_5'], run['simqc_5'])],
            [None, '10%', '{:.3f} / {:.3f}'.format(run['denqc_50_10'], run['simqc_10'])]]

    t = Table(data, style=[
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
        ('BACKGROUND', (0, 0), (2, 0), colors.lightgrey),
        ('SPAN', (0, 0), (2, 0)),
        ('SPAN', (0, 1), (1, 1)),
        ('SPAN', (0, 2), (0, 4)),
        ('VALIGN', (0, 2), (0, 2), 'MIDDLE'),
        ('FONTSIZE', (0, 2), (0, 2), 9)
    ], colWidths=[0.75*p_width/4, 0.25*p_width/4, p_width/4])

    t.wrapOn(canvas, p_width / 2, 50)
    t.drawOn(canvas, 1.5*cm, y)

    y -= 20
    p = Paragraph('* Reads taken into account to compute the QC indicators.', legend_style)
    frame(1.5*cm, y, p_width/2, 20, left=5, top=5).addFromList([p], canvas)

    ######################
    # LINE
    ######################
    y -= 5
    canvas.line(page_width/4, y, 0.75*page_width, y)

    ######################
    # LOCAL QC
    ######################
    y -= 40
    p = Paragraph('<b>Read count intensity profile illustrated '
                  'in the context of its corresponding local QC indicators (heatmap)</b>. '
                  'On the upper figure, genes are represented by green-colored rectangles. '
                  'Overlapping genes are represented by a deeper green and TSS are displayed as a small dark bar. '
                  'The lower figure is a zoom of the center of the first figure.', legend_style)
    frame(1.5*cm, y, p_width, 42, top=5).addFromList([p], canvas)

    imgs = {}
    string = 'localqc.rep{}.region'.format(run['rep'])
    localqcdir = os.path.join(rundir, 'localqc')

    # Index the regions by their original order
    for f in os.listdir(localqcdir):
        if f.startswith(string):
            idx = int(f.replace(string, '')[:-4])
            imgs[idx] = os.path.join(localqcdir, f)

    imgs_width = 0.95 * p_width
    margin_left = 0.05 * p_width

    colorbar = os.path.join(localqcdir, 'colorbar.png')
    try:
        im = Image.open(colorbar)
    except IOError:
        insert_colorbar = False
    else:
        insert_colorbar = True
        colb_w, colb_h = im.size
        del im

    num_page = 1

    for idx, path in sorted(imgs.items()):
        if idx == 3:
            # We placed the first region (low+high resolution) on the first page
            # Now, move to the next page
            insert_footer(canvas, p_width, db_version=dbversion)
            canvas.showPage()
            insert_header(canvas, page_height, page_width)
            y = page_height - 2.5*cm - 10
            num_page += 1

        if num_page == 2 and idx % 2 != 0:
            # We only want high resolution images on page 2
            # These regions have a even number in the list
            continue

        im = Image.open(path)
        img_w, img_h = im.size
        del im

        if insert_colorbar:
            img_h *= imgs_width / img_w
        else:
            img_h *= p_width / img_w

        y -= img_h

        if y <= 1.5*cm:
            num_page += 1
            if num_page == 3:
                break

        if insert_colorbar:
            canvas.drawImage(path, 1.5*cm + margin_left, y, width=imgs_width, height=img_h)
            colb_w *= img_h / colb_h
            canvas.drawImage(colorbar, 1.5*cm, y, width=colb_w, height=img_h)
            insert_colorbar = False
        else:
            canvas.drawImage(path, 1.5*cm, y, width=p_width, height=img_h)

        y -= 10

    # Last page's footer
    insert_footer(canvas, p_width, db_version=dbversion)

    #####################
    # TARGET
    #####################
    if run['target'] is not None:
        imgs = [os.path.join(rundir, 'public.{}.png'.format(disp)) for disp in (2.5, 5, 10)]
        if len([path for path in imgs if os.path.isfile(path)]) == 3:
            canvas.showPage()
            insert_header(canvas, page_height, page_width)
            y = page_height - 2.5*cm - 10

            for path in imgs:
                im = Image.open(path)
                img_w, img_h = im.size
                img_h *= 0.38
                img_w *= 0.38
                del im
                x = (page_width - img_w) / 2
                y -= img_h
                canvas.drawImage(path, x, y, width=img_w, height=img_h)

            # Legend
            y -= 35
            p = Paragraph('<b>Comparison of the assessed quality grade</b> with those computed '
                          'for publicly available datasets currently hosted '
                          'in the NGS-QC database that correspond to the same antibody target.', legend_style)
            frame(x+20, y, width=img_w-25, height=40, left=5).addFromList([p], canvas)

            insert_footer(canvas, p_width, db_version=dbversion)
    canvas.save()
    return dest


def init(summary, dest, src, quartiles=None, dbversion=None):
    """
    Parse a summary file and generate a NGS-QC report out of it
    :param summary:         Path to summary file
    :param dest:            Output PDF file
    :param src:             NGSQC's working directory (where images are)
    :param quartiles:       Dictionary of quartiles
    :param dbversion:       Version of the NGS-QC DB
    :return:
    """
    run = parse_summary(summary)

    if run is None:
        return None

    return ngsqc_report(run, dest, src, quartiles, dbversion)


def insert_colorbar_gv(canvas, top, img, img_h):
    """
    Insert the dispersion colorbar
    :param canvas:  Canvas to write on
    :param top:     Page top
    :param img:     Image path
    :param img_h:  Image desired height (since regions and colorbar have the same height)
    :return:
    """

    im = Image.open(img)
    width, height = im.size
    img_w = width * img_h / height
    del im
    canvas.drawImage(img, 25, top - 75 - img_h, width=img_w, height=img_h)


def gv_report(dest, tmp_dir):
    canvas = Canvas(dest, pagesize=A4)
    (page_width, page_height) = defaultPageSize
    top = page_height
    p_width = page_width - 1.5 * 2 * cm
    canvas.setTitle("LocalQCs Report")
    canvas.setAuthor("Hinrich Gronemeyer group")

    pdfmetrics.registerFont(TTFont('Enriqueta-Bold',
                                   os.path.join(os.path.dirname(__file__), 'fonts', 'Enriqueta-Bold.ttf')))

    ####################################
    # Header
    ####################################
    insert_header(canvas, page_height, page_width)

    ####################################
    # Footer
    ####################################
    num_page = 1
    insert_footer(canvas, p_width, page=num_page)

    imgs = {}
    for f in os.listdir(tmp_dir):
        if f.endswith('.png') and f.startswith('region'):
            idx = int(f.replace('region', '').replace('.png', ''))
            imgs[idx] = os.path.join(tmp_dir, f)

    colorbar = os.path.join(tmp_dir, 'colorbar.png')

    ####################################
    # Place images
    ####################################
    y_pos = top - 75

    for idx, img in sorted(imgs.items()):
        im = Image.open(img)
        # Image's real size
        img_w, img_h = im.size
        del im

        # Resize image to to fit to the page keeping the same ratio
        img_h *= p_width / img_w
        y_pos -= img_h

        if y_pos <= 1.5*cm:
            # Not enough space: start a new page
            insert_colorbar_gv(canvas, top, colorbar, img_h)
            canvas.showPage()
            num_page += 1
            insert_header(canvas, page_height, page_width)
            insert_footer(canvas, p_width, page=num_page)
            y_pos = top - 75 - img_h

        canvas.drawImage(img, 2 * cm, y_pos, width=p_width, height=img_h)
        y_pos -= 15  # Margin bottom

    insert_colorbar_gv(canvas, top, colorbar, img_h)
    canvas.save()
