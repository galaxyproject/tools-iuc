#!/usr/bin/env python3
import argparse
import csv
import logging
import sys


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    options = _set_options()
    data, headers = _read_map_file(options.map)
    html = _print_html(data, headers, options.out)
    index_file = options.out + '/index.html'
    fh = open(index_file, mode='w')
    fh.write(html)
    fh.close()


def _get_google_script_headers(data, headers, out_dir):
    html = '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>' + "\n"
    html += '<script type="text/javascript">' + "\n"
    html += 'google.charts.load(\'current\', {\'packages\':[\'table\']});' + "\n"
    chart_names, java_scripts = _get_google_js(data, headers, out_dir)
    for i in range(0, len(chart_names)):
        html += 'google.charts.setOnLoadCallback(' + chart_names[i].replace('-', '_') + ');' + "\n"
        html += 'function ' + chart_names[i].replace('-', '_') + '() {' + "\n"
        html += java_scripts[i] + "\n"
        html += '}' + "\n"
    html += '</script>' + "\n"
    return html


def _get_google_js(data, headers, out_dir):
    java_scripts = []
    chart_names = []
    for cdd in data:
        chart_names.append(cdd['cdd_id'] + '_' + cdd['description'])
        js = 'var data = new google.visualization.DataTable();' + "\n"
        mat, head = _parse_csv(out_dir + '/' + cdd['cluster_nb_reads_files'])
        for el in head:
            if el == '#OTU_name':
                js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
            elif el == 'taxonomy':
                js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
            elif el == 'contigs_list' or el == 'seq_list':
                js += 'data.addColumn(\'string\', \'' + el + '\');' + "\n"
            else:
                js += 'data.addColumn(\'number\', \'' + el + '\');' + "\n"
        js += 'data.addRows([' + "\n"
        for j in range(0, len(mat)):
            js += '[\'' + mat[j][head[0]] + '\''
            for i in range(1, len(head) - 2):
                js += ',' + mat[j][head[i]]
            js += ',\'' + mat[j][head[len(head) - 2]] + '\''
            js += ',\'' + mat[j][head[len(head) - 1]] + '\''
            js += ']'
            if j != (len(mat) - 1):
                js += ','
            js += "\n"
        js += ']);' + "\n"
        js += 'var table = new google.visualization.Table(document.getElementById(\'' + (cdd['cdd_id'] + '_' + cdd['description']).replace('-', '_') + '_div' + '\'));' + "\n"
        js += 'table.draw(data, {showRowNumber: false, width: \'70%\', height: \'70%\'});' + "\n"
        java_scripts.append(js)
    return chart_names, java_scripts


def _parse_csv(file):
    fh = open(file)
    reader = csv.reader(fh, delimiter="\t")
    data = list(reader)
    headers = data[0]
    matrix = []
    for i in range(1, len(data)):
        dict = {}
        for j in range(0, len(data[i])):
            if data[i][j] == '':
                dict[headers[j]] = None
            elif data[i][j] == 'null':
                dict[headers[j]] = None
            else:
                dict[headers[j]] = data[i][j]
        matrix.append(dict)
    return matrix, headers


def _print_html(data, headers, out_dir):
    html = '<html>' + "\n"
    html += '<head>' + "\n"
    html += '<title>' + 'rps2tree' + '</title>'
    html += _get_google_script_headers(data, headers, out_dir)
    html += '</head>' + "\n"
    html += '<div style="text-align:center">' + "\n"
    html += '<h1 align=center>rps2tree</h1>' + "\n"
    html += '<body>' + "\n"
    html += _print_data(data, headers)
    html += '</body>' + "\n"
    html += '</div>' + "\n"
    html += '</html>' + "\n"
    return html


def _print_data(data, headers):
    html = ''
    for cdd in data:
        html += '<h2>' + cdd['cdd_id'] + ' ' + cdd['description'] + '</h2>' + "\n"
        html += '<p>' + cdd['full_description'] + '</br>' + '</p>' + "\n"
        html += '<div id="' + (cdd['cdd_id'] + '_' + cdd['description']).replace('-', '_') + '_div' + '"></div>' + "\n"
        html += '</br>' + "\n"
        html += '</br>' + "\n"
        html += '<img src=' + cdd['tree_files'] + ' href="' + cdd['tree_files'] + '">' + "\n"
        html += '</br>' + "\n"
        html += '<a href="' + cdd['align_files'] + '">' + cdd['align_files'] + '</a>' + "\n"
        html += '</br>' + "\n"
        html += '<a href="' + cdd['cluster_files'] + '">' + cdd['cluster_files'] + '</a>' + "\n"
        html += '</br>' + "\n"
        html += '<a href="' + cdd['cluster_nb_reads_files'] + '">' + cdd['cluster_nb_reads_files'] + '</a>' + "\n"
        html += '</br>' + "\n"
        html += '<a href="' + cdd['pairwise_files'] + '">' + cdd['pairwise_files'] + '</a>' + "\n"
        html += '</br>' + "\n"
        html += '</br>' + "\n"
        html += '<hr>' + "\n"
    return html


def _read_map_file(file):
    reader = csv.reader(file, delimiter="\t")
    data = list(reader)
    headers = data[0]
    headers[0] = headers[0][1:]
    map_obj = []
    for i in range(1, len(data)):
        dict = {}
        if len(data[i]) != len(headers):
            sys.exit('line and headers not the same length.')
        for j in range(0, len(headers)):
            dict[headers[j]] = data[i][j]
        map_obj.append(dict)
    return map_obj, headers


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--map', help='The map file produced by rps2tree.pl script.', action='store', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o', '--out', help='The title for the HTML page.', action='store', type=str, default='./')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
