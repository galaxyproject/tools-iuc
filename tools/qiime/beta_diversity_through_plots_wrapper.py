#!/usr/bin/env python

# create a html output for beta_diverity_through_plots.py
# to replace beta_diversity_through_plot_html.py in qiime 1.8
# Leroi Laura laura.leroi@ifremer.fr

import optparse
import os
import sys
import subprocess

# Create html template
def html_maker(metrics):

    tag=""
    plots = metrics.split(',')
    for i in range(0, len(plots)):
        tag = tag+"""\t\t<option value='%s_emperor_pcoa_plot/index.html'>%s</option>\n""" % (plots[i], plots[i])
    template_html= """
<!doctype html>
<html lang="en">
<head>
        <meta charset="utf-8">
        <script  type="text/javascript" src="/galaxy/static/scripts/libs/jquery/jquery.js"></script>
</head>
<body>
    <label for="outputs">Select a beta-diversity metric :</label>
        <select id="outputs" name="outputs">
%s
        </select>
    <br/>
    <br/>
        <iframe id="plot" name="plot" style="position:absolute; width: 90%%; height: 90%%;" src="" frameborder="0"></iframe>
        <script>
        $(document).ready(function(){
                var str = "";
                $( "#outputs option:selected" ).each(function() {
                        str += $( this ).val() + " ";
                        $("#plot").attr("src", str);
                });

        });
        $( "#outputs" ).change(function () {
                var str = "";
                $( "#outputs option:selected" ).each(function() {
                        str += $( this ).val() + " ";
                        $("#plot").attr("src", str);
                });
        }).change();
</script>
</body>
</html>    """ % (tag)
    return template_html


def html_matrix_maker(metrics):

    tag=""
    plots = metrics.split(',')
    for i in range(0, len(plots)):
        tag = tag+"""\t\t<a href='%s_dm.txt'>%s_dm</a><br/>\n""" % (plots[i], plots[i])

    template_html= """
<!doctype html>
<html lang="en">
<body>
        <label>Metrics matrix : </label><br/>
%s
</body>
</html> """ % (tag)
    return template_html


def main():

    #On importe l environnemnent envqiime
    pipe = subprocess.Popen(". /home12/caparmor/bioinfo/GALAXY_DEV/ourenv/envqiime; python -c 'import os; print \"newenv = %r\" % os.environ'", stdout=subprocess.PIPE, shell=True)
    exec(pipe.communicate()[0])
    os.environ.update(newenv)

    #New parser
    parser = optparse.OptionParser()
    #beta diversity options
    parser.add_option('--metrics',dest='metrics')
    parser.add_option('--biom',dest='biom')
    parser.add_option('--mapping',dest='mapping')
    parser.add_option('--tree',dest='tree')
    parser.add_option('--depth',dest='depth')
    parser.add_option('--masterhtml',dest='masterhtml')
    parser.add_option('--outputfolder',dest='outputfolder')
    parser.add_option('--matrixhtml',dest='matrixhtml')

    (options,args) = parser.parse_args()

    #New hmtl template
    master = open('./master.html', 'w')
    master.write(html_maker(options.metrics))
    master.close()
    #On lance le fichier avec les options precisees pour l'affichage html
    os.system("mv master.html %s" % options.masterhtml)

    #parameter file, which specifies changes to the default behavior
    parameters = open('./parameters.txt','w')
    parameters.write("beta_diversity:metrics %s" % options.metrics)
    parameters.close()

    if options.depth:
        if options.tree :
            os.system("beta_diversity_through_plots.py -i %s -m %s -t %s -e %s -o %s -p ./parameters.txt >& file_log.txt" % (options.biom, options.mapping, options.tree, options.depth, options.outputfolder))
        else :
            os.system("beta_diversity_through_plots.py -i %s -m %s -e %s -o %s -p ./parameters.txt >& file_log.txt" % (options.biom, options.mapping, options.depth, options.outputfolder))
    else:
        if options.tree :
            os.system("beta_diversity_through_plots.py -i %s -m %s -t %s -o %s -p ./parameters.txt >& file_log.txt" % (options.biom, options.mapping, options.tree, options.outputfolder))
        else :
            os.system("beta_diversity_through_plots.py -i %s -m %s -o %s -p ./parameters.txt >& file_log.txt" % (options.biom, options.mapping, options.outputfolder))

    os.system("mkdir %s/txt" % (options.outputfolder))
    os.system("mv %s/*_dm.txt %s/txt" %(options.outputfolder,options.outputfolder))
    os.system("compress_path.py -i %s/txt -o %s" % (options.outputfolder,options.matrixhtml))

if __name__=="__main__":
    main()

