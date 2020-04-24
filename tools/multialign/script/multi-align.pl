#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Getopt::Long qw(VersionMessage);
use File::Find qw(finddepth);
our $VERSION = '1.0';

if(@ARGV)
{
    my ($ref, @files, $fastq_n, $dir, $fastq, $help, $html_out, $build_index, $ma, $mis, $max_procs);
    # Uncomment the following line to run manually
    #my $max_procs = 1; 
    $ma = 0;
    my @tab = ();
    
    GetOptions (
    "fastq=s" => \$fastq,
    "dir=s" => \$dir,
    "ref=s" => \$ref,
    "ma=i" => \$ma,
    "html:s" => \$html_out,
    "help" => \$help,
    "mis=i" => \$mis,
    "threads:1" => \$max_procs,
    "version" => sub { VersionMessage(0); }
    );
    
    # Uncomment the following line to run manually
    #mkdir $dir or die "Cannot create $dir: $!\n";
    my $file = $dir.'/report.txt';
    open my $report, '>', $file or die "Cannot open $file $!\n";
    
    my @references =(); my @ref_names = ();
    my $i = -1; my $string = '';
    open (my $in, $ref) or die "Cannot open fasta file: $!\n";
    
    while(<$in>) 
    {
        if ($_ =~ /^>/)
        {
            push @ref_names, $1 if ($_ =~ /^>(.*?)\s|\n/);
            $references[$i] = $references[$i].$string if $string ne '';
            push @references, $_."\n";
            $string = '';
            $i++;
        }
        else 
        {
            chomp $_;
            $string = $string.$_;
        }
    }
    close $in;
    $references[$i] = $references[$i].$string if $string ne '';
    
    my ($name,$path,$suffix) = fileparse($fastq,'.fastq','.ref', '.fq','.dat');
    my $Gviz = $dir; 
    my $j = 0;
    
    push (@tab, dirname(__FILE__)."/align_count_stranded.pl");
    push (@tab, dirname(__FILE__)."/G_windows.pl");
    push (@tab, dirname(__FILE__)."/bootstrap/");
    
    foreach my $re (@references)
    {
        my $sam = $dir.'/'.$name.'_'.$ref_names[$j].'_sorted.sam';
        open (my $tmp_re, ">$ref_names[$j]") or die "Cannot open reference file: $!\n";
        print $tmp_re $re;
        close $tmp_re;
    
        my $fileAlignCountStrandedLocated = $tab[0];
        chomp $fileAlignCountStrandedLocated;
        my $fileGwindowsLocated= $tab[1];
        chomp $fileGwindowsLocated;
    
        #Construct index
        system("bwa index '$ref_names[$j]'");
        
        #Mapping
        system("$fileAlignCountStrandedLocated --fin $fastq --ref $ref_names[$j] --dir $dir --mis $mis --ma $ma  --noD --rep --p $max_procs  --aligner BWA 2>>$file");
        #print STDERR ("INFO COMMANDE: "."bash -c '$fileAlignCountStrandedLocated --fin $fastq --ref $ref_names[$j] --dir $dir --mis $mis --ma $ma  --noD --rep --p $max_procs  --aligner BWA'"."\n");
        
        system("$fileGwindowsLocated  --sam $sam --dir $Gviz  --win auto  --p $max_procs --ma $ma 2>>$file");
        #print STDERR ("INFO COMMANDE: "."bash -c '$fileGwindowsLocated  --sam $sam --dir $Gviz  --win auto  --p $max_procs --ma $ma'"."\n");
   
        #cleaning index
        system("rm '$ref_names[$j]'*");
        $j++;
    }
    
    $name = $name."_";
    
    # Collect counts
    my $folderGwindowsLocated = $tab[2];
    chomp $folderGwindowsLocated;
    
    system("cp -R $folderGwindowsLocated/* $dir");
    #print STDERR ("INFO COMMANDE : "."bash -c 'cp -R $folderGwindowsLocated/* $dir'"."\n");
    
    system("cat $dir/*_reads_counts.txt > tmp_count && rm $dir/*_reads_counts.txt");
    #print STDERR ("INFO COMMANDE : "."bash -c 'cat $dir/*_reads_counts.txt > tmp_count && rm $dir/*_reads_counts.txt'"."\n");
    my %raw_count = ();
    
    open (my $tc, 'tmp_count') or die "Cannot open tmp_count: $!\n";
    my @tmp = <$tc>;
    
    for (my $i = 0; $i <= $#tmp; $i++)
    {
        if ($i % 2 == 1){
            chomp $tmp[$i];
            my @split_r = split /\t/, $tmp[$i];
            $raw_count{$split_r[0]} = [$split_r[1],$split_r[2],0,0,0,0];
        }
    }
    close $tc;
    system("rm tmp_count");
    
    system("cat $dir/*_reads_counts_sens.txt > tmp_count_sens && rm $dir/*_reads_counts_sens.txt");
    #print STDERR ("INFO COMMANDE : "."bash -c 'cat $dir/*_reads_counts_sens.txt > tmp_count_sens && rm $dir/*_reads_counts_sens.txt'"."\n");
    
    open ( $tc, 'tmp_count_sens') or die "Cannot open tmp_count_sens: $!\n";
    @tmp = <$tc>;
    for (my $i = 0; $i <= $#tmp; $i++)
    {
        unless ($tmp[$i] =~ /^ID/)
        {
            chomp $tmp[$i];
            my @split_r = split /\t/, $tmp[$i];
            $raw_count{$split_r[0]}[2] =  $split_r[1];
            $raw_count{$split_r[0]}[3] =  $split_r[2];
        }
    }
    close $tc;
    system('rm tmp_count_sens');
    
    system("cat $dir/*_reads_counts_reverse.txt > tmp_count_reverse && rm $dir/*_reads_counts_reverse.txt");
    #print STDERR ("INFO COMMANDE : "."bash -c 'cat $dir/*_reads_counts_reverse.txt > tmp_count_reverse && rm $dir/*_reads_counts_reverse.txt'"."\n");
    
    open ( $tc, 'tmp_count_reverse') or die "Cannot open tmp_count_reverse: $!\n";
    @tmp = <$tc>;
    for (my $i = 0; $i <= $#tmp; $i++)
    {
        unless ($tmp[$i] =~ /^ID/)
        {
            chomp $tmp[$i];    
            my @split_r = split /\t/, $tmp[$i];
            $raw_count{$split_r[0]}[4] =  $split_r[1];
            $raw_count{$split_r[0]}[5] =  $split_r[2];
        }
    }
    close $tc;
    system("rm tmp_count_reverse");
    
    open (my $tcn, ">".$dir.'/count.txt') or die "Cannot open count: $!\n";    
    print $tcn "ID\treads count\trpkm\tsens reads count\trpkm\treverse reads count\trpkm\n";
    
    while (my ($k,$v) = each %raw_count)
    {
        print $tcn $k."\t".$v->[0]."\t".$v->[1]."\t".$v->[2]."\t".$v->[3]."\t".$v->[4]."\t".$v->[5]."\n";
    }
    close $tcn;
    
    # Sort summary count file
    system("cp $dir/count.txt $dir/countBckp.txt ; cat $dir/countBckp.txt | (sed -u 1q; sort) > $dir/count.txt ; rm $dir/countBckp.txt");
    #print STDERR ("INFO COMMANDE : "."bash -c 'cp $dir/count.txt $dir/countBckp.txt ; cat $dir/countBckp.txt | (sed -u 1q; sort) > $dir/count.txt ; rm $dir/countBckp.txt'"."\n");
    
    # Uncomment the following line to run manually
    #open (my $h,">".$dir."/".$html_out) or die "Cannot open $html_out: $!\n";
    # Comment the following line to run manually
    open (my $h,">"."/".$html_out) or die "Cannot open $html_out: $!\n";
    header($h);
    navbar($h);
    print $h "<div class=\"container\">  <div class=\"featurette\">  <p class=\"featurette-p\"> <A HREF=\"count.txt\"># of reads in each sequence</A></p></div></div>\n";
    print $h "<div class=\"container\"><p><a class=\"btn\" href=\"distri.html\">size distribution per sequence</a></p></div>\n";
    print $h "<div class=\"container\"><p><a class=\"btn\" href=\"align.html\">alignment on each sequence</a></p></div>\n";
    carousel($h,\@ref_names);
    footer($h);
    close $h;
    
    open (my $dis, ">".$dir."/distri.html") or die "Cannot open distri: $!n";
    header($dis);
    printDistri($dis,\@ref_names);
    footer($dis);
    close $dis;
    
    open (my $ali, ">".$dir."/align.html") or die "cannot open align: $!\n";
    header($ali);
    printAlign($ali,\@ref_names);
    footer($ali);
    close $ali;
    
    my $torm = $dir."/*_rejected_*";
    system("rm $torm");
    
    sub printAlign
    {
        my ($h, $tab) = @_;
        my ($fastq, $sam);
        my $cmp = 0;
        print $h "<div class=\"container\">\n";
        print $h "<div class=\"row text-center\">";
        foreach my $k (@{$tab})
        {
            $fastq = $name.$k.'.fastq'; 
            system("gzip '$dir'/'$fastq'");
            print STDERR "gzip '$dir'/'$fastq'";
            $fastq = $fastq.'.gz';
            $sam = $dir.'/'.$name.$k.'_sorted.sam';
            my $sam2 = $dir.'/'.$name.$k.'_sorted_mapped.bam';
            system("samtools view -@ $max_procs -Sbh -F4 $sam -o $sam2 && rm $sam");
            $sam2 = $name.$k.'_sorted_mapped.bam';
            print $h "</div><div class=\"row text-center\">" if $cmp != 0 && $cmp % 2 == 0;
            print $h "
            <div class=\"span6\">
            <h2>$k</h2>
            <p class=\"featurette-p\"><a href=\"$sam2\">bam file</a></p>
            <p class=\"featurette-p\"><a href=\"$fastq\">fastq.gz file</a></p>
            </div>
            ";
            $cmp++;
        }
        print $h "</div></div>";
    }
    
    sub printDistri
    {
        my ($h, $tab) = @_;
        my ($png, $txt);
        my $cmp = 0;
        print $h "<div class=\"container\">\n";
        print $h "<div class=\"row text-center\">";
        foreach my $k (@{$tab})
        {
            $txt = $name.$k.'_distribution.txt';
            $png = $name.$k.'_distribution.png';
            if ( -e $dir.'/'.$png )
            {
                print $h "</div><div class=\"row text-center\">" if $cmp != 0 && $cmp % 2 == 0;
                print $h "
                <div class=\"span6\">
                <h2>$k</h2>
                <p> <img src=\"$png\"/></p>
                <p class=\"featurette-p\"><a href=\"$txt\">text file</a></p>
                </div>
                ";
                $cmp++;
            }
        }
        print $h "</div></div>";
    }
    
    sub carousel
    {
        my ($file, $non_unique) = @_;
        my $ac = 0;
        print $file "
        <div id=\"page\">
        <div id=\"container\">
        <div class=\"each-gallery\">
        <div id=\"gallery\" class=\"content\">
        <div id=\"controls0\" class=\"controls\"></div>
        <div class=\"slideshow-container\">
        <div id=\"loading0\" class=\"loader\"></div>
        <div id=\"slideshow0\" class=\"slideshow\"></div>
        </div>
        <div id=\"caption0\" class=\"caption-container\"></div>
        </div>
        <div id=\"thumbs0\" class=\"navigation\">
        <ul class=\"thumbs noscript\">
        ";
        foreach my $u (@{$non_unique})
        {
            if (-e $dir.'/'.$u.'.jpg')
            {
                 print $file "
                 <li>
                 <a class=\"thumb\"  href=\"$u.jpg\" title=\"$u\">$u</a>
                 <div class=\"caption\">
                 <div class=\"download\">
                 <a href=\"$u.bed\">Bedfile</a>
                 </div>
                 </div>
                 </li>
                 ";
            }
        }
        print $file "
        </ul>
        </div>
        <div style=\"clear: both;\"></div></div>
        </div>
        </div>
        ";
    }
    
    sub navbar
    {
        my $file = shift;
        print $file "
        <div class=\"navbar navbar-inverse navbar-fixed-top\">
        <div class=\"navbar-inner\">
        <div class=\"container\">
        <button type=\"button\" class=\"btn btn-navbar\" data-toggle=\"collapse\" data-target=\".nav-collapse\">
        <span class=\"icon-bar\"></span>
        <span class=\"icon-bar\"></span>
        <span class=\"icon-bar\"></span>
        </button>
        <a class=\"brand\" href=\"report.txt\">Report</a>
        </div>
        </div>
        </div>";
    }
    
    sub header
    {
        my $file = shift;
        print $file "
        <!DOCTYPE html>
        <html lang=\"en\">
        <head>
        <meta charset=\"utf-8\">
        <title>pipeline</title>
        <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
        <meta name=\"description\" content=\"\">
        <meta name=\"author\" content=\"\">
        <!-- Le styles -->
        <link href=\"bootstrap.css\" rel=\"stylesheet\">
        <style type=\"text/css\">
        body {
          padding-top: 60px;
          padding-bottom: 40px;
        }
        div#page {
            width: 940px;
            background-color: #fff;
            margin: 0 auto;
            text-align: left;
            border-color: #fff;
            border-style: none solid solid;
            border-width: medium 1px 1px;
        }
        div.content {
        display: none;
        float: right;
        width: 550px;
        }
        div.content a, div.navigation a 
        {
            text-decoration: none;
            color: #777;
        }
        div.content a:focus, div.content a:hover, div.content a:active 
        {
            text-decoration: underline;
        }
        div.controls 
        {
            margin-top: 5px;
            height: 23px;
        }
        div.controls a 
        {
            padding: 5px;
        }
        div.ss-controls 
        {
            float: left;
        }
        div.nav-controls 
        {
            float: right;
        }
        div.slideshow-container 
        {
            position: relative;
            clear: both;
            height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
        }
        div.loader 
        {
            position: absolute;
            top: 0;
            left: 0;
            background-image: url('loader.gif');
            background-repeat: no-repeat;
            background-position: center;
            width: 550px;
            height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
        }
        div.slideshow {   }
        div.slideshow span.image-wrapper {
            display: block;
            position: absolute;
            top: 0;
            left: 0;
        }
        div.slideshow a.advance-link 
        {
            display: block;
            width: 550px;
            height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
            line-height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
            text-align: center;
        }
        div.slideshow a.advance-link:hover, div.slideshow a.advance-link:active, div.slideshow a.advance-link:visited 
        {
            text-decoration: none;
        }
        div.slideshow img 
        {
            vertical-align: middle;
            border: 1px solid #ccc;
        }
        div.image-title 
        {
            font-weight: bold;
            font-size: 1.4em;
        }
        div.image-desc 
        {
            line-height: 1.3em;
            padding-top: 12px;
        }
        div.navigation {    }
        ul.thumbs 
        {
            clear: both;
            margin: 0;
            padding: 0;
        }
        ul.thumbs li 
        {
            float: none;
            padding: 0;
            margin: 0;
            list-style: none;
        }
        a.thumb 
        {
            padding: 0;
            display: inline;
            border: none;
        }
        ul.thumbs li.selected a.thumb 
        {
            color: #000;
            font-weight: bold;
        }
        a.thumb:focus 
        {
            outline: none;
        }
        ul.thumbs img 
        {
            border: none;
            display: block;
        }
        div.pagination 
        {
            clear: both;
        }
        div.navigation div.top 
        {
            margin-bottom: 12px;
            height: 11px;
        }
        div.navigation div.bottom 
        {
            margin-top: 12px;
        }
        div.pagination a, div.pagination span.current, div.pagination span.ellipsis 
        {
            display: block;
            float: left;
            margin-right: 2px;
            padding: 4px 7px 2px 7px;
            border: 1px solid #ccc;
        }
        div.pagination a:hover 
        {
            background-color: #eee;
            text-decoration: none;
        }
        div.pagination span.current 
        {
            font-weight: bold;
            background-color: #000;
            border-color: #000;
            color: #fff;
        }
        div.pagination span.ellipsis 
        {
            border: none;
            padding: 5px 0 3px 2px;
        }
        div.download 
        {
            float: right;
        }
        div.caption-container 
        {
            position: relative;
            clear: left;
            height: 75px;
        }
        span.image-caption 
        {
            display: block;
            position: absolute;
            width: 550px;
            top: 0;
            left: 0;
        }
        div.caption 
        {
            padding: 12px;
        }
        /* Featurettes ------------------------- */
        .featurette 
        {
            padding-top: 20px; /* Vertically center images part 1: add padding above and below text. */
            overflow: hidden; /* Vertically center images part 2: clear their floats. */
            text-align: center;
        }
        .featurette-p
        {
            text-align: left;
        }
        .featurette-image 
        {
            margin-top: 10px; /* Vertically center images part 3: negative margin up the image the same amount of the padding to center it. */
            width: 600px;
            height: auto;
        }
        </style>
        <link href=\"bootstrap-responsive.css\" rel=\"stylesheet\">
        </head>
        <body>
        ";
    }
    
    sub footer 
    {
        my $file = shift;
        print $file "
        <!-- FOOTER -->
        <!-- Le javascript
        ================================================== -->
        <!-- Placed at the end of the document so the pages load faster -->
        <script src=\"jquery.js\"></script>
        <script src=\"bootstrap-transition.js\"></script>
        <script src=\"bootstrap-alert.js\"></script>
        <script src=\"bootstrap-modal.js\"></script>
        <script src=\"bootstrap-dropdown.js\"></script>
        <script src=\"bootstrap-scrollspy.js\"></script>
        <script src=\"bootstrap-tab.js\"></script>
        <script src=\"bootstrap-tooltip.js\"></script>
        <script src=\"bootstrap-popover.js\"></script>
        <script src=\"bootstrap-button.js\"></script>
        <script src=\"bootstrap-collapse.js\"></script>
        <script src=\"bootstrap-carousel.js\"></script>
        <script src=\"bootstrap-typeahead.js\"></script>
        <script src=\"holder.js\"></script>
        <script type=\"text/javascript\" src=\"jquery-1.3.2.js\"></script>
        <script type=\"text/javascript\" src=\"jquery.galleriffic.js\"></script>
        <script type=\"text/javascript\" src=\"jquery.opacityrollover.js\"></script>
        <script type=\"text/javascript\">
        jQuery(document).ready(function(\$) 
        {
            // We only want these styles applied when javascript is enabled
            \$('div.navigation').css({'width' : '300px', 'float' : 'left'});
            \$('div.content').css('display', 'block');
            \$(\".each-gallery\").each(function(i)
            {
                // Initially set opacity on thumbs and add
                // additional styling for hover effect on thumbs
                var onMouseOutOpacity = 0.67;
                \$('#thumbs + i + ul.thumbs li').opacityrollover(
                {
                    mouseOutOpacity:   onMouseOutOpacity,
                    mouseOverOpacity:  1.0,
                    fadeSpeed:         'fast',
                    exemptionSelector: '.selected'
                });
                // Initialize Advanced Galleriffic Gallery
                var gallery = \$('#thumbs'+i).galleriffic(
                    {
                        delay:                     2500,
                        numThumbs:                 22,
                        preloadAhead:              10,
                        enableTopPager:            true,
                        enableBottomPager:         true,
                        maxPagesToShow:            7,
                        imageContainerSel:         '#slideshow'+ i,
                        controlsContainerSel:      '#controls' + i,
                        captionContainerSel:       '#caption' + i,
                        loadingContainerSel:       '#loading' + i,
                        renderSSControls:          true,
                        renderNavControls:         true,
                        playLinkText:              'Play',
                        pauseLinkText:             'Pause',
                        prevLinkText:              '&lsaquo; Previous',
                        nextLinkText:              'Next &rsaquo;',
                        nextPageLinkText:          'Next &rsaquo;',
                        prevPageLinkText:          '&lsaquo; Prev',
                        enableHistory:             false,
                        autoStart:                 false,
                        syncTransitions:           true,
                        defaultTransitionDuration: 900,
                        onSlideChange:             function(prevIndex, nextIndex) 
                        {
                             // 'this' refers to the gallery, which is an extension of \$('#thumbs')
                            this.find('ul.thumbs').children()
                            .eq(prevIndex).fadeTo('fast', onMouseOutOpacity).end()
                            .eq(nextIndex).fadeTo('fast', 1.0);
                        },
                        onPageTransitionOut:       function(callback) 
                        {
                            this.fadeTo('fast', 0.0, callback);
                        },
                        onPageTransitionIn:        function() 
                        {
                            this.fadeTo('fast', 1.0);
                        }
                  });
            });
        });
        </script>
        </body>
        </html>
        ";
    }
} else 
{
    print "multi-allign version $VERSION
    Usage:
    multi-align.pl --fastq <fastq file> --ref <reference genome> --dir <results directory> --html <results.html> [options]
    Arguments:
    --dir\t\t\tFolder where results will be stored
    --fastq <fastq file>\tFastq file to process
    --html\t\t\tMain HTML file where results will be displayed
    --ref <reference>\t\tFasta file containing the reference genome
    Options:
    --ma <INT>\t\t\tNumber of reads (default: 0)
    --mis <INT>\t\t\tMaximal genome mismatches (default: 0)
    ";
}
