#!/usr/bin/env perl
# load modules
use warnings;
use Getopt::Std;
use File::Basename;
use Number::Format;
use Cwd;

my %opts;

# number format
my $de = new Number::Format(-thousands_sep =>',',-decimal_point => '.');

##########
## opts ##
##########
## input files
# b : path to input (b)am file
# t : path to input (t)arget regions in BED format
## output files
# o : report pdf (o)utput file
## entries in the report
# r : Coverage per (r)egion (boolean)
# s : (s)ubregion coverage if average < specified (plots for positions along target region) (boolean)
# S : (S)ubregion coverage for ALL failed exons  => use either s OR S or you will have double plots. 
# A : (A)ll exons will be plotted. 
# L : (L)ist failed exons instead of plotting
# m : (m)inimal Coverage threshold
# f : fraction of average as threshold
# n : sample (n)ame.
 

getopts('b:t:o:rsSALm:n:f:', \%opts) ;

my $tmp = getcwd();
# make output directory in (tmp) working dir
our $wd = "$tmp/Coverage.".int(rand(1000));
while (-d $wd) {
	$wd = "$tmp/Coverage.".int(rand(1000));
}
system("mkdir $wd");

## variables
our %commandsrun = ();
my ($thresh,$frac,$pdffile,$samplename,$totalmapped);


if (!exists($opts{'b'}) || !-e $opts{'b'}) {
	die('Bam File not found');
}
if (!exists($opts{'t'}) || !-e $opts{'t'}) {
	die('Target File (BED) not found');
}

if (exists($opts{'m'})) {
	$thresh = $opts{'m'};
}
else {
	$thresh = 40;
}

if (exists($opts{'f'})) {
	$frac = $opts{'f'};
}
else {
	$frac = 0.2;
}

if (exists($opts{'o'})) {
	$pdffile = $opts{'o'};
}
else {
	$pdffile = "$wd/CoverageReport.pdf";
}


# 1. Global Summary => default
&GlobalSummary($opts{'b'}, $opts{'t'});

# 2. Coverage per position
&SubRegionCoverage($opts{'b'}, $opts{'t'});
our %filehash;
if (exists($opts{'s'}) || exists($opts{'S'}) || exists($opts{'A'}) || exists($opts{'L'})) {
	system("mkdir $wd/SplitFiles");
	## get position coverages
	## split input files
	open IN, "$wd/Targets.Position.Coverage";
	my $fileidx = 0;
	my $currreg = '';
	while (<IN>) {
		my $line = $_;
		chomp($line);
		my @p = split(/\t/,$line);
		my $reg = $p[0].'-'.$p[1].'-'.$p[2]; #.$p[3];
		my $ex = $p[3];
		if ($reg ne $currreg) {
			## new exon open new outfile
			if ($currreg ne '') {
				## filehandle is open. close it
				close OUT; 
			}
			if (!exists($filehash{$reg})) {
				$fileidx++;
				$filehash{$reg}{'idx'} = $fileidx;
				$filehash{$reg}{'exon'} = $ex;
				open OUT, ">> $wd/SplitFiles/File_$fileidx.txt";
				$currreg = $reg;
			}
			else {
				open OUT, ">> $wd/SplitFiles/File_".$filehash{$reg}{'idx'}.".txt";
				$currreg = $reg;
			}
		}
		## print the line to the open filehandle. 	
		print OUT "$line\n";
	}
	close OUT;
	close IN;

}

## sort output files according to targets file
if (exists($opts{'r'}) ) {
	my %hash = ();
	open IN, "$wd/Targets.Global.Coverage";
	while (<IN>) {
            chomp;
            my @p = split(/\t/,$_) ;
            $hash{$p[3]} = $_;
	}
	close IN;
	open OUT, ">$wd/Targets.Global.Coverage";
	open IN, $opts{'t'};
	while (<IN>) {
            chomp;
            my @p = split(/\t/,$_) ;
            print OUT $hash{$p[3]} . "\n";
	}
	close IN;
	close OUT;
}


####################################
## PROCESS RESULTS & CREATE PLOTS ##
####################################
system("mkdir $wd/Report");

system("mkdir $wd/Rout");
system("mkdir $wd/Plots");

$samplename = $opts{'n'};
$samplename =~ s/_/\\_/g;
$samplename =~ s/\..+$//g;

# 0. Preamble
## compose preamble
open OUT, ">$wd/Report/Report.tex";
print OUT '\documentclass[a4paper,10pt]{article}'."\n";
print OUT '\usepackage[left=2cm,top=1.5cm,right=1.5cm,bottom=2.5cm,nohead]{geometry}'."\n";
print OUT '\usepackage{longtable}'."\n";
print OUT '\usepackage{fancyhdr}'."\n";
print OUT '\usepackage{fontspec}' . "\n";
print OUT '\usepackage{color}'."\n";
print OUT '\definecolor{grey}{RGB}{160,160,160}'."\n";
print OUT '\definecolor{darkgrey}{RGB}{100,100,100}'."\n";
print OUT '\definecolor{red}{RGB}{255,0,0}'."\n";
print OUT '\definecolor{orange}{RGB}{238,118,0}'."\n";
print OUT '\setlength\LTleft{0pt}'."\n";
print OUT '\setlength\LTright{0pt}'."\n";
print OUT '\begin{document}'."\n";
print OUT '\pagestyle{fancy}'."\n";
print OUT '\fancyhead{}'."\n";
print OUT '\renewcommand{\footrulewidth}{0.4pt}'."\n";
print OUT '\renewcommand{\headrulewidth}{0pt}'."\n";
print OUT '\fancyfoot[R]{\today\hspace{2cm}\thepage\ of \pageref{endofdoc}}'."\n";
print OUT '\fancyfoot[C]{}'."\n";
print OUT '\fancyfoot[L]{Coverage Report for ``'.$samplename.'"}'."\n";
print OUT '\let\oldsubsubsection=\subsubsection'."\n";
print OUT '\renewcommand{\subsubsection}{%'."\n";
print OUT '  \filbreak'."\n";
print OUT '  \oldsubsubsection'."\n";
print OUT '}'."\n";
# main title
print OUT '\section*{Coverage Report for ``'.$samplename.'"}'."\n";
close OUT; 

# 1. Summary Report
# Get samtools flagstat summary of BAM file	
my $flagstat = `samtools flagstat $opts{'b'}`;
my @s = split(/\n/,$flagstat);
# Get number of reads mapped in total 
## updated on 2012-10-1 !!
$totalmapped = $s[2];
$totalmapped =~ s/^(\d+)(\s.+)/$1/;

if ( not $totalmapped){
    exit;
}

# count columns
my $head = `head -n 1 $wd/Targets.Global.Coverage`;
chomp($head);
my @cols = split(/\t/,$head);
my $nrcols = scalar(@cols);
my $covcol = $nrcols - 3;
# get min/max/median/average coverage => values
my $covs = `cut -f $covcol $wd/Targets.Global.Coverage`;
my @coverages = split(/\n/,$covs);
my ($eavg,$med,$min,$max,$first,$third,$ontarget) = arraystats(@coverages);
my $spec = 0;
$spec=sprintf("%.1f",($ontarget / $totalmapped)*100);
# get min/max/median/average coverage => boxplot in R
open OUT, ">$wd/Rout/boxplot.R";
print OUT 'coverage <- read.table("../Targets.Global.Coverage",as.is=TRUE,sep="\t",header=FALSE)'."\n";
print OUT 'coverage <- coverage[,'.$covcol.']'."\n";
print OUT 'png(file="../Plots/CoverageBoxPlot.png", bg="white", width=240, height=480)'."\n";
print OUT 'boxplot(coverage,range=1.5,main="Target Region Coverage")'."\n";
print OUT 'graphics.off()'."\n";
close OUT;
system("cd $wd/Rout && Rscript boxplot.R") == 0 || die "Could not run boxplot.R with following error '$!'\n";

## global nt coverage plot
## use perl to make histogram (lower memory)
open IN, "$wd/Targets.Position.Coverage";
my %dens;
my $counter = 0;
my $sum = 0;
my $avg = 0;

while (<IN>) {
	chomp();
	my @p = split(/\t/);
	$sum += $p[-1];
	$counter++;
	if (defined($dens{$p[-1]})) {
		$dens{$p[-1]}++;
	}
	else {
		$dens{$p[-1]} = 1;
	}
}
$avg = $sum/$counter;
close IN;
open OUT, ">$wd/Rout/hist.txt";
if (!defined($dens{'0'})) {
	$dens{'0'} = 0;
}
foreach (keys(%dens)) {
	print OUT "$_;$dens{$_}\n";
}
close OUT;
open OUT, ">$wd/Rout/ntplot.R";
# read coverage hist in R to plot
print OUT 'coverage <- read.table("hist.txt" , as.is = TRUE, header=FALSE,sep=";")'."\n";
print OUT 'mincov <- '."$thresh \n";
print OUT "avg <- round($avg,digits=8)\n";
print OUT "colnames(coverage) <- c('cov','count')\n";
print OUT 'coverage$cov <- coverage$cov / avg'."\n";
print OUT 'rep <- which(coverage$cov > 1)'."\n";
print OUT 'coverage[coverage$cov > 1,1] <- 1'."\n";
print OUT 'values <- coverage[coverage$cov < 1,]'."\n";
print OUT 'values <- rbind(values,c(1,sum(coverage[coverage$cov == 1,"count"])))'."\n";
print OUT 'values <- values[order(values$cov),]'."\n";
print OUT 'prevcount <- 0'."\n";
# make cumulative count data frame
print OUT 'for (i in rev(values$cov)) {'."\n";
print OUT '  values[values$cov == i,"count"] <- prevcount + values[values$cov == i,"count"]'."\n";
print OUT '  prevcount <- values[values$cov == i,"count"]'."\n";
print OUT '}'."\n";
print OUT 'values$count <- values$count / (values[values$cov == 0,"count"] / 100)'."\n";
# get some values to plot lines.
print OUT 'mincov.x <- mincov/avg'."\n";
print OUT 'if (mincov/avg <= 1) {'."\n";
print OUT '  ii <- which(values$cov == mincov.x)'."\n";
print OUT '  if (length(ii) == 1) {'."\n";
print OUT '    mincov.y <- values[ii[1],"count"]'."\n";
print OUT '  } else {'."\n";
print OUT '    i1 <- max(which(values$cov < mincov.x))'."\n";
print OUT '    i2 <- min(which(values$cov > mincov.x))'."\n";
print OUT '    mincov.y <- ((values[i2,"count"] - values[i1,"count"])/(values[i2,"cov"] - values[i1,"cov"]))*(mincov.x - values[i1,"cov"]) + values[i1,"count"]'."\n";
print OUT '  }'."\n";
print OUT '}'."\n";
# open output image and create plot
print OUT 'png(file="../Plots/CoverageNtPlot.png", bg="white", width=540, height=480)'."\n";
print OUT 'par(xaxs="i",yaxs="i")'."\n";
print OUT 'plot(values$cov,values$count,ylim=c(0,100),pch=".",main="Cumulative Normalised Base-Coverage Plot",xlab="Normalizalised Coverage",ylab="Cumulative Nr. Of Bases")'."\n";
print OUT 'lines(values$cov,values$count)'."\n";
print OUT 'if (mincov.x <= 1) {'."\n";
print OUT '  lines(c(mincov.x,mincov.x),c(0,mincov.y),lty=2,col="darkgreen")'."\n";
print OUT '  lines(c(0,mincov.x),c(mincov.y,mincov.y),lty=2,col="darkgreen")'."\n";
print OUT '  text(1,(95),pos=2,col="darkgreen",labels="Threshold: '.$thresh.'x")'."\n";
print OUT '  text(1,(91),pos=2,col="darkgreen",labels=paste("%Bases: ",round(mincov.y,2),"%",sep=""))'."\n";
print OUT '} else {'."\n";
print OUT '  text(1,(95),pos=2,col="darkgreen",labels="Threshold ('.$thresh.'x) > Average")'."\n";
print OUT '  text(1,(91),pos=2,col="darkgreen",labels="Plotting impossible")'."\n";
print OUT '}'."\n";
print OUT 'frac.x <- '."$frac\n";
print OUT 'ii <- which(values$cov == frac.x)'."\n";
print OUT 'if (length(ii) == 1) {'."\n";
print OUT '  frac.y <- values[ii[1],"count"]'."\n";
print OUT '} else {'."\n";
print OUT '  i1 <- max(which(values$cov < frac.x))'."\n";
print OUT '  i2 <- min(which(values$cov > frac.x))'."\n";
print OUT '  frac.y <- ((values[i2,"count"] - values[i1,"count"])/(values[i2,"cov"] - values[i1,"cov"]))*(frac.x - values[i1,"cov"]) + values[i1,"count"]'."\n";
print OUT '}'."\n";
print OUT 'lines(c(frac.x,frac.x),c(0,frac.y),lty=2,col="red")'."\n";
print OUT 'lines(c(0,frac.x),c(frac.y,frac.y),lty=2,col="red")'."\n";
#iprint OUT 'text((frac.x+0.05),(frac.y - 2),pos=4,col="red",labels=paste(frac.x," x Avg.Cov : ",round(frac.x * avg,2),"x",sep="" ))'."\n";
#print OUT 'text((frac.x+0.05),(frac.y-5),pos=4,col="red",labels=paste("%Bases: ",round(frac.y,2),"%",sep=""))'."\n";
print OUT 'text(1,86,pos=2,col="red",labels=paste(frac.x," x Avg.Cov : ",round(frac.x * avg,2),"x",sep="" ))'."\n";
print OUT 'text(1,82,pos=2,col="red",labels=paste("%Bases: ",round(frac.y,2),"%",sep=""))'."\n";

print OUT 'graphics.off()'."\n";

close OUT;
system("cd $wd/Rout && Rscript ntplot.R") == 0 || die "Could not run ntplot.R with following error '$!'\n";
## PRINT TO .TEX FILE
open OUT, ">>$wd/Report/Report.tex";
# average coverage overviews
print OUT '\subsection*{Overall Summary}'."\n";
print OUT '{\small ';
# left : boxplot
print OUT '\begin{minipage}{0.3\linewidth}\centering'."\n";
print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/CoverageBoxPlot.png}'."\n";
print OUT '\end{minipage}'."\n";
# right : cum.cov.plot
print OUT '\hspace{0.6cm}'."\n";
print OUT '\begin{minipage}{0.65\linewidth}\centering'."\n";
print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/CoverageNtPlot.png}'."\n";
print OUT '\end{minipage} \\\\'."\n";
## next line
print OUT '\begin{minipage}{0.48\linewidth}'."\n";
print OUT '\vspace{-1.2em}'."\n";
print OUT '\begin{tabular}{ll}'."\n";
# bam statistics
print OUT '\multicolumn{2}{l}{\textbf{\underline{Samtools Flagstat Summary}}} \\\\'."\n";
foreach (@s) {
	$_ =~ m/^(\d+)\s(.+)$/;
	my $one = $1;
	my $two = $2;
	$two =~ s/\s\+\s0\s//;
	$two = ucfirst($two);
	$one =~ s/%/\\%/g;
	# remove '+ 0 ' from front
	$two =~ s/\+\s0\s//;
	# remove trailing from end
	$two =~ s/(\s\+.*)|(:.*)/\)/;
	$two =~ s/%/\\%/g;
	$two =~ s/>=/\$\\ge\$/g;
	$two = ucfirst($two);
	print OUT '\textbf{'.$two.'} & '.$one.' \\\\'."\n";
}
print OUT '\end{tabular}\end{minipage}'."\n";
print OUT '\hspace{1.5cm}'."\n";
# target coverage statistics
print OUT '\begin{minipage}{0.4\linewidth}'."\n";
#print OUT '\vspace{-4.8em}'."\n";
print OUT '\begin{tabular}{ll}'."\n";
print OUT '\multicolumn{2}{l}{\textbf{\underline{Target Region Coverage}}} \\\\'."\n";
print OUT '\textbf{Number of Target Regions} & '.scalar(@coverages).' \\\\'."\n";
print OUT '\textbf{Minimal Region Coverage} & '.$min.' \\\\'."\n";
print OUT '\textbf{25\% Region Coverage} & '.$first.' \\\\'. "\n";
print OUT '\textbf{50\% (Median) Region Coverage} & '.$med.' \\\\'. "\n";
print OUT '\textbf{75\% Region Coverage} & '.$third.' \\\\'. "\n";
print OUT '\textbf{Maximal Region Coverage} & '.$max.' \\\\'. "\n";
print OUT '\textbf{Average Region Coverage} & '.int($eavg).' \\\\'. "\n";
print OUT '\textbf{Mapped On Target} & '.$spec.' \\\\'."\n";
print OUT '\multicolumn{2}{l}{\textbf{\underline{Target Base Coverage }}} \\\\'."\n";
print OUT '\textbf{Number of Target Bases} & '.$counter.' \\\\'."\n";
print OUT '\textbf{Average Base Coverage} & '.int($avg).' \\\\'. "\n";
print OUT '\textbf{Non-Covered Bases} & '.$dens{'0'}.' \\\\'."\n";
#print OUT '\textbf{Bases Covered $ge$ '.$frac.'xAvg.Cov} & '.
print OUT '\end{tabular}\end{minipage}}'."\n";
close OUT;

# 2. GLOBAL COVERAGE OVERVIEW PER GENE
my @failedexons;
my @allexons;
my @allregions;
my @failedregions;
if (exists($opts{'r'}) || exists($opts{'s'}) || exists($opts{'S'})) {
	# count columns
	my $head = `head -n 1 $wd/Targets.Global.Coverage`;
	chomp($head);
	my @cols = split(/\t/,$head);
	my $nrcols = scalar(@cols);
	my $covcol = $nrcols - 3;
	# Coverage Plots for each gene => barplots in R, table here.
	open IN, "$wd/Targets.Global.Coverage";
	my $currgroup = '';
	my $startline = 0;
	my $stopline = 0;
	my $linecounter = 0;
	my $nrcol=0;
	while (<IN>) {
		$linecounter++;
		chomp($_);
		my @c = split(/\t/,$_);
		push(@allregions,$c[0].'-'.$c[1].'-'.$c[2]);
		my $group = $c[3];
		## coverage failure?
		if ($c[$nrcol-1] < 1 || $c[$covcol-1] < $thresh) {
			push(@failedexons,$group);
			push(@failedregions,$c[0].'-'.$c[1].'-'.$c[2]);
		}
		## store exon
		push(@allexons,$group);
		## extract and check gene
		$group =~ s/^(\S+)[\|\s](.+)/$1/;

                my $scale;
		if ($group ne $currgroup ) {
		    if ($currgroup ne '') {
			# new gene, make plot. 
			open OUT, ">$wd/Rout/barplot.R";
			print OUT 'coveragetable <- read.table("../Targets.Global.Coverage",as.is=TRUE,sep="\t",header=FALSE)'."\n";
			print OUT 'coverage <- coveragetable[c('.$startline.':'.$stopline.'),'.$covcol.']'."\n";
			print OUT 'entries <- coveragetable[c('.$startline.':'.$stopline.'),4]'."\n";
			print OUT 'entries <- sub("\\\\S+\\\\|","",entries,perl=TRUE)'."\n";
			print OUT 'coverage[coverage < 1] <- 1'."\n";
			print OUT 'colors <- c(rep("grey",length(coverage)))'."\n";
			# coverage not whole target region => orange
			print OUT 'covperc <- coveragetable[c('.$startline.':'.$stopline.'),'.$nrcols.']'."\n";
			print OUT 'colors[covperc<1] <- "orange"'."\n";
			# coverage below threshold => red
			print OUT 'colors[coverage<'.$thresh.'] <- "red"'."\n";

			if ($stopline - $startline > 20) {
				$scale = 2;
			}
			else {
				$scale = 1;
			}
			my $width = 480 * $scale;
			my $height = 240 * $scale;
			print OUT 'png(file="../Plots/Coverage_'.$currgroup.'.png", bg="white", width='.$width.', height='.$height.')'."\n";
			print OUT 'ylim = c(0,max(max(log10(coverage),log10('.($thresh+20).'))))'."\n";
			print OUT 'mp <- barplot(log10(coverage),col=colors,main="Exon Coverage for '.$currgroup.'",ylab="Log10(Coverage)",ylim=ylim)'."\n";
			print OUT 'text(mp, log10(coverage) + '.(0.4/$scale).',format(coverage),xpd = TRUE,srt=90)'."\n";
			print OUT 'text(mp,par("usr")[3]-0.05,labels=entries,srt=45,adj=1,xpd=TRUE)'."\n";
			print OUT 'abline(h=log10('.$thresh.'),lwd=4,col=rgb(255,0,0,100,maxColorValue=255))'."\n";
			print OUT 'graphics.off()'."\n";
			close OUT;
			system("cd $wd/Rout && Rscript barplot.R") == 0 || die "Could not run barplot.R with following error '$!'\n";
			if ($scale == 1) {
				push(@small,'\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$currgroup.'.png}');
			}
			else {
				push(@large,'\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$currgroup.'.png}');
			}

		    }
		    $currgroup = $group;
		    $startline = $linecounter;
		}
		$stopline = $linecounter;
	}
	close IN;
	if ($currgroup ne '') {
		# last gene, make plot. 
		open OUT, ">$wd/Rout/barplot.R";
		print OUT 'coveragetable <- read.table("../Targets.Global.Coverage",as.is=TRUE,sep="\t",header=FALSE)'."\n";
		print OUT 'coverage <- coveragetable[c('.$startline.':'.$stopline.'),'.$covcol.']'."\n";
		print OUT 'entries <- coveragetable[c('.$startline.':'.$stopline.'),4]'."\n";
		print OUT 'entries <- sub("\\\\S+\\\\|","",entries,perl=TRUE)'."\n";
		print OUT 'coverage[coverage < 1] <- 1'."\n";
		print OUT 'colors <- c(rep("grey",length(coverage)))'."\n";
		print OUT 'colors[coverage<'.$thresh.'] <- "red"'."\n";

		if ($stopline - $startline > 20) {
			$scale = 2;
		}
		else {
			$scale = 1;
		}
		my $width = 480 * $scale;
		my $height = 240 * $scale;
		print OUT 'png(file="../Plots/Coverage_'.$currgroup.'.png", bg="white", width='.$width.', height='.$height.')'."\n";
		print OUT 'ylim = c(0,max(max(log10(coverage),log10('.($thresh+20).'))))'."\n";
		print OUT 'mp <- barplot(log10(coverage),col=colors,main="Exon Coverage for '.$currgroup.'",ylab="Log10(Coverage)", ylim=ylim)'."\n";
		print OUT 'text(mp, log10(coverage) + log10(2),format(coverage),xpd = TRUE,srt=90)'."\n";
		print OUT 'text(mp,par("usr")[3]-0.1,labels=entries,srt=45,adj=1,xpd=TRUE)'."\n";
		print OUT 'abline(h=log10('.$thresh.'),lwd=4,col=rgb(255,0,0,100,maxColorValue=255))'."\n";
		print OUT 'graphics.off()'."\n";
		close OUT;
		system("cd $wd/Rout && Rscript barplot.R") == 0 || die "Could not run barplot.R with following error '$!'\n";
		if ($scale == 1) {
			push(@small,'\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$currgroup.'.png}');
		}
		else {
			push(@large,'\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$currgroup.'.png}');
		}
	}
	## print to TEX
	open OUT, ">>$wd/Report/Report.tex";
	print OUT '\subsection*{Gene Summaries}'."\n";
	print OUT '\underline{Legend:} \\\\'."\n";
	print OUT '{\color{red}\textbf{RED:} Coverage did not reach set threshold of '.$thresh.'} \\\\'."\n";
	print OUT '{\color{orange}\textbf{ORANGE:} Coverage was incomplete for the exon. Overruled by red.} \\\\' ."\n";
	my $col = 1;
	foreach (@small) {
		if ($col > 2) {
			$col = 1;
			print OUT "\n";
		}
		print OUT '\begin{minipage}{0.5\linewidth}\centering'."\n";
		print OUT $_."\n";
		print OUT '\end{minipage}'."\n";
		$col++;
	}
	## new line 
	if ($col == 2) {
		print OUT '\\\\'." \n";
	}
	foreach(@large) {
		print OUT $_."\n";
	}
	close OUT;
	
}

# 3. Detailed overview of failed exons (globally failed)
if (exists($opts{'s'})) {
	# count columns
	my $head = `head -n 1 $wd/Targets.Position.Coverage`;
	chomp($head);
        my $subtitle;
	my @cols = split(/\t/,$head);
	my $nrcols = scalar(@cols);
	my $covcol = $nrcols;
	my $poscol = $nrcols -1;
	# tex section header
	open TEX, ">>$wd/Report/Report.tex";
	print TEX '\subsection*{Failed Exon Plots}'."\n";
	$col = 1;
	print TEX '\underline{NOTE:} Only exons with global coverage $<$'.$thresh.' or incomplete coverage were plotted \\\\'."\n";
	foreach(@failedregions) {
		if ($col > 2) {
			$col = 1;
			print TEX "\n";
		}
		# which exon
		my $region = $_;
		my $exon = $filehash{$region}{'exon'};
		# link exon to tmp file
		my $exonfile = "$wd/SplitFiles/File_".$filehash{$region}{'idx'}.".txt";
		## determine transcript orientation and location
		my $firstline = `head -n 1 $exonfile`;
		my @firstcols = split(/\t/,$firstline);
		my $orient = $firstcols[5];
		my $genomicchr = $firstcols[0];
		my $genomicstart = $firstcols[1];
		my $genomicstop = $firstcols[2];
		if ($orient eq '+') {
			$bps = $genomicstop - $genomicstart + 1;
			$subtitle = "Region 0-$bps: $genomicchr:".$de->format_number($genomicstart)."+".$de->format_number($genomicstop);
		}
		else {
			$bps = $genomicstop - $genomicstart + 1;
			$subtitle = "Region 0-$bps: $genomicchr:".$de->format_number($genomicstart)."-".$de->format_number($genomicstop);
		}
		# print Rscript
		open OUT, ">$wd/Rout/exonplot.R";
		print OUT 'coveragetable <- read.table("'.$exonfile.'",as.is=TRUE,sep="\t",header=FALSE)'."\n";
		print OUT 'coverage <- coveragetable[,'.$covcol.']'."\n";
		print OUT 'coverage[coverage < 1] <- 1'."\n";
		print OUT 'positions <- coveragetable[,'.$poscol.']'."\n";
		
		my $width = 480 ;
		my $height = 240 ;
		my $exonstr = $exon;
		$exonstr =~ s/\s/_/g;
		$exon =~ s/_/ /g;
		$exon =~ s/\|/ /g;
		print OUT 'png(file="../Plots/Coverage_'.$exonstr.'.png", bg="white", width='.$width.', height='.$height.')'."\n";
		print OUT 'ylim = c(0,log10(max(max(coverage),'.($thresh+10).')))'."\n";
		if ($orient eq '-') {
			print OUT 'plot(positions,log10(coverage),type="n",main="Coverage for '.$exon.'",ylab="log10(Coverage)",ylim=ylim,xlab="Position",xlim=rev(range(positions)),sub="(Transcribed from minus strand)")'."\n";
			print OUT 'mtext("'.$subtitle.'")'."\n";
		}
		else {
			print OUT 'plot(positions,log10(coverage),type="n",main="Coverage for '.$exon.'",ylab="log10(Coverage)",ylim=ylim,xlab="Position",sub="(Transcribed from plus strand)")'."\n";
			print OUT 'mtext("'.$subtitle.'")'."\n";
		}
		print OUT 'lines(positions,log10(coverage))'."\n";
		print OUT 'abline(h=log10('.$thresh.'),lwd=4,col=rgb(255,0,0,100,maxColorValue=255))'."\n";
		print OUT 'failedpos <- positions[coverage<'.$thresh.']'."\n";
		print OUT 'failedcov <- coverage[coverage<'.$thresh.']'."\n";
		print OUT 'points(failedpos,log10(failedcov),col="red",pch=19)'."\n";
		print OUT 'graphics.off()'."\n";
		close OUT;
		# run R script
		system("cd $wd/Rout && Rscript exonplot.R")== 0 || die "Could not run exonplot.R with following error '$!'\n";
		# Add to .TEX
		print TEX '\begin{minipage}{0.5\linewidth}\centering'."\n";
		print TEX '\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$exonstr.'.png}'."\n";
		print TEX '\end{minipage}'."\n";
		$col++;
	}
}

## plot failed (subregion) or all exons
if (exists($opts{'S'}) || exists($opts{'A'})) {
	# count columns
	my $head = `head -n 1 $wd/Targets.Position.Coverage`;
	chomp($head);
	my @cols = split(/\t/,$head);
	my $nrcols = scalar(@cols);
	my $covcol = $nrcols;
	my $poscol = $nrcols -1;
	# tex section header
	open TEX, ">>$wd/Report/Report.tex";
	print TEX '\subsection*{Failed Exon Plots}'."\n";
	if (exists($opts{'S'})) {
		print TEX '\underline{NOTE:} ALL exons were tested for local coverage $<$'.$thresh.' \\\\'."\n";
	}
	elsif (exists($opts{'A'})) {
		print TEX '\underline{NOTE:} ALL exons are plotted, regardless of coverage \\\\'."\n";
	}
	$col = 1;
	foreach(@allregions) {
		if ($col > 2) {
			$col = 1;
			print TEX "\n";
		}
		# which exon
		my $region = $_;
		my $exon = $filehash{$region}{'exon'};
		# grep exon to tmp file
		my $exonfile = "$wd/SplitFiles/File_".$filehash{$region}{'idx'}.".txt";
		## determine transcript orientation.
                my $firstline = `head -n 1 $exonfile`;
                my @firstcols = split(/\t/,$firstline);
                my $orient = $firstcols[5];
		my $genomicchr = $firstcols[0];
		my $genomicstart = $firstcols[1];
		my $genomicstop = $firstcols[2];
                my $subtitle;
                
		if ($orient eq '+') {
			my $bps = $genomicstop - $genomicstart + 1;
			$subtitle = "Region 0-$bps: $genomicchr:".$de->format_number($genomicstart)."+".$de->format_number($genomicstop);

		}
		else {
			my $bps = $genomicstop - $genomicstart + 1;
			$subtitle = "Region 0-$bps: $genomicchr:".$de->format_number($genomicstart)."-".$de->format_number($genomicstop);

		}

		# check if failed
		if (exists($opts{'S'})) {
			my $cs = `cut -f $covcol '$exonfile' `;
			my @c = split(/\n/,$cs);
			@c = sort { $a <=> $b } @c;
			if ($c[0] >= $thresh) {
				# lowest coverage > threshold => skip
				next;
			}
		}
		# print Rscript
		open OUT, ">$wd/Rout/exonplot.R";
		print OUT 'coveragetable <- read.table("'.$exonfile.'",as.is=TRUE,sep="\t",header=FALSE)'."\n";
		print OUT 'coverage <- coveragetable[,'.$covcol.']'."\n";
		print OUT 'coverage[coverage < 1] <- 1'."\n";
		print OUT 'positions <- coveragetable[,'.$poscol.']'."\n";
		my $width = 480 ;
		my $height = 240 ;
		my $exonstr = $exon;
		$exonstr =~ s/\s/_/g;
		$exon =~ s/_/ /g;
		$exon =~ s/\|/ /g;
		print OUT 'png(file="../Plots/Coverage_'.$exonstr.'.png", bg="white", width='.$width.', height='.$height.')'."\n";
		print OUT 'ylim = c(0,log10(max(max(coverage),'.($thresh+10).')))'."\n";
		if ($orient eq '-') {
			print OUT 'plot(positions,log10(coverage),type="n",main="Coverage for '.$exon.'",ylab="log10(Coverage)",ylim=ylim,xlab="Position",xlim=rev(range(positions)),sub="(Transcribed from minus strand)")'."\n";
			print OUT 'mtext("'.$subtitle.'")'."\n";
                }
                else {
			print OUT 'plot(positions,log10(coverage),type="n",main="Coverage for '.$exon.'",ylab="log10(Coverage)",ylim=ylim,xlab="Position",sub="(Transcribed from plus strand)")'."\n";
			print OUT 'mtext("'.$subtitle.'")'."\n";
                }
		
		print OUT 'lines(positions,log10(coverage))'."\n";
		print OUT 'abline(h=log10('.$thresh.'),lwd=4,col=rgb(255,0,0,100,maxColorValue=255))'."\n";
		print OUT 'failedpos <- positions[coverage<'.$thresh.']'."\n";
		print OUT 'failedcov <- coverage[coverage<'.$thresh.']'."\n";
		print OUT 'points(failedpos,log10(failedcov),col="red",pch=19)'."\n";
		print OUT 'graphics.off()'."\n";
		close OUT;
		# run R script
		system("cd $wd/Rout && Rscript exonplot.R") == 0 || die "Could not run exonplot.R with following error '$!'\n";
		# Add to .TEX
		print TEX '\begin{minipage}{0.5\linewidth}\centering'."\n";
		print TEX '\includegraphics[width=\textwidth,keepaspectratio=true]{../Plots/Coverage_'.$exonstr.'.png}'."\n";
		print TEX '\end{minipage}'."\n";
		$col++;
	}
}
## list failed exons
if (exists($opts{'L'})) {
    # count columns
        my $subtitle;
	my $head = `head -n 1 $wd/Targets.Position.Coverage`;
	chomp($head);
	my @cols = split(/\t/,$head);
	my $nrcols = scalar(@cols);
	my $covcol = $nrcols;
	my $poscol = $nrcols -1;
	## hash to print
	# tex section header
	open TEX, ">>$wd/Report/Report.tex";
	print TEX '\subsection*{List of Failed Exons}'."\n";
	print TEX '\underline{NOTE:} ALL exons were tested for local coverage $<$'.$thresh.' \\\\'."\n";
	print TEX '{\footnotesize\begin{longtable}[l]{@{\extracolsep{\fill}}llll}'."\n".'\hline'."\n";
	print TEX '\textbf{Target Name} & \textbf{Genomic Position} & \textbf{Avg.Coverage} & \textbf{Min.Coverage} \\\\'."\n".'\hline'."\n";
	print TEX '\endhead'."\n";
	print TEX '\hline '."\n".'\multicolumn{4}{r}{{\textsl{\footnotesize Continued on next page}}} \\\\ '."\n".'\hline' ."\n". '\endfoot' . "\n". '\endlastfoot' . "\n";
 
	$col = 1;
	open IN, "$wd/Targets.Global.Coverage";
	while (<IN>) {
		chomp($_);
		my @p = split(/\t/,$_);
		my $region = $p[0].'-'.$p[1].'-'.$p[2];
		my $exon = $filehash{$region}{'exon'};
		# grep exon to tmp file
		my $exonfile = "$wd/SplitFiles/File_".$filehash{$region}{'idx'}.".txt";
		## determine transcript orientation.
                my $firstline = `head -n 1 $exonfile`;
                my @firstcols = split(/\t/,$firstline);
                my $orient = $firstcols[5];
		my $genomicchr = $firstcols[0];
		my $genomicstart = $firstcols[1];
		my $genomicstop = $firstcols[2];

		if ($orient eq '+') {
			my $bps = $genomicstop - $genomicstart + 1;
			$subtitle = "$genomicchr:".$de->format_number($genomicstart)."+".$de->format_number($genomicstop);

		}
		else {
			my $bps = $genomicstop - $genomicstart + 1;
			$subtitle = "$genomicchr:".$de->format_number($genomicstart)."-".$de->format_number($genomicstop);
		}

		# check if failed
		my $cs = `cut -f $covcol '$exonfile' `;
		my @c = split(/\n/,$cs);
		my ($avg,$med,$min,$max,$first,$third,$ontarget) = arraystats(@c);

		if ($min >= $thresh) {
			# lowest coverage > threshold => skip
			next;
		}

		# print to .tex table 
		if (length($exon) > 30) {
			$exon = substr($exon,0,27) . '...';
		}
		$exon =~ s/_/ /g;
		$exon =~ s/\|/ /g;

		print TEX "$exon & $subtitle & ".int($avg)." & $min ".'\\\\'."\n"; 
	}
	close IN;
	print TEX '\hline'."\n";
	print TEX '\end{longtable}}'."\n";
	close TEX;
}


## Close document
open OUT, ">>$wd/Report/Report.tex";
print OUT '\label{endofdoc}'."\n";
print OUT '\end{document}'."\n";
close OUT;
system("cd $wd/Report && tectonic Report.tex") == 0 || die "Could not run tectonic with following error '$!'\n";

## mv report to output file
system("cp -f $wd/Report/Report.pdf '$pdffile'");
##create tar.gz file
system("mkdir $wd/Results");
system("cp -Rf $wd/Plots $wd/Results/");
system("cp -Rf $wd/Report/ $wd/Results/");
if (-e "$wd/Targets.Global.Coverage") {
	system("cp -Rf $wd/Targets.Global.Coverage $wd/Results/");
}
if (-e "$wd/Targets.Position.Coverage") {
	system("cp -Rf $wd/Targets.Position.Coverage $wd/Results/");
}


exit;

###############
## FUNCTIONS ##
###############
sub arraystats{
	my @array = @_;
	my $count = scalar(@array);
	@array = sort { $a <=> $b } @array;
	# median
	my $median = 0;
	if ($count % 2) { 
		$median = $array[int($count/2)]; 
	} else { 
		$median = ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
	# average 
	my $sum = 0;
	foreach (@array) { $sum += $_; } 
	my $average = $sum / $count;
	# quantiles (rounded)
	my $quart = int($count/4) ;	
	my $first = $array[$quart];
	my $third = $array[($quart*3)];
	my $min = $array[0];
	my $max = $array[($count-1)];
	return ($average,$median,$min,$max,$first,$third,$sum);
}

sub GlobalSummary {
	my ($bam,$targets) = @_;

	my $command = "cd $wd && coverageBed -abam $bam -b $targets > $wd/Targets.Global.Coverage";
	if (exists($commandsrun{$command})) {
		return;
	}
	system($command);
	$commandsrun{$command} = 1;
}

sub CoveragePerRegion {
	my ($bam,$targets) = @_;
	my $command = "cd $wd && coverageBed -abam $bam -b $targets > $wd/Targets.Global.Coverage";
	if (exists($commandsrun{$command})) {
		return;
	}
	system($command);
	$commandsrun{$command} = 1;
}

sub SubRegionCoverage {
	my ($bam,$targets) = @_;
	my $command = "cd $wd && coverageBed -abam $bam -b $targets -d > $wd/Targets.Position.Coverage";
	system($command);
	$commandsrun{$command} = 1;
}

