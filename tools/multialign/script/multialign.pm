#!/usr/bin/env perl

package multialign;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter); 
our @EXPORT_OK = qw(&align &histogram &histogram_Gviz &rpms_rpkm &size &windows); 
use Carp;
use File::Basename;
use Statistics::R;
use POSIX;
use String::Random;
use Parallel::ForkManager;

sub align
{
    my ($index, $mismatches, $fastq, $fastq_accepted, $fastq_accepted_unique, $fastq_rejected, $sorted_sam_accepted, $duplicate_hashR, $best_hit_number_hashR, $number_of_cpus) = @_ ;
    my ($reads, $mappers, $mappersUnique, @garbage, %size_num, %size_num_spe, %number, %numberSens, %numberReverse, %unique_number, %numberNM, %numberM, %size);
    $mappers = $mappersUnique = $reads = 0;
    my $sam = new String::Random;
    $sam = $sam->randpattern("CCcccccc");
    $number_of_cpus = Sys::CPU::cpu_count() unless defined($number_of_cpus);
    open (O,$fastq) or die "Cannot open $fastq: $!\n";
    my @fastq = <O>; $reads = scalar(@fastq) / 4;
    my $bam_sorted = $sam."sorted";
    my $bam_sorted_bam = $sam."sorted.bam"; push (@garbage, $bam_sorted_bam);
    `bwa aln -t $number_of_cpus -n $mismatches '$index' $fastq 2> /dev/null | bwa samse '$index' /dev/stdin $fastq > $sam 2> /dev/null`;
    `samtools view -Shb $sam 2> /dev/null | samtools sort /dev/stdin | samtools view -h  > $sorted_sam_accepted`;
    unlink(@garbage);
    open(my $fic,$sorted_sam_accepted) or die "Cannot open $sam file: $!\n";
    open(my $accepted, ">".$fastq_accepted) || die "Cannot open $fastq_accepted: $! \n";
    open(my $unique, ">".$fastq_accepted_unique) || die "Cannot open $fastq_accepted_unique: $! \n";
    open(my $rejected, ">".$fastq_rejected) || die "Cannot open $fastq_rejected: $! \n";
    my $sequence = '';
    while(<$fic>)
    {
        chomp $_;
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                my $N = $1;
                $N =~ tr/\|:\//___/;
                $size{$N} = $2;
                $unique_number{$N} = 0;
                $number{$N} = 0;
                $numberNM{$N} = 0;
                $numberM{$N} = 0;
            }
            next;
        }
        my @line = split (/\t/,$_);
        $line[2] =~ tr/\|:\//___/;   
        $sequence = $line[9];
        if ($line[1] & 16)
        {
            $sequence =reverse($sequence);
            $sequence =~ tr/atgcuATGCU/tacgaTACGA/;
        }
        if ($line[1] == 16 || $line[1] == 0)
        {
            my $len = length($sequence);
            $size_num{$len} ++;
            $size_num_spe{$line[2]}{$len}++;
            $mappers ++;
            ${$best_hit_number_hashR}{$sequence} = $1  if  ($line[13] =~ /X0:i:(\d*)/ ||  $line[14] =~/X0:i:(\d*)/ );
            ${$duplicate_hashR}{$sequence}++;
            $number{$line[2]}++;
            $numberSens{$line[2]}++ if $line[1] == 0 ;
            $numberReverse{$line[2]}++ if $line[1] == 16 ;
            print $accepted "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
            if ($line[11] eq "XT:A:U")
            {
                $unique_number{$line[2]}++;
                $mappersUnique++;
                print $unique "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
            }
            if ($_ =~ /.*XM:i:(\d+).*/)
            {
                if ($1 == 0)
                {
                    $numberNM{$line[2]}++;
                }
                else
                {
                    $numberM{$line[2]}++;
                }
            }
        }
        else
        {
            ${$best_hit_number_hashR}{$sequence} = 0;
            ${$duplicate_hashR}{$sequence}++;
            print $rejected "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
        }
    }
    close $fic; close $accepted; close $unique; close $rejected;
    print STDERR "-----------------------------\n";
    print STDERR "$index\n$fastq\n";
    print STDERR "\treads: $reads\n";
    print STDERR "\tmappers: $mappers\n";
    print STDERR "\tunique mappers: $mappersUnique\n";
    unlink $sam;
    return (\%number, \%unique_number, \%numberSens, \%numberReverse, \%size, $mappers, $mappersUnique, \%size_num, \%size_num_spe, \%numberNM,\%numberM );
}

#input $Hash ref keys = abscissa values = ordonate -png file name
sub histogram 
{
    my ($size_hashR, $out_png, $size)  = @_;
    my (@abs, @ord);
    my $i = 0;
    foreach my $k (sort {$a <=> $b} keys %{$size_hashR})
    {
        my $percentage = $size_hashR->{$k} * 100 / $size;
        $abs[$i] = $k ; $ord[$i] = $percentage; $i++;
    }
    my $abs = join (",", @abs );
    my $ord = join (",", @ord );
    if (scalar(@abs) != 0)
    {
        my $R = Statistics::R->new();
        $R->startR;
        $R->send(
            qq`
            library(ggplot2)
            percentage <- c($ord)
            size <- c($abs)
            min <- min(size) - 1
            max <- max(size) + 1
            dat <- data.frame(size,percentage)
            c <- ggplot(dat,aes(size,percentage))+ geom_bar(stat="identity") + scale_x_continuous(breaks=min:max)
            ggsave(filename="$out_png", dpi=600)
            `);
        $R->stopR();
    }
}

sub histogram_Gviz 
{
    my ($name, $dir, $chro, $start, $end, $sens, $anti) = @_;
    my $svg = $dir."/".$name.".jpg";
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(qq `library(grid)`);
    $R->send(qq `library(Gviz)`);
    $R->send(qq `jpeg(filename="$svg", width = 480, height = 360,quality = 100)`);
    $R->send(qq `options(ucscChromosomeNames=FALSE)`);
    $R->send(qq `chromosome = c($chro)`);
    $R->send(qq `start = c($start)`);
    $R->send(qq `end = c($end)`);
    $R->send(qq `sens = c($sens)`);
    $R->send(qq `anti = c($anti)`);
    $R->send(qq `anti = -anti`);
    $R->send(qq `T1 = data.frame(chromosome,start,end,sens,anti)`);
    $R->send(qq `names(T1) =c("chromosome","start","end","sens","anti")`);
    $R->send(qq `Track = DataTrack(T1,background.title = "black",genome= "dmel", name = "$name",type= "hist", groups = c("sens","anti"))`);
    $R->send(qq `axis = GenomeAxisTrack()`);
    $R->send(qq `plotTracks(list(Track,axis))`);
    $R->send(qq `dev.off()`);
    $R->stopR();
}

sub rpms_rpkm
{
    my ($counthashR, $sizehashR, $mapped, $out_file, $piRNA_number, $miRNA_number, $bonafide_number) =@_;
    open(my $out, ">".$out_file) or die "Cannot open normalized file: $! \n";
    print $out "ID\treads counts\tRPKM";
    print $out "\tper million of piRNAs" if ($piRNA_number != 0);
    print $out "\tper million of miRNAs" if ($miRNA_number != 0);
    print $out "\tper million of bonafide reads" if ($bonafide_number != 0);
    print $out "\n";
    
    while( my ($k, $v) = each %{$counthashR})
    {
        my ($rpkm, $pirna, $mirna, $bonafide) = (0,0,0,0);
        $rpkm = ($v * 1000000000) / (${$sizehashR}{$k} * $mapped) if (${$sizehashR}{$k} * $mapped) != 0;
        print $out $k."\t".$v."\t"; printf $out "%.2f",$rpkm;
        if ($piRNA_number != 0 )
        {
            $pirna = ($v  * 1000000) / $piRNA_number;
            printf $out "\t%.2f",$pirna;
        }
        if ($miRNA_number != 0 )
        {
            $mirna = ($v  * 1000000) / $miRNA_number;
            printf $out "\t%.2f",$mirna;
        }
        if ($bonafide_number != 0 )
        {
            $bonafide = ($v  * 1000000) / $bonafide_number;
            printf $out "\t%.2f",$bonafide;
        }
        print $out "\n";
    }
    close $out;
}

sub size
{
    my ($min, $max, $in_file, $out_file, $sizeHashR, $duplicateHashR) = @_;
    my ($numreads, $size, $cmp, $ok, $line) = (0, 0, 0, 0, 0);
    my @fastq;
    open (my $in, $in_file) or die "Cannot open fastq file: $!\n";
    open (my $out, ">".$out_file) or die "Cannot open output fastq file: $!\n";
    while(<$in>)
    {
        chomp $_;
        $cmp++; $line++;
        if ($cmp == 1)
        {
            die "file do not contain a @ at line $line\n" unless ($_ =~ /^\@/ );
            $ok = 0; @fastq = ();
            push(@fastq,$_);
        }
        elsif ($cmp == 2)
        {
            die "unrecognized symbol at line $line\n" unless ($_ =~ /[atcgATCGnN]+/ || $_ =~ /^$/ );
            push(@fastq,$_);
            $size = length($_);
            if ($size >= $min && $size <= $max)
            {
                $numreads++;
                ${$sizeHashR}{$size}+=1;
                ${$duplicateHashR}{$_}+=1 if (defined($duplicateHashR));
                $ok = 1;
            }
        }
        elsif ($cmp == 3 )
        {
            die "file do not contain a + at line $line\n" unless $_ =~ /^\+/;
            push(@fastq,$_);
        }
        elsif ($cmp == 4 )
        {
            push(@fastq,$_);
            $cmp = 0;
            if ($ok == 1)
            {
                foreach my $t (@fastq)
                {
                    print $out $t."\n";
                }
            }
         }
    }
    close $in; close $out;
    return $numreads;
}

sub windows
{
    my ($sam, $windows, $dir, $uni, $dirU, $max_procs, $mapIn) = @_;
    my $pm = Parallel::ForkManager->new($max_procs);
    my (%Sens, %Anti, %stepK, $step, $first_window,$interval_number, $last_window, $mappers, $p, $pu);
    my (%uSens, %uAnti, $mappersU) if $uni == 1;
    $mappers = 0;
    $mappersU = 0 if $uni ==1;
    open(my $s,$sam) or die "Cannot open $sam file: $!\n";
    while(<$s>)
    {
        chomp $_; 
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                my $N = $1;
                $N =~ tr/\|:\//___/;
                if ($windows eq 'auto')
                {
                    $step = floor (sqrt ($2));
                }
                else { $step = floor ($2 / $windows); }
                if ($step != 0 )
                {
                    $interval_number = floor(($2-1) / $step); $stepK{$N} = $step
                }
                else {$interval_number = 0; $stepK{$1} = $2}
                $Sens{$N} =[]; $#{$Sens{$N}} = $interval_number;
                $Anti{$N} =[]; $#{$Anti{$N}} = $interval_number;
                for my $i (0 .. $#{$Sens{$N}})
                {
                    ${$Sens{$N}}[$i] = 0;
                    ${$Anti{$N}}[$i] = 0;
                }
                if ($uni == 1)
                {
                    $uSens{$N} =[]; $#{$uSens{$N}} = $interval_number;
                    $uAnti{$N} =[]; $#{$uAnti{$N}} = $interval_number;
                    for my $i (0 .. $#{$uSens{$N}})
                    {
                        ${$uSens{$N}}[$i] = 0;
                        ${$uAnti{$N}}[$i] = 0;
                    }
                }
            }
            next;
        }
        my @line =split(/\t/,$_);
        if ($line[1] == 0 or $line[1] == 16)
        {
            $line[2]=~ tr/\|:\//___/;
            $first_window = 0; $last_window = 0;
            $mappers++;
            $mappersU++ if ($uni == 1 && $line[11] eq "XT:A:U");
            $first_window = floor( ($line[3]-1)  / $stepK{$line[2]}) unless $stepK{$line[2]} == 0;
            $last_window = floor( ($line[3] + length ($line[9]) - 2) / $stepK{$line[2]}) unless $stepK{$line[2]} == 0;
            if ($line[1] == 0)
            {
                $p = $Sens{$line[2]};
                $pu = $uSens{$line[2]} if $uni == 1;
            }
            else
            {
                $p = $Anti{$line[2]};
                $pu = $uAnti{$line[2]} if $uni == 1;
            }
            $first_window = 0 if $first_window < 0;
            $last_window = $#{$p} if $last_window >  $#{$p};
            for (my $i = $first_window; $i <= $last_window; $i++)
            {
                ${$p}[$i] += 1;
                ${$pu}[$i] += 1 if ($uni == 1 && $line[11] eq "XT:A:U");
            }
        }
    }
    close $sam; 
    if (defined($mapIn) && $mapIn != 0 )
    {
        $mappers = $mapIn;
        $mappersU = $mapIn if $uni == 1;
    } 
    foreach my $k (keys %Sens)
    {
        $pm->start() and next;
        my ($Stot, $Atot, $S, $A, $begin, $end, $rpkmS, $rpkmA, $SUtot, $AUtot) =(0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0);
        my (@Hchromosome, @Hstart, @Hend, @Hsens, @Hanti);
        my (@HsensU, @HantiU, $rpkmSU, $rpkmAU, $AU, $SU);
        #my $k2 = $k;
        #$k =~ tr/\|:\//___/;
        open (my $out_uni, '>'.$dirU.'/'.$k.'.bed');
        open (my $out, '>'.$dir.'/'.$k.'.bed');
        $Stot = 0; $Atot = 0; $SUtot = 0; $AUtot = 0;
        for (my $t = 0; $t <= $#{$Sens{$k}}; $t++)
        {
            $S = 0; $A = 0;$rpkmS = 0; $rpkmA = 0;
            $begin = $t * $stepK{$k} + 1;
            $end = ($t + 1) * $stepK{$k};
            push(@Hchromosome, $k);
            push(@Hstart, $begin);
            push(@Hend, $end);
            $S = $Sens{$k}->[$t]; $Stot += $S;
            $A = $Anti{$k}->[$t]; $Atot += $A;
            $rpkmS = ($S * 1000000000) / ($mappers *  $stepK{$k} ) if ($mappers *  $stepK{$k} ) != 0 ;
            $rpkmA = ($A * 1000000000) / ($mappers *  $stepK{$k}) if  ($mappers *  $stepK{$k} ) != 0 ;
            push(@Hsens,$rpkmS); push(@Hanti,$rpkmA);
            print $out "$k\t$begin\t$end\t$rpkmS\t$rpkmA\n";
            if ( $uni == 1)
            {
                $SU = 0; $AU = 0; $rpkmSU = 0; $rpkmAU = 0;
                $SU = $uSens{$k}->[$t]; $SUtot += $SU;
                $AU = $uAnti{$k}->[$t]; $AUtot += $AU;
                $rpkmSU = ($SU * 1000000000) / ($mappersU * $stepK{$k}) if  ($mappersU *  $stepK{$k} ) != 0 ;
                $rpkmAU = ($AU * 1000000000) / ($mappersU * $stepK{$k}) if ($mappersU *  $stepK{$k} ) != 0 ;
            push(@HsensU,$rpkmSU); push(@HantiU,$rpkmAU);
            print $out_uni "$k\t$begin\t$end\t$rpkmSU\t$rpkmAU\n"; 
            }
        }
        @Hchromosome = map { "'$_'" } @Hchromosome;
        my $Hchromosome = join ',', @Hchromosome;
        my $Hstart = join ',', @Hstart;
        my $Hend = join ',', @Hend;
        my $Hsens = join ',', @Hsens;
        my $Hanti = join ',', @Hanti;
        if ($Stot != 0 or $Atot != 0)
        {
            histogram_Gviz($k, $dir, $Hchromosome, $Hstart, $Hend, $Hsens, $Hanti );
            if ($uni==1)
            {
                if ($SUtot != 0 or $AUtot != 0)
                {
                    my $HsensU = join ',', @HsensU;
                    my $HantiU = join ',', @HantiU;
                    histogram_Gviz($k, $dirU, $Hchromosome, $Hstart, $Hend, $HsensU, $HantiU );
                }
                else{unlink $dirU.'/'.$k.'.bed';}
            }
        }
        else{unlink $dir.'/'.$k.'.bed'; unlink $dirU.'/'.$k.'.bed' if $uni == 1;}
        $pm->finish(); # pass an exit code to finish
    }
    $pm->wait_all_children;
}
1;
