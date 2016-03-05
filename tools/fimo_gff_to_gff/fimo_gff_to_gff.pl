#! /usr/bin/perl

die "FIMO_GFF_File\tOutput_Path\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";

##gff-version 3
#chr10:265210-265270(-)  fimo    nucleotide_motif        25      36      40.2    +       .       Name=1;ID=1-1-chr10:265210-265270(-);pvalue=9.48e-05;qvalue=0.00885;sequence=ACTTACCCTCAT;
#chr10:295039-295099(+)  fimo    nucleotide_motif        25      36      55.3    +       .       Name=1;ID=1-1-chr10:295039-295099(+);pvalue=2.97e-06;qvalue=0.00107;sequence=TGTTACCCGTTC;
#chr10:576747-576807(-)  fimo    nucleotide_motif        25      36      56.2    +       .       Name=1;ID=1-1-chr10:576747-576807(-);pvalue=2.37e-06;qvalue=0.00107;sequence=CGTTACCCGACC;

#chr1        genetrack        .        123950        123970        22        +        .        stddev=0.0
#chr1        genetrack        .        565745        565765        12        +        .        stddev=0.0
#chr1        genetrack        .        565793        565813        44        +        .        stddev=0.298065387468

@COORD = ();
@ID_NUM = ();
$line = "";
while($line = <IN>) {
        chomp($line);
        next if($line =~ /gff-version/);
        @array = split(/\t/, $line);
        @CHR = split(/\:/, $array[0]);
        @gff_COORD = split(/\(/, $CHR[1]);
        @START_array = split(/\-/, $gff_COORD[0]);
        $fimo_DIR = "+";
        if($gff_COORD[1] =~ "-") { $fimo_DIR = "-"; }

        $DIR = $array[6];
        $SCORE = $array[5];

        @NAME = split(/\;/, $array[8]);
        $NEW = 0;
        for($x = 0; $x <= $#ID_NUM; $x++) {
                if($ID_NUM[$x] eq $NAME[0]) {
                        $NEW = 1;
                        $x = $#ID_NUM + 1;
                }
        }
        if($NEW == 0) { push(@ID_NUM, $NAME[0]); }

        $START = $START_array[0] + $array[3] - 1;
        $STOP = $START_array[0] + $array[4] - 1;

        if($fimo_DIR eq "-") {
                if($DIR eq "+") { $DIR = "-"; }
                else { $DIR = "+"; }
        }

        $newline = "$CHR[0]\tfimo\tmotif\t$START\t$STOP\t$SCORE\t$DIR\t.\t$CHR[0]\_$START\_$STOP\_$DIR";
        $EXISTS = 0;
        for($x = 0; $x <= $#COORD; $x++) {
                if($newline eq $COORD[$x]{'line'}) {
                        $EXISTS = 1;
                }
        }
        if($EXISTS == 0) {
                push(@COORD, {chr => $CHR[0], start => $START, stop => $STOP, dir => $DIR, score =>$SCORE, id => $NAME[0], line => $newline});
        }
}
close IN;
@SORT = sort { $$b{'score'} <=> $$a{'score'} } @COORD;

for($x = 0; $x <= $#ID_NUM; $x++) {
        @FILENAME = split(/\=/, $ID_NUM[$x]);
        $FILE = "MOTIF$FILENAME[1]";
        open(OUT, ">$output/$FILE.gff") or die "Can't open $output/$FILE.gff for writing!\n";
        for($y = 0; $y <= $#SORT; $y++) {
                if($SORT[$y]{'id'} eq $ID_NUM[$x]) {
                        print OUT $SORT[$y]{'line'},"\n";
                }
        }
        close OUT;
}
