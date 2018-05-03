#!/usr/bin/env perl

use strict ;
use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;
use Digest::MD5  qw(md5 md5_hex md5_base64) ;

my $i ;
my @readFileList ;
my $numOfThread = 1 ;
my $kmerSize = 23 ;
my $bloomFilterSize = 100000000 ;
my $WD = dirname( abs_path( $0 ) ) ;
my $stage = 0 ;

my $usage = "Usage: perl ./run_rcorrector.pl [OPTIONS]\n".
		"OPTIONS:\n".
		"Required parameters:\n".
		"\t-s seq_files: comma separated files for single-end data sets\n".
		"\t-1 seq_files_left: comma separated files for the first mate in the paried-end data sets\n".
		"\t-2 seq_files_right: comma separated files for the second mate in the paired-end data sets\n".
		"\t-i seq_files_interleaved: comma sperated files for interleaved paired-end data sets\n".
		"Other parameters:\n".
		"\t-k kmer_length (<=32, default: 23)\n".
		"\t-od output_file_directory (default: ./)\n".
		"\t-t number_of_threads (default: 1)\n".
		#"\t-trim allow trimming (default: false)\n".
		"\t-maxcorK INT: the maximum number of correction within k-bp window (default: 4)\n".
		"\t-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold, lower for more divergent genome (default: 0.95)\n".
		"\t-ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)\n".
		"\t-stdout: output the corrected reads to stdout (default: not used)\n".
		"\t-verbose: output some correction information to stdout (default: not used)\n".
		"\t-stage INT: start from which stage (default: 0)\n".
		"\t\t0-start from begining(storing kmers in bloom filter);\n".
		"\t\t1-start from count kmers showed up in bloom filter;\n".
		"\t\t2-start from dumping kmer counts into a jf_dump file;\n".
		"\t\t3-start from error correction." ;

if ( scalar( @ARGV ) == 0 )
{
	die "$usage\n" ;
}

my $jellyfishBin = "jellyfish" ;
if ( -e "$WD/jellyfish/bin/jellyfish" )
{
	$jellyfishBin = "$WD/jellyfish/bin/jellyfish" ;
}

my $fileArguments ;
my @singleFileList ;
my @firstMateFileList ;
my @secondMateFileList ;
my @interleavedFileList ;

my @rcorrectorArguments ;
my @gzippedFileList ;

sub AddJellyFishReadFile
{
	if ( $_[0] =~ /gz$/ )
	{
		push @gzippedFileList, $_[0] ;
	}
	else
	{
		push @readFileList, $_[0] ;
	}
}

for ( $i = 0 ; $i < scalar(@ARGV) ; ++$i )
{
	if ( $ARGV[$i] eq "-r" )
	{
		AddJellyFishReadFile( $ARGV[ $i + 1] ) ;

		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-p" )
	{
		AddJellyFishReadFile( $ARGV[ $i + 1 ] ) ;
		AddJellyFishReadFile( $ARGV[ $i + 2 ] ) ;

		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;
		push @rcorrectorArguments, $ARGV[$i + 2] ;

		$i += 2 ;
	}
	elsif ( $ARGV[ $i ] eq "-s" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			AddJellyFishReadFile( $cols[ $j ] ) ;
			push @singleFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-1" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			AddJellyFishReadFile( $cols[ $j ] ) ;
			push @firstMateFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[ $i ] eq "-2" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			AddJellyFishReadFile( $cols[ $j ] ) ;
			push @secondMateFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-i" )
	{
		my @cols = split /,/, $ARGV[$i + 1] ;
		my $j ;
		for ( $j = 0 ; $j < scalar( @cols ) ; ++$j )
		{
			AddJellyFishReadFile( $cols[ $j ] ) ;
			push @interleavedFileList, $cols[ $j ] ;
		}
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-t" )
	{
		$numOfThread = $ARGV[$i+1] ;

		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-k" )
	{
		$kmerSize = $ARGV[$i + 1] ;

		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-ek" )
	{
		$bloomFilterSize = $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-maxcor" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-maxcorK" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;

		++$i ;
	}
	elsif ( $ARGV[$i] eq "-stdout" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
	}
	elsif ( $ARGV[$i] eq "-verbose" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
	}
	elsif ( $ARGV[$i] eq "-stage" )
	{
		$stage = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-od" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-wk" )
	{
		push @rcorrectorArguments, $ARGV[$i] ;
		push @rcorrectorArguments, $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown argument ".$ARGV[$i]."\n" ;
	}
}
#`echo $numOfThread > tmp.out `

if ( $kmerSize > 32 )
{
	die "The kmer size can not be greater than 32.\n" ;
}

# Build the input file arguments for
for ( my $i = 0 ; $i < @singleFileList ; ++$i )
{
	$fileArguments = $fileArguments." -r ".$singleFileList[$i] ;
}


die "The number of files from -1,-2 should be the same" if ( scalar( @firstMateFileList ) != scalar( @secondMateFileList ) ) ;
for ( my $i = 0 ; $i < @firstMateFileList ; ++$i )
{
	$fileArguments = $fileArguments." -p ".$firstMateFileList[$i]." ".$secondMateFileList[$i] ;
}

for ( my $i = 0 ; $i < @interleavedFileList ; ++$i )
{
	$fileArguments = $fileArguments." -i ".$interleavedFileList[$i] ;
}

# build the file list for jellyfish
my $jellyfishFiles = "" ;
for ( my $i = 0 ; $i < @readFileList ; ++$i )
{
	$jellyfishFiles .= $readFileList[$i]." " ;
}

for ( my $i = 0 ; $i < @gzippedFileList ; ++$i )
{
	$jellyfishFiles .= "<(gzip -cd ".$gzippedFileList[$i].") " ;
}

my $crc = md5_hex( $jellyfishFiles ) ;

if ( $stage <= 0 )
{
	print STDERR ( "Put the kmers into bloom filter\n" ) ;
	print STDERR ( "$jellyfishBin bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp_$crc.bc $jellyfishFiles\n" ) ;
	die "Failed at stage 0.\n" if ( system( "bash -c \"$jellyfishBin bc -m $kmerSize -s $bloomFilterSize -C -t $numOfThread -o tmp_$crc.bc $jellyfishFiles\"" ) != 0 ) ;
}

if ( $stage <= 1 )
{
	print STDERR ( "Count the kmers in the bloom filter\n" ) ;
	print STDERR ( "$jellyfishBin count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp_$crc.bc -o tmp_$crc.mer_counts $jellyfishFiles\n" ) ;
	die "Failed at stage 1.\n" if ( system( "bash -c \"$jellyfishBin count -m $kmerSize -s 100000 -C -t $numOfThread --bc tmp_$crc.bc -o tmp_$crc.mer_counts $jellyfishFiles\"" ) != 0 ) ;
}

if ( $stage <= 2 )
{
	print STDERR ( "Dump the kmers\n" ) ;
	print STDERR ( "$jellyfishBin dump -L 2 tmp_$crc.mer_counts > tmp_$crc.jf_dump\n" ) ;
	die "Failed at stage 2.\n" if ( system( "$jellyfishBin dump -L 2 tmp_$crc.mer_counts > tmp_$crc.jf_dump" ) != 0 )
}

if ( $stage <= 3 )
{
	print STDERR ( "Error correction\n" ) ;
	print STDERR ( "$WD/rcorrector @rcorrectorArguments $fileArguments -c tmp_$crc.jf_dump\n" ) ;
	die "Failed at stage 3.\n" if ( system( "$WD/rcorrector @rcorrectorArguments $fileArguments -c tmp_$crc.jf_dump" ) != 0 ) ;
}

system( "rm tmp_$crc.bc tmp_$crc.mer_counts tmp_$crc.jf_dump" );
