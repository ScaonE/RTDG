#!/usr/bin/perl
use warnings; use strict;
use Getopt::Long;

if (@ARGV < 1) { die "$!\n\tUsage\t: assemblyStats.pl contig/scaffs.fa(single-line) [Options]\n\tOptions\t:\n\t\t-genome\t<int>\t: Estimated genome size (bp) to compute NG50|90\n\t\t-size1\t<int>\t: Will reports the number of sequences >= size1 (bp)\n\t\t-size2\t<int>\t: Will reports the number of sequences >= size2 (bp)\n\t\t-graph\t\t: Produces 2 files for R plots\n\tOutput\t:\n\t\tsequences\n\t\tbases\n\t\tmedian\n\t\tmean\n\t\tmin\n\t\tmax\n\t\tN50\n\t\tNG50\n\t\tGC content\n\t\tsequences with N's\n\t\tN's\n\t\tseq > size1\n\t\tseq > size2
\n"; }

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

### Initialize variables & retrieve options
my (%hash, @allLengths) = ()x2;
my ($sumNucl, $headers, $sumN, $seqN, $GC, $mean, $genome, $size1, $size2, $graph) = (0)x11;
GetOptions ("genome=i"	=> \$genome,
			"size1=i"	=> \$size1,
			"size2=i"	=> \$size2,
			"graph"		=> \$graph); # TODO : CREATE R FILES # Counter-TODO : QUAST does it nicely already

### Iterate once through the contigs/scaffs file
open(INPUT, "<", $ARGV[0]) or die "Error : $!\n";
	while(<INPUT>) {
		if ($_ =~ /^>/) { $headers ++; }
		else {
			chomp $_; $_ = (uc $_); my $lgt = length($_);
			$hash{$lgt} ++; # Hash keys <=> to all existing sequences lengths. Values <=> to the number of time sequences with such lengths occurs
			$sumNucl += $lgt; # Calc total number of nucleotides in assembly
			my $temp = ($_ =~ tr/N/N/);
			if ($temp) { $sumN += $temp; $seqN ++; } # If N(s) occurs within current sequence, sum it into $sumN and increment $seqN
			$GC += ($_ =~ tr/GC//); # Count all G|C in sequences
			my $delta = ($lgt-$mean); # This delta allow to stream compute the sequences mean
			$mean += ($delta/$headers);
		}
	}
close INPUT;

## Set NYXX thresholds <=> the number of bp to represent XX% of the assembly/genome
my $N50thresh = $sumNucl*((100-50)/100);
my $NG50thresh = $genome*((100-50)/100);
## Other variables
my $reachN50 = $sumNucl;
my $reachNG50 = $genome;
my ($N50, $NG50, $supSize1, $supSize2, $max, $min, $GCpct) = (0)x7;

foreach (sort {$b <=> $a} (keys(%hash))) {
## While $reachYXX values are > to their respective thresholds, decrement them
	if ($reachN50 >= $N50thresh) { $reachN50 -= ($_)*($hash{$_}); $N50 = $_; }
	if ($reachNG50 >= $NG50thresh) { $reachNG50 -= ($_)*($hash{$_}); $NG50 = $_; }
## Retrieve number of sequences > size1/2 & push every sequence length into an array for further median calcuation
	if ($_ >= $size1) { $supSize1 += $hash{$_}; }
	if ($_ >= $size2) { $supSize2 += $hash{$_}; }
	for (my $i=0; $i<$hash{$_}; $i++) {	push (@allLengths, $_); }
}
undef %hash;

## Compute some statistics
$max = $allLengths[0];
$min = $allLengths[-1];
my $median = median(@allLengths);
$GCpct = ($GC*100)/$sumNucl;
$mean = $mean;

print"$headers\t\tnumber of sequences\n$sumNucl\tnumber of total nts\n$median\t\tmedian seq length\n$mean\t\tmean\n$min\t\tmin\n$max\t\tmax\n$N50\t\tN50\n$GCpct\t\tGC%\n$seqN\t\tnumber of seq with N's\n$sumN\t\ttotal number of N's\n$supSize1\t\tnumber of seq > size1\n$supSize2\t\t> size2\n\n";
#print"$headers\n$sumNucl\n$median\n$mean\n$min\n$max\n$N50\n$NG50\n$GCpct\n$seqN\n$sumN\n$supSize1\n$supSize2\n";

__END__
#############################################################################
### Computes basic statistics for a contigs/scaffs.fasta file without ref ###
### ~~~ Erwan SCAON - IRISA Rennes - Latest modification : 01.03.2013 ~~~ ###
#############################################################################
