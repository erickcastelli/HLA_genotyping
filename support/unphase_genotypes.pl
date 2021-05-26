#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Basename;

our ($opt_i,$opt_o);
getopts('i:o:'); 

if (($opt_i eq "") || (! -e $opt_i)) {print "\nInput file not detected!\n"; help();}
if (($opt_o eq "") || (! -e dirname $opt_o)) {print "\nOutput folder not detected!\n"; help();}

my $in = $opt_i;
my $out = $opt_o;

open (IN, $in);
open (OUT, ">$out");

while (<IN>)
{
	if (substr($_,0,1) eq "#") {print OUT $_;next;}
	my $line = $_;
	$line =~ s/\|/\//g;
	print OUT $line;
}
close (IN);
close (OUT);

sub help {
	print "\n\nPerl script to convert | to / in VCF files\n";
	print "How to use it:\nperl $0 -i input_vcf -o output_vcf\n\n";
	exit;
}
