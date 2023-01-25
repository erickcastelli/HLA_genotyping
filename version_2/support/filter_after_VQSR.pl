#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Basename;

our($opt_r, $opt_v);
getopts('r:v:');


if (($opt_r eq "") || (! -e $opt_r)) {help();}
if (($opt_v eq "") || (! -e $opt_v)) {help();}

print "Loading reference...\n";
my $base = $opt_r;
my %db = ();
open (IN, $base) or die;
while (<IN>)
{
	chomp $_;
	if (substr($_,0,1) eq "#") {next;}
	my @data = split("\t",$_);
	$db{$data[1]} = $data[4];
}
close (IN);



print "Processing vcf...\n";
my $vcf = $opt_v;
my $dir = dirname($vcf);
my $base = basename($vcf);
$base =~ s/\.vcf/\.filter.vcf/;
my $out = $dir . "/" . $base;

my $forced = 0;
open (IN, $vcf) or die;
open (OUT, ">$out") or die;
while (<IN>)
{
	chomp $_;
	if (substr($_,0,2) eq "##") {print OUT $_ . "\n"; next;}
	if (substr($_,0,2) eq "#C") 
	{
		print OUT $_ . "\n"; 
		next;
	}
	my $line = $_;
	my @data = split("\t",$line);
	if ($data[6] =~ /PASS/) {print OUT $line . "\n"; next;}
	
	if (exists $db{$data[1]}) {
		my %alts = ();
		my @base_alt = split(",",$db{$data[1]});
		foreach(@base_alt) {$alts{$_} = 1;}
		my @select_alt = split(",",$data[4]);
		my $help = 0;
		foreach (@select_alt)
		{
			if ($_ eq "*") {next;}
			if ($_ eq ".") {next;}
			if (exists $alts{$_}) {$help = 1; last;}
		}
		
		if ($help eq 1) 
		{
			print "Position $data[1] forced to PASS ...$forced\n";
			$data[6] = "PASS";
			my $new = join("\t",@data);
			print OUT $new . "\n";
			$forced++;
			next;
		}
		next;
	}
}
close (IN);
close (OUT);


print "Forced: $forced\n";

sub help()
{
	print "-v vcf_to_be_processed -r reference_vcf (usually, MHC.all.vcf)\n";
	exit;
}