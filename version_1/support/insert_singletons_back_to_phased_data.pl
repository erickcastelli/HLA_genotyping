#!/usr/bin/perl
#Code by Erick C. Castelli
#Version 1.0

use strict;
use Getopt::Std;

our ($opt_l, $opt_p, $opt_r, $opt_o);

getopts('l:p:r:o:');

if (! -e $opt_l) {help();}
if (! -e $opt_p) {help();}
if (! -e $opt_r) {help();}
if ($opt_o eq "") {help();}

my $log = $opt_l;
my $rbp = $opt_r;
my $phased = $opt_p;
my $out = $opt_o;

my @logdata = ();
open (IN, $log);
while (<IN>)
{
	chomp $_;
	my @data = split("\t",$_);
	push(@logdata,$_);
}
close (IN);



my @rbp_data = ();
open (IN, $rbp);
while (<IN>)
{
	chomp $_;
	if ($_ eq "") {next;}
	if (substr($_,0,2) eq "##") {next;}
	push(@rbp_data,$_);
}
close (IN);




my %phased_data = ();
my %snps = ();
my @samples = ();
my @head = ();
my @sample_db = ();
my %sample_hash = ();

open (IN, $phased);
while (<IN>)
{
	chomp $_;
	if ($_ eq "") {next;}
	if (substr($_,0,2) eq "##") {push(@head,$_);next;}
	if (substr($_,0,2) eq "#C") 
	{
		push(@head,$_);
		@samples = split("\t",$_);
		for (my $a = 9; $a < scalar(@samples); $a++)
		{
			push(@sample_db, $samples[$a]);
			$sample_hash{$samples[$a]} = 1;
		}
		next;
	}
	my @data = split("\t",$_);
	my $pos = $data[1];
	my $snp = join("\t",@data[0..8]);
	$snps{$pos} = $snp;
	for (my $a = 9; $a < scalar(@data); $a++)
	{
		$phased_data{$pos}{$samples[$a]} = $data[$a];
	}
}
close (IN);



foreach (@logdata)
{
	my @data = split("\t",$_);
	my $pos = $data[0];

	my $sing = $data[1];
	my $sample = $data[3];
	if ($sample_hash{$sample} ne 1) {next;}
	print "$pos\t$sample\tallele $sing\t";
	
	my @samples = ();
	my $index = 0;
	my $info = "";
	my $var = "";
	my $relate_to = "";
	
	
	foreach (@rbp_data)
	{
		if (substr($_,0,2) eq "#C") 
		{
			@samples = split("\t",$_);
			for ($index = 0; $index < scalar(@samples);$index++)
			{
				if ($samples[$index] eq $sample) {last;}
			}
			next;
		}
		my @data = split("\t",$_);
		if ($data[1] ne $pos) {next;}
		$info = $data[$index];
		if ($info ne "")
		{
			my @tmp = split(":",$info);
			@tmp = split(",",$tmp[4]);
			if (scalar(@tmp) eq 2) 
			{
				my @tmp2 = split("-",$tmp[0]);
				$relate_to = $tmp2[0];
			}
		}
		$var = join("\t",@data[0..7]);
		last;		
	}

	my @vectors = (); 
	foreach (@rbp_data)
	{
		if (substr($_,0,1) eq "#") 
		{
			next;
		}
		my @data = split("\t",$_);
		if ($data[1] eq $pos) {next;}
		if ($data[$index] =~ /$relate_to/)
		{
			push(@vectors,$data[1] . ":" . $data[$index]);
		}
	}
	if (($relate_to eq "") || ($relate_to eq ".")) {@vectors = ();}

	my $status = "";
	if (scalar(@vectors) eq 0) {$status = "unphased";}
	if (scalar(@vectors) > 0) {$status = "phased";}
	if ($info eq "") {$status = "no data"};	
	if (exists $phased_data{$pos}) {$status .= ",exist";}	
	if (! exists $phased_data{$pos}) {$status .= ",new";}

	
	if ($status =~ /no data/) {print "$status\n";next;}

	if ($status =~ /new/)
	{
		$snps{$pos} = $var . "\tGT";
		foreach (@sample_db) {$phased_data{$pos}{$_} = "0|0";}
	}

	
	if ($status eq "unphased,exist") {
		my @infos = split(":",$info);
		my $gen = $infos[0];
		$phased_data{$pos}{$sample} = $gen;
	}
	
	if ($status eq "phased,new") {
		my @infos = split(":",$info);
		my @alleles = split("/",$infos[0]);
		my $gen = "";
		if ($infos[4] eq "$relate_to\-1,$relate_to\-2") {$gen = $alleles[0] . "|" . $alleles[1]};
		if ($infos[4] eq "$relate_to\-2,$relate_to\-1") {$gen = $alleles[1] . "|" . $alleles[0]};
		my $genvector = "";
		my $vector = "";
		foreach (@vectors)
		{
			my @infos = split(":",$_);
			$vector = $infos[0];
			my @alleles = split("/",$infos[1]);
			if ($infos[5] eq "$relate_to\-1,$relate_to\-2") {$genvector = $alleles[0] . "|" . $alleles[1]};
			if ($infos[5] eq "$relate_to\-2,$relate_to\-1") {$genvector = $alleles[1] . "|" . $alleles[0]};
			last;
		}	
		
		my $ref = $phased_data{$vector}{$sample};
		
#		print "\t$ref\t$gen\t$genvector\t$vector\t";
		
		if ($genvector eq $ref) {
			$phased_data{$pos}{$sample} = $gen;
		}
		if ($genvector ne $ref) {
			$gen = reverse $gen;
			$phased_data{$pos}{$sample} = $gen;
		}	
	}
	
	if ($status eq "unphased,new") {
		my @infos = split(":",$info);
		my $gen = $infos[0];
		$phased_data{$pos}{$sample} = $gen;
	}
		
	print $status . "\n";
#	last;
	
}

open (OUT, ">$out");
foreach (@head)
{
	print OUT $_ . "\n";
}

my @positions = sort {$a <=> $b} keys %snps;
foreach (@positions)
{
	print OUT $snps{$_};
	my $pos = $_;
	foreach (@sample_db)
	{
		print OUT "\t$phased_data{$pos}{$_}";
	}
	print OUT "\n";
}






sub help {
	print "\nPerl script to reintroduce unphased singletons back to a phased VCF\n";
	print "How to use it:\n";
	print "-l singleton log file (produced by remove_unphased_singleton_from_normalized_rbp.pl)\n";
	print "-p phased VCF (the results.vcf from Phasex)\n";
	print "-r VCF from RBP (the original VCF file outputted by ReadBackedPhasing)\n";
	print "-o output VCF (the VCF file to be generated)\n\n";
	exit;
}
