#!/usr/bin/perl
#Code by Erick C. Castelli
#Version 1.0

use Getopt::Std;
use strict;

our ($opt_i);
getopts('i:');

my $file = $opt_i;

if (! -e $file) {help();}

my $out = substr($file,0,length($file)-4) . ".nosingleton.vcf";
my $sing = substr($file,0,length($file)-4) . ".unphased_singletons.vcf";
my $log = substr($file,0,length($file)-4) . ".singletons.log";


my %snp = ();
my %vectors= ();
my %phased = ();
my %sample_data = ();
my @samples = ();

open (IN, $file);
while (<IN>)
{
	chomp $_;
	my $line = $_;
	if ($_ eq "") {next;}
	if (substr($_,0,2) eq "##") {next;}
	if (substr($_,0,2) eq "#C") {@samples = split("\t",$_); next;}
	my @data = split("\t",$_);

	my $key = "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]";
	
	for (my $a = 9; $a < scalar(@data); $a++)
	{
		my @fields = split(":",$data[$a]);
		my @alleles = ();
		my @phase_status = split("-",$fields[4]);
		my $phase = "";
		if ($fields[4] eq ".") {$phase = ".";} else {$phase = $phase_status[0];}
		
		if ($fields[0] =~ /\//) {@alleles = split("/",$fields[0]);}
		if ($fields[0] =~ /\|/) {@alleles = split("\\|",$fields[0]);}
		$snp{$key}{$alleles[0]}++;
		$snp{$key}{$alleles[1]}++;
		$vectors{$samples[$a]}{$phase}++;
		$vectors{$samples[$a]}{$phase}++;
		$phased{$key}{$alleles[0]} = $phase;
		$phased{$key}{$alleles[1]} = $phase;
		$sample_data{$key}{$alleles[0]}{$samples[$a]} = 1;
		$sample_data{$key}{$alleles[1]}{$samples[$a]} = 1;
	}
}
close (IN);


open (OUT, ">$out");
open (SING, ">$sing");
open (LOG, ">$log");

open (IN, $file);

while (<IN>)
{
	chomp $_;
	my $line = $_;
	if ($_ eq "") {next;}
	if (substr($_,0,2) eq "##") {print OUT $_ . "\n"; print SING $_ . "\n"; next;}
	if (substr($_,0,2) eq "#C") {print OUT $_ . "\n"; print SING $_ . "\n"; @samples = split("\t",$_); next;}
	my @data = split("\t",$_);
	
	my $alt = $data[4];
	my $snpid = $data[2];
	my $key = "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]";
	
	my @alleles = keys $snp{$key};
	my $status = 1000;
	my $info = "phased";
	foreach (@alleles)
	{
		if ($_ eq ".") {next;}
		if ($snp{$key}{$_} eq 1)
		{
			$status = $_;
			last;
		}
	}
	
	my $infosample = "";
	if ($status ne 1000)
	{
		my $vector_singleton = $phased{$key}{$status};
		if ($vector_singleton eq ".") 
		{
			my @sample = keys $sample_data{$key}{$status};
			$infosample = $sample[0];
			$info = "unphased";
		}
		if ($vector_singleton ne ".") 
		{
			my @sample = keys $sample_data{$key}{$status};
			$infosample = $sample[0];
			my $count = $vectors{$sample[0]}{$vector_singleton};
			if ($count > 1) {$info = "phased";} else {$info = "unphased";}
		}	
	}

	if ($status ne 1000)
	{
		print LOG $key . "\tsingleton\t$info\t$infosample\n"
	}
	

	if ($alt eq "\*") {
		my $new = join("\t",@data[0..8]);
		for (my $a = 9; $a < scalar(@data); $a++)
		{
			my @fields = split(":",$data[$a]);
			$fields[4] = ".";
			my $newdata = join(":",@fields);
			$new .= "\t$newdata"
		}
		$line = $new;
	}
	

	my @snpids = split(";",$snpid);
	if (scalar(@snpids) > 1)
	{
		@data = split("\t",$line);
		my $new = $data[0] . "\t" . $data[1] . "\t" . $snpids[0] . "\t" . $data[3] . "\t" . $data[4] . "\t" . $data[5] . "\t" . $data[6] . "\t" . $data[7] . "\t" . $data[8];
		for (my $a = 9; $a < scalar(@data); $a++)
		{
			my @fields = split(":",$data[$a]);
			$fields[4] = ".";
			my $newdata = join(":",@fields);
			$new .= "\t$newdata"
		}
		$line = $new;
	}
	
	
	if ($info eq "unphased")
	{
		print SING $line . "\n";;
	}
	else {
		print OUT $line . "\n";
	}
}
close (IN);

close OUT;
close (SING);



sub help()
{
	print "\nPerl script to detect and remove unphased singletons from VCF files produced by ReadBackedPhasing\n";
	print "How to use it:\nperl $0 -i vcf_file\n\n";
	exit;
}
