#!/usr/bin/perl
# from the hla-mapper workflow, calling SNPs and alleles for HLA genes
# Erick C. Castelli, erick.castelli@unesp.br
# version 1.0b

use strict;
use File::Basename;
use Getopt::Std;
use threads;
use Thread::Semaphore;


our ($opt_v, $opt_p, $opt_b, $opt_o, $opt_r, $opt_t, $opt_w, $opt_x);

getopts('v:p:b:o:r:t:w:x:');


my $vcf = $opt_v;
my $ped = $opt_p;
my $bams = $opt_b;
my $out = $opt_o;
my $reference = $opt_r;
my $threads = $opt_t;



# YOU MAY CONFIGURE THE SCRIPT WITH YOUR LOCAL INFO
=c
my $vcf = "vcffile.vcf";
my $ped = "";
my $bams = "/home/User/bams/";
my $out = "/home/User/out";
my $reference = "chr6.fasta";
my $threads = 10;
=cut


my $perl = `which perl`;
chomp $perl;
if (! -e $perl) {print "Need Perl to proceed!\n";}

my $bcf = `which bcftools`;
chomp $bcf;
if ($opt_x ne "") {$bcf = $opt_x;}
if (! -e $bcf) {print "Need bcftools to proceed!\n";}

my $whatshap = `which whatshap`;
chomp $whatshap;
if ($opt_w ne "") {$whatshap = $opt_w;}
if (! -e $whatshap) {print "Need WhatsHap to proceed!\n";}

my $format = "PS";

if ($threads eq "") {$threads = "4";}


our $message = "";
if (! -e $vcf) {$message .= "Please indicate a VCF\n"; help();}
if (! -e $reference) {$message .= "Please provide a reference fasta\n"; help();}
if ($vcf eq "") {$message .= "Please indicate a VCF\n"; help();}
if (($ped ne "") && (! -e $vcf)) {$message .= "This PED file is invalid!\n"; help();}
if ($bams eq "") {$message .= "Please provide a folder with all BAM files!\n"; help();}
if (! -e $bams) {$message .= "This BAM folder is invalid!\n"; help();}
if (! -e $out) {mkdir $out;}
if (! -e $out) {$message .= "Can't create the output structure!\n"; help();}

mkdir "$out/splited_vcfs";
mkdir "$out/phased_vcfs";






my %ped_data = ();
my %ped_line = ();
my %ped_samples = ();
if ($ped ne "") {
	print "Loading PED...\n";
	open (IN, $ped) or die;
	while (<IN>)
	{
		chomp $_;
		if ($_ eq "") {next;}
		if ($_ =~ /Maternal/) {next;}
		my @data = split("\t",$_);
		$ped_line{$data[1]} = $_;
		$ped_samples{$data[1]} = 1;
		$ped_data{$data[1]}{m} = $data[2];
		if ($data[2] ne 0) {$ped_samples{$data[2]} = 1;}
		$ped_data{$data[1]}{f} = $data[3];
		if ($data[3] ne 0) {$ped_samples{$data[3]} = 1;}
	}
	close (IN);
}








my %used = ();
my %db = ();
my @head = ();
my @samples = ();
my $index = 0;
my $index_detect = 0;
my @snps = ();

print "Loading VCF ... this might take a while ...\n";
open (IN, $vcf) or die;
while (<IN>)
{
	chomp $_;
	if ($_ eq "") {next;}
	if (substr($_,0,2) eq "##") {push(@head,$_); next;}
	if (substr($_,0,2) eq "#C") 
	{
		@samples = split("\t",$_);
		foreach (@samples) 
		{
			if ($_ eq "FORMAT") {$index_detect = 1; last;}
			$index++;
		}
		if ($index_detect eq 0) {$message .= "This VCF file does not have a FORMAT field\n"; help();}
		$index++;
		next;
	}

	my @data = split("\t",$_);
	if ($used{$data[1]} eq 1) {$message .= "There are multiple lines with the same position\n"; help();}
	
	for (my $a = $index; $a < scalar(@data); $a++)
	{
		my $value = $data[$a];
		if (exists $ped_samples{$samples[$a]}) {$db{$data[1]}{$samples[$a]} = $value; next;}
		my @fields = split(":",$value);
		if ($fields[0] =~ /\./) {next;}
		my @alleles = split("/",$fields[0]);
		if ($alleles[0] eq $alleles[1]) {next;}
		$db{$data[1]}{$samples[$a]} = $value;
	}
	
	my $var = $data[0];
	for (my $a = 1; $a < $index; $a++)
	{
		$var .= "\t$data[$a]";
	}
	push(@snps,$var);
}
close (IN);
undef %used;
#goto joining;




print "Loading BAMs ...\n";
my @bam = <$bams/*.bam>;
if (scalar(@bam) eq 0) {@bam = <$bams/*/*.bam>;}
if (scalar(@bam) eq 0) {@bam = <$bams/*/*/*.bam>;}
if (scalar(@bam) eq 0) {$message .= "There are no BAM files in this folder\n"; help();}

my @sample_list = ();
my %sample_bam = ();
for (my $a = $index; $a < scalar(@samples); $a++)
{
	push(@sample_list,$samples[$a]);
}

foreach (@sample_list)
{
	my $sample = $_;
	foreach (@bam)
	{
		if ($_ =~ /$sample/) {$sample_bam{$sample} = $_; last;}
	}
}

my $fail = "";
foreach (@sample_list)
{
	if (! exists $sample_bam{$_}) {$fail .= ",$_";}
}
if ($fail ne "") {$message .= "Can't find BAM for samples $fail\n\n"; exit;}








print "Generating temporary VCFs ... \n";
my @cmds = ();
my %used = ();

if ($ped ne "") {
	my $group = 1;
	my @ped_samples = sort keys %ped_data;
	foreach (@ped_samples)
	{
		my $ped_out = "$out/splited_vcfs/group_" . $group . ".ped";
		open (OUT, ">$ped_out");
		print OUT $ped_line{$_};
		close (OUT);
		
		my $vcf_out = "$out/splited_vcfs/group_" . $group . ".vcf";
		open (OUT, ">$vcf_out");
		
		foreach (@head)
		{
			print OUT $_ . "\n";
		}
		for (my $a = 0; $a < $index; $a++)
		{
			print OUT $samples[$a] . "\t";
		}
		
		print OUT $_;
		$used{$_} = 1;
		if ($ped_data{$_}{m} ne "0") {print OUT "\t$ped_data{$_}{m}";$used{$ped_data{$_}{m}} = 1;}
		if ($ped_data{$_}{f} ne "0") {print OUT "\t$ped_data{$_}{f}";$used{$ped_data{$_}{f}} = 1;}
		print OUT "\n";
		my $sample = $_;

		
		foreach (@snps)
		{
			my @data = split("\t",$_);
			print OUT $_;

			print OUT "\t" . $db{$data[1]}{$sample};
			
			if ($ped_data{$sample}{m} ne "0") 
			{
				print OUT "\t" . $db{$data[1]}{$ped_data{$sample}{m}};
			}
			if ($ped_data{$sample}{f} ne "0") 
			{
				print OUT "\t" . $db{$data[1]}{$ped_data{$sample}{f}};
			}
			print OUT "\n";
		}
		close (OUT);
		
		my $cmd = "$whatshap phase --indels --ped $out/splited_vcfs/group_$group\.ped -o $out/phased_vcfs/group_$group.vcf --reference $reference --tag $format $out/splited_vcfs/group_$group\.trim.clean.vcf";
		$cmd .= " " . $sample_bam{$sample};
		if ($ped_data{$sample}{m} ne "0") {$cmd .= " " . $sample_bam{$ped_data{$sample}{m}};}
		if ($ped_data{$sample}{f} ne "0") {$cmd .= " " . $sample_bam{$ped_data{$sample}{f}};}
		$cmd .= " 2> $out/phased_vcfs/group_$group.log";
		push(@cmds,$cmd);

		$group++;
	}
}





for ($a = $index; $a < scalar(@samples); $a++)
{
	my $sample = $samples[$a];
	if (exists $used{$sample}) {next;}
	
	my $sample_index = 0;
	my $count = 0;
	foreach (@samples) {if ($_ eq $sample) {$sample_index = $count; last;} $count++;}
	
	
	my $vcf_out = "$out/splited_vcfs/" . $sample . ".vcf";
	open (OUT, ">$vcf_out");
	
	foreach (@head)
	{
		print OUT $_ . "\n";
	}
	for (my $a = 0; $a < $index; $a++)
	{
		print OUT $samples[$a] . "\t";
	}
	
	print OUT $sample;
	print OUT "\n";
	foreach (@snps)
	{
		my @data = split("\t",$_);
		my $var = $data[0] . "\t" . $data[1];
		
		if (! exists $db{$data[1]}{$sample}) {next;}
	
		print OUT $_;
		print OUT "\t" . $db{$data[1]}{$sample};
		print OUT "\n";
	}
	close (OUT);
	
	
	my $cmd = "$whatshap phase --indels -o $out/phased_vcfs/$sample\.vcf --reference $reference --tag $format $out/splited_vcfs/" . $sample . ".trim.clean.vcf";
	
	$cmd .= " " . $sample_bam{$sample};
	$cmd .= " 2> $out/phased_vcfs/$sample.log";
	push(@cmds,$cmd);
}
undef %used;
%db = ();
undef %db;






print "Trimming temporary VCFs ... \n";
my $bcf = `which bcftools`;
chomp $bcf;

my @vcfs = <$out/splited_vcfs/*.vcf>;
foreach (@vcfs)
{
	if ($_ =~ /\.trim/){next;}
	my $base = $_;
	$base =~ s/\.vcf//;
	my $cmd = "$bcf view --trim-alt-alleles --min-ac 1 $_ > $base\.trim.vcf";
	system($cmd);
	open (IN, "$base\.trim.vcf");
	open (OUT, ">$base\.trim.clean.vcf");
	
	while (<IN>)
	{
		if (substr($_,0,1) eq "#") {print OUT $_; next;}
		my @data = split("\t",$_);
		if ($data[4] eq ".") {next;}
		if ($data[4] eq "*") {next;}
		print OUT $_;
	}
	close IN;
	close OUT;
}






print "Running WhatsHap with $threads parallel runs ...\n";
open (OUT,">$out/commands.txt");
foreach (@cmds)
{
	print OUT $_ . "\n";
}
close (OUT);




our $max_segments = scalar(@cmds);
#Running WhatsHap in parallel
my $sem = Thread::Semaphore->new($threads); # max num of threads
my $current_segment = 0;

my @threads = map {
	# request a thread slot, waiting if none are available:
	$current_segment++; 
	$sem->down;
	
	threads->create(\&run_cmd, $cmds[($current_segment-1)], $current_segment)
	
} 1..$max_segments;
$_->join for @threads;





joining:
print "Loading phased data ...\n";

my %phased = ();
my @subsamples = ();
my @files = <$out/phased_vcfs/*.vcf>;
my $count=0;
foreach (@files)
{
	my $file = $_;
	open (IN, $file);
	while (<IN>)
	{
		chomp $_;
		if ($_ eq "") {next;}
		if (substr($_,0,2) eq "#C") 
		{
			@subsamples = split("\t",$_);
			next;
		}
		
		my @data = split("\t",$_);
		
		my $ref = $data[3];
		my @alts = split(",",$data[4]);
		my @alleles = ();
		push(@alleles,$ref);
		foreach (@alts)
		{
			push(@alleles, $_);
		}
		
		my $pos = $data[1];
		
		for (my $a = $index; $a < scalar(@data); $a++)
		{
			my @fields = split(":",$data[$a]);
			if ($fields[0] =~ /\//) {next;}
			my @base = split("\\|",$fields[0]);
			my $gen = $alleles[$base[0]] . "|" . $alleles[$base[1]];
			for (my $b = 1; $b < scalar(@fields); $b++)
			{
				$gen .= ":$fields[$b]";
			}
			$phased{$subsamples[$a]}{$pos} = $gen;
		}
		
	}
	close (IN);
	$count++;
}






print "Writing phased VCF ...\n";
open (OUT,">$out/whatshap.vcf");
open (IN, $vcf) or die;

while (<IN>)
{
	chomp $_;
	if ($_ eq "") {next;}
	if (substr($_,0,1) eq "#") {print OUT $_ . "\n"; next;}

	my @data = split("\t",$_);
	print OUT $data[0];
	for (my $a = 1; $a < $index; $a++){
		print OUT "\t" . $data[$a];
	}
	print OUT ":$format";
	
	my $pos = $data[1];
	
	my %alleles = ();
	$alleles{$data[3]} = 0;
	my @alts = split(",",$data[4]);
	my $count = 1;
	foreach (@alts)
	{
		$alleles{$_} = $count;
		$count++;
	}
	
	for (my $a = $index; $a < scalar(@samples); $a++){
		my $sample = $samples[$a];
		my $value = $data[$a];
		if (exists $phased{$sample}{$pos}) 
		{
			my @data = split(":",$phased{$sample}{$pos});
			my @bases = split("\\|",$data[0]);
			my $gen = $alleles{$bases[0]} . "|" . $alleles{$bases[1]};
			for (my $b = 1; $b < scalar(@data); $b++)
			{
				$gen .= ":$data[$b]";
			}
			$value = $gen;
		}
		else {$value .= ":.";}
		print OUT "\t$value";
	}
	print OUT "\n";
}

print "All done!\n";




sub help {
	print "\n$message\n";
	
	print "-v vcf_file\n";
	print "-p ped_file [optional]\n";
	print "-b folder_with_the_bam_files\n";
	print "-o output_folder\n";
	print "-r chr6_reference_sequence\n";
	print "-t number_of_threads [default: 4]\n";
	print "-w path_to_whatshap [default: whatshap]\n";
	print "-x path_to_bcftools [default: bcftools]\n";
	
	print "\n\n";

	
	exit;
}





sub run_cmd {
	 	my ($cmd, $run) = @_;
		print "Starting run number $run of $max_segments ...\n";
		system ($cmd);
		$sem->up;
	}
