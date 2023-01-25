#!/usr/bin/perl
# from the hla-mapper workflow, calling SNPs and alleles for HLA genes
# Erick C. Castelli, erick.castelli@unesp.br
# version 1.0b

use strict;
use threads;
use File::Basename;
use Getopt::Std;


my $sc = basename $0;
my $dir = dirname $0;

our ($opt_p,$opt_v,$opt_o,$opt_r,$opt_x,$opt_a);

getopts('p:v:o:r:x:');

my $profile = $opt_p;
my $vcf = $opt_v;
my $output = $opt_o . "/";
my $ref = $opt_r;
my $alt = $opt_a;


my $vcfx = `which vcfx`;
chomp $vcfx;
if ($opt_x ne "") {$vcfx = $opt_x;}


# YOUR MAY CONFIGURE THE SCRIPT WITH YOUR LOCAL VARIABLES
=c
my $profile = "HLA-A";
my $vcf = "final.vcf.gz";
my $output = "/home/lab/out/";
my $ref = "chr6.fasta";
my $alt = 0;
my $vcfx = "vcfx";
=cut




my $bed_cds = "";
my $bed_gen = "";
my $start = "";
my $end = "";
my $chr = "";
my $nuc = "";
my $gen = "";
my $prot = "";
my $invert = 0;
my $start = "";
my $end = "";
my $bed_exons = "";
our @warning = ();
my $master_start = 0;
my $master_end = 0;


if ($profile ne "")
{
	$nuc = "$dir/imgt_db/HLA/$profile\_nuc.fasta";
	$gen = "$dir/imgt_db/HLA/$profile\_gen.fasta";
	$prot = "$dir/imgt_db/HLA/$profile\_prot.fasta";
	$bed_cds = "$dir/bed/$profile\.CDS.bed";
	$bed_gen = "$dir/bed/$profile\.genomic.bed";
	$bed_exons = "$dir/bed/$profile\.exons.bed";
}


if (($profile ne "") && ($alt eq 1))
{
	$nuc = "$dir/imgt_db/HLA/$profile\_nuc.fasta";
	$gen = "$dir/imgt_db/HLA/$profile\_gen.fasta";
	$prot = "$dir/imgt_db/HLA/$profile\_prot.fasta";
	$bed_cds = "$dir/bed/$profile\.alt.CDS.bed";
	$bed_gen = "$dir/bed/$profile\.genomic.bed";
	$bed_exons = "$dir/bed/$profile\.alt.exons.bed";
}


if (($profile eq "HLA-B") || ($profile eq "HLA-C"))
{
	$invert = 1;
}
if (($profile eq "HLA-DRB1") || ($profile eq "HLA-DMA"))
{
	$invert = 1;
}
if (($profile eq "TAP2") || ($profile eq "TAP1"))
{
	$invert = 1;
}
if (($profile eq "HLA-DMB") || ($profile eq "HLA-DOA"))
{
	$invert = 1;
}
if (($profile eq "HLA-DQB1") || ($profile eq "HLA-DPA1"))
{
	$invert = 1;
}
if (($profile eq "HLA-DOB"))
{
	$invert = 1;
}





if (-e $bed_gen)
{
	open (IN, $bed_gen);
	my $line = <IN>;
	chomp $line;
	my @data = split("\t",$line);
	$start = $data[1];
	$end = $data[2];
	$chr = $data[0];
	close (IN);
}

if ($master_start > 0) {$start = $master_start;} 
if ($master_end> 0) {$end = $master_end;} 


my $valid = 1;
if (! -e $vcfx) {push(@warning,"Need vcfx to proceed!");$valid = 0;}
if (! -e $vcf) {push(@warning,"Invalid VCF file!");$valid = 0;}
if ($vcf !~ /.gz/) {push(@warning,"VCF must be compressed with BGZIP and indexed with tabix!");$valid = 0;}
if (! -e "$vcf.tbi") {push(@warning,"VCF must be compressed with BGZIP and indexed with tabix!");$valid = 0;}
if (! -e $nuc) {push(@warning,"Invalid NUC file!");$valid = 0;}
if (! -e $ref) {push(@warning,"Invalid reference file!");$valid = 0;}
if (! -e $gen) {push(@warning,"Invalid GEN file!");$valid = 0;}
if (! -e $prot) {push(@warning,"Invalid PROT file!");$valid = 0;}
if (! -e $bed_cds) {push(@warning,"Invalid CDS BED file!");$valid = 0;}
if (! -e $bed_exons) {push(@warning,"Invalid Exon BED file!");$valid = 0;}
if ($start <= 0) {push(@warning,"Invalid genomic BED or Start point!");$valid = 0;}
if ($end <= 0) {push(@warning,"Invalid genomic BED or End point!");$valid = 0;}
if ($valid eq 0) {help();}



mkdir ("$output");
if (! -e $output) {push(@warning,"Could not create output structure!");help();}

mkdir ("$output/$profile");
if (! -e "$output/$profile") {push(@warning,"Could not create output structure!");help();}
$output = $output . "/" . $profile;


print "\n";
print $sc . "\n";
print "Extracting $profile from VCF ...\n";
my $cmd = "bcftools view --threads 4 '$vcf' $chr:$start" . "-" . "$end" . " > '$output/$profile.vcf' 2>'$output/null.txt'";
system($cmd);
unlink "$output/null.txt";
if (! -e "$output/$profile.vcf") {print "Something went wrong!\n\n"; exit;}
if (-s "$output/$profile.vcf" eq 0) {print "Something went wrong!\n\n"; exit;}


print "Running vcfx...\n";
my ($thr1) = threads->create(sub {
				my $cmd = "$vcfx fasta input=$output/$profile.vcf start=$start end=$end reference=$ref output=$output/$profile\.genomic.fas --quiet";
				my $out = `$cmd`;
				});
my ($thr2) = threads->create(sub {
				my $cmd = "$vcfx transcript input=$output/$profile.vcf reference=$ref bed=$bed_cds output=$output/$profile\.cds.fas --quiet";
				my $out = `$cmd`;
				});
my ($thr3) = threads->create(sub {
				my $cmd = "$vcfx transcript input=$output/$profile.vcf reference=$ref bed=$bed_exons output=$output/$profile\.exon.fas --quiet";
				my $out = `$cmd`;
				});
	
$thr1->join();
$thr2->join();
$thr3->join();



if ($invert eq 1)
{
	my $cmd = "revseq -sequence $output/$profile\.genomic.fas -outseq $output/$profile\.genomic.tmp.fas";
	my $out = `$cmd`;
	my $cmd = "revseq -sequence $output/$profile\.cds.fas -outseq $output/$profile\.cds.tmp.fas";
	my $out = `$cmd`;
	my $cmd = "revseq -sequence $output/$profile\.exon.fas -outseq $output/$profile\.exon.tmp.fas";
	my $out = `$cmd`;
	
	unlink "$output/$profile\.genomic.fas";
	unlink "$output/$profile\.cds.fas";
	
	my %tmp = ();
	open (IN, "$output/$profile\.genomic.tmp.fas");
	my $id = "";
	while (<IN>)
	{
		chomp $_;
		if ($_ eq ""){next;}
		if (substr($_,0,1) eq ">") {my @data = split(" ",substr($_,1)); $id = $data[0]; next;}
		$tmp{$id} .= uc $_;
	}
	close (IN);
	open (OUT, ">$output/$profile\.genomic.fas");
	my @ids = sort keys %tmp;
	foreach (@ids)
	{
		print OUT ">$_\n$tmp{$_}\n";
	}
	close (OUT);
	undef %tmp;
	
	my %tmp = ();
	open (IN, "$output/$profile\.cds.tmp.fas");
	my $id = "";
	while (<IN>)
	{
		chomp $_;
		if ($_ eq ""){next;}
		if (substr($_,0,1) eq ">") {my @data = split(" ",substr($_,1)); $id = $data[0]; next;}
		$tmp{$id} .= uc $_;
	}
	close (IN);
	open (OUT, ">$output/$profile\.cds.fas");
	my @ids = sort keys %tmp;
	foreach (@ids)
	{
		print OUT ">$_\n$tmp{$_}\n";
	}
	close (OUT);
	undef %tmp;
	
	
	
	my %tmp = ();
	open (IN, "$output/$profile\.exon.tmp.fas");
	my $id = "";
	while (<IN>)
	{
		chomp $_;
		if ($_ eq ""){next;}
		if (substr($_,0,1) eq ">") {my @data = split(" ",substr($_,1)); $id = $data[0]; next;}
		$tmp{$id} .= uc $_;
	}
	close (IN);
	open (OUT, ">$output/$profile\.exon.fas");
	my @ids = sort keys %tmp;
	foreach (@ids)
	{
		print OUT ">$_\n$tmp{$_}\n";
	}
	close (OUT);
	undef %tmp;
	
	unlink "$output/$profile\.exon.tmp.fas";
	unlink "$output/$profile\.cds.tmp.fas";
	unlink "$output/$profile\.genomic.tmp.fas";
	
	
	
}






print "Translating CDS...\n";
my $cmd = "transeq -sequence $output/$profile\.exon.fas -outseq $output/$profile\.prot.tmp";
my $out = `$cmd`;
open (IN, "$output/$profile\.prot.tmp");
my $id = "";
my %tmp = ();
my @alleles = ();
while (<IN>)
{
	chomp $_;
	if (substr($_,0,1) eq ">")
	{
		my @data = split("_",$_);
		$id = $data[0] . "_" . $data[1];
		push(@alleles,$id);
		next;
	}
	$tmp{$id} .= $_; 
}
open (OUT, ">$output/$profile\.prot.fas");
foreach (@alleles)
{
	my @seq = split("\\*",$tmp{$_});
	my $seqnew = $seq[0];
	if ($seqnew =~ /N/) {$seqnew .= "X";}
	print OUT "$_\n$seqnew\n";
}
close (OUT);
undef %tmp;
unlink "$output/$profile\.prot.tmp";





print "Counting sequences and naming ...\n";

my $thr1 = threads->create("name_cds", $nuc);
my $thr2 = threads->create("name_prot", $prot);
my $thr3 = threads->create("name_gen", $gen);
$thr1->join();
$thr2->join();
$thr3->join();







print "Creating db ...\n";

my %gen = ();
open (IN, "$output/$profile\.genomic.counted.fas");
open (OUT, ">$output/$profile\.genomic.groups.txt");
while (<IN>)
{
	chomp $_;
	my @ids = split(";",substr($_,1));
	my $size = $ids[(scalar(@ids) - 2)];
	my @list = ();
	for (my $a = 0; $a < (scalar(@ids) - 2); $a++)
	{
		push(@list,$ids[$a]);
	}

	@list = sort @list;
	if (scalar(@list) > 1) {
		my $line = join(";",@list);
		print OUT $list[0] . "\t" . $line . "\n";
	}
	
	my $id = $list[0] . ";" . $size;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$gen{$seq} = $id;
}
close (IN);
close (OUT);

my %cds = ();
open (IN, "$output/$profile\.cds.counted.fas");
while (<IN>)
{
	chomp $_;
	my @ids = split(";",substr($_,1));
	my $size = $ids[(scalar(@ids) - 2)];
	my $id = $ids[0] . ";" . $size;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$cds{$seq} = $id;
}
close (IN);


my %prot = ();
open (IN, "$output/$profile\.prot.counted.fas");
while (<IN>)
{
	chomp $_;
	my @id = split(";",substr($_,1));
	my $size = $id[(scalar(@id) - 2)];
	$id = $id[0] . ";" . $size;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$prot{$seq} = $id;
}
close (IN);




my %db = ();
open (IN, "$output/$profile\.genomic.fas");
while (<IN>)
{
	chomp $_;
	my $id = substr($_,1);
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$db{$id}{gen} = $gen{$seq};
}
close (IN);


open (IN, "$output/$profile\.cds.fas");
while (<IN>)
{
	chomp $_;
	my $id = substr($_,1);
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$db{$id}{cds} = $cds{$seq};
}
close (IN);

open (IN, "$output/$profile\.prot.fas");
while (<IN>)
{
	chomp $_;
	my $id = substr($_,1);
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$seq =~ s/X//;
	$db{$id}{prot} = $prot{$seq};
}
close (IN);


my @ids = sort keys %db;
open (OUT, ">$output/$profile\.db.txt");
print OUT "Sample\tVector\tGenomic\tCDS\tProtein\n";
foreach (@ids)
{
	my @data = split("_",$_);
	print OUT $data[0] . "\t" . $data[1];
	print OUT "\t$db{$_}{gen}";
	print OUT "\t$db{$_}{cds}";
	print OUT "\t$db{$_}{prot}";
	print OUT "\n";
}







sub name_prot
{
	my $new = 1;
	my %used = ();
	my $prot = $_[0];
	my %prot_named = ();


	my %refs_tmp = ();
	my $id = "";
	open (IN, $prot) or die;
	while (<IN>)
	{
		chomp $_;
		if ($_ eq "") {next;}
		if (substr($_,0,1) eq ">")
		{
			my @data = split(" ",$_);
			$id = $data[1];
			next;
		}
		$refs_tmp{$id} .= $_;
	}
	close (IN);
	
	
	my %refs_tmp_2 = ();
	my @alleles = keys %refs_tmp;
	foreach (@alleles)
	{
		my $seq = uc $refs_tmp{$_};
		$refs_tmp_2{$seq} .= ";$_";
	}
	
	
	
	my @seqs = keys %refs_tmp_2;
	my @seqs_all = keys %refs_tmp_2;
	my @valid = ();
	foreach (@seqs)
	{
		my $query = $_;
		my $ignore = 0;
		foreach (@seqs_all)
		{
			my $subj = $_;
			if ($query eq $subj) {next;}
			if ($subj =~ /$query/) {$ignore = 1;}
		}
		if ($ignore eq 0)
		{
			push(@valid,$query);
		}
	}
	

	
	
	my %refs = ();
	my %trunc = ();

	foreach (@valid)
	{
		my $seq = uc $_;
		my $id = $refs_tmp_2{$seq};
		if ($seq =~ /X/) {$trunc{$id} = 2; } else {$trunc{$id} = 1;}
		$seq =~ s/X//;
		$refs{$seq} = "$id";
	}
	undef %refs_tmp;
	undef %refs_tmp_2;
	undef @valid;
	
		
	
	my @refs_list = keys %refs;
	
	my %sequences = ();
	open (IN, "$output/$profile\.prot.fas");
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {next;}
		if ($_ eq "") {next;}
		my $seq = $_ . "\n";
		$seq =~ s/X\n//;
		$seq =~ s/\n//;
		$sequences{uc $seq}++;
	}
	close (IN);
	
	my @sequence_list = keys %sequences;
	
	my %used = ();
	foreach(@sequence_list)
	{
		my $seq = $_;
		my $id = "";
		$id = $refs{$seq};
				
		
		if ($id eq "") 
		{

			foreach (@refs_list)
			{
				my $ref = $_;
				my $ref_id = $refs{$ref};
				if ($ref eq "") {next;}
				
			
				if ($trunc{$ref_id} eq 2) {
					if (($seq =~ /$ref/) && (substr($seq,length($seq)-10, 10) eq substr($ref,length($ref)-10, 10))) 
					{
						$id .= ";$ref_id.L";
					}
					next;
				}

			
				if ($trunc{$ref_id} eq 1) {
					if ($seq =~ /$ref/) 
					{
						my @names = split(";",$refs{$ref});
						my $change = "";
						foreach (@names)
						{
							if ($_ eq "") {next;}
							if ($trunc{$_} eq 1) {next;}
							$change .= ";$_" . ".L";
						}
						$id .= $change;
					}
					if ($ref =~ /$seq/) 
					{
						my @names = split(";",$refs{$ref});
						
						my $change = "";
						foreach (@names)
						{
							if ($_ eq "") {next;}
							if ($trunc{$_} eq 1) {next;}
							$change .= ";$_" . ".S";
						}
						$id .= $change;
					}
				}

				
			}
	
		}
	
		if ($id eq "") 
		{
			$id = ";newprot_$new";
			$new++;
		}
		
		if ($id !~ /newprot/) {
			my @names = split(";",$id);
			my @alpha = sort @names;
			my @full = ();
			my @dif = ();
			foreach (@alpha) {if ($_ =~ /.L/) {next;} if ($_ =~ /.S/) {next;} push(@full,$_);} 
			foreach (@alpha) {if ($_ =~ /.L/) {push(@dif,$_);}}
			foreach (@alpha) {if ($_ =~ /.S/) {push(@dif,$_);}}
			
			my $new = "";
			foreach (@full) {if ($_ eq ""){next;} $new .= ";$_";}
			foreach (@dif) {if ($_ eq ""){next;} $new .= ";$_";}
			$id = $new;
			
			my @first = split(";",$id);
			my $first = $first[1];
			my @fields = split(":",$first);
			my $prot_name = $fields[0];
			if ($fields[1] ne "") {$prot_name .= ":$fields[1]";}
			if (($first =~ /.L/) && ($prot_name !~ /.L/)) {$prot_name .= ".L";}
			if (($first =~ /.S/) && ($prot_name !~ /.S/)) {$prot_name .= ".S";}
			
			my $new = "";
			$used{$prot_name}++;
			$new = ";$prot_name";
			if ($used{$prot_name} >= 2) {$new .= ".$used{$prot_name}";}
			$new .= $id;
			$id = $new;
			
		}
		$prot_named{$seq} = substr($id,1);
	}
	
	open (OUT, ">$output/$profile\.prot.counted.fas");
	foreach(@sequence_list)
	{
		my $count = $sequences{$_};
		print OUT ">" . $prot_named{$_} . ";" .length($_) .  "aa;$count\n";
		print OUT $_ . "\n";
	}
	close (OUT);
}


sub name_cds
{
	my $new = 1;
	my %used = ();
	my $cds = $_[0];
	my %cds_named = ();
	
	
	my %refs_tmp = ();
	my $id = "";
	open (IN, $cds) or die;
	while (<IN>)
	{
		chomp $_;
		if ($_ eq "") {next;}
		if (substr($_,0,1) eq ">")
		{
			my @data = split(" ",$_);
			$id = $data[1];
			next;
		}
		$refs_tmp{$id} .= $_;
	}
	close (IN);
	
	
	my %refs_tmp_2 = ();
	my @alleles = keys %refs_tmp;
	foreach (@alleles)
	{
		my $seq = uc $refs_tmp{$_};
		$refs_tmp_2{$seq} .= ";$_";
	}
	
	
	
	my @seqs = keys %refs_tmp_2;
	my @seqs_all = keys %refs_tmp_2;
	my @valid = ();
	foreach (@seqs)
	{
		my $query = $_;
		my $ignore = 0;
		foreach (@seqs_all)
		{
			my $subj = $_;
			if ($query eq $subj) {next;}
			if ($subj =~ /$query/) {$ignore = 1;}
		}
		if ($ignore eq 0)
		{
			push(@valid,$query);
		}
	}
	
	
	
	
	my %refs = ();
	foreach (@valid)
	{
		my $seq = uc $_;
		my $id = $refs_tmp_2{$seq};
		$refs{$seq} = "$id";
	}
	undef %refs_tmp;
	undef %refs_tmp_2;
	undef @valid;
	
	
	
	my @refs_list = keys %refs;
	
	my %sequences = ();
	open (IN, "$output/$profile\.cds.fas");
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {next;}
		if ($_ eq "") {next;}
		$sequences{uc $_}++;
	}
	close (IN);
	
	my @sequence_list = keys %sequences;
	
	my %used = ();
	foreach(@sequence_list)
	{
		my $seq = $_;
		my $id = "";
		$id = $refs{$seq};
		
		
		if ($id eq "") 
		{
			
			foreach (@refs_list)
			{
				my $ref = $_;
				my $ref_id = $refs{$ref};
				if ($ref eq "") {next;}
				
				
				if ($seq =~ /$ref/) 
				{
					my @names = split(";",$refs{$ref});
					my $change = "";
					foreach (@names)
					{
						if ($_ eq "") {next;}
						$change .= ";$_" . ".L";
					}
					$id .= $change;
				}
				if ($ref =~ /$seq/) 
				{
					my @names = split(";",$refs{$ref});
					
					my $change = "";
					foreach (@names)
					{
						if ($_ eq "") {next;}
						$change .= ";$_" . ".S";
					}
					$id .= $change;
				}
			}
			
		}
		
		if ($id eq "") 
		{
			$id = ";newcds_$new";
			$new++;
		}
		
		if ($id !~ /newcds/) {
			my @names = split(";",$id);
			my @alpha = sort @names;
			my @full = ();
			my @dif = ();
			foreach (@alpha) {if ($_ =~ /.L/) {next;} if ($_ =~ /.S/) {next;} push(@full,$_);} 
			foreach (@alpha) {if ($_ =~ /.L/) {push(@dif,$_);}}
			foreach (@alpha) {if ($_ =~ /.S/) {push(@dif,$_);}}
			
			my $new = "";
			foreach (@full) {if ($_ eq ""){next;} $new .= ";$_";}
			foreach (@dif) {if ($_ eq ""){next;} $new .= ";$_";}
			$id = $new;
			
			my @first = split(";",$id);
			my $first = $first[1];
			my @fields = split(":",$first);
			my $cds_name = $fields[0];
			if ($fields[1] ne "") {$cds_name .= ":$fields[1]";}
			if ($fields[2] ne "") {$cds_name .= ":$fields[2]";}
			if (($first =~ /.L/) && ($cds_name !~ /.L/)) {$cds_name .= ".L";}
			if (($first =~ /.S/) && ($cds_name !~ /.S/)) {$cds_name .= ".S";}
			
			my $new = "";
			$used{$cds_name}++;
			$new = ";$cds_name";
			if ($used{$cds_name} >= 2) {$new .= ".$used{$cds_name}";}
			$new .= $id;
			$id = $new;
			
		}
		$cds_named{$seq} = substr($id,1);
	}
	
	open (OUT, ">$output/$profile\.cds.counted.fas");
	foreach(@sequence_list)
	{
		my $count = $sequences{$_};
		print OUT ">" . $cds_named{$_} . ";" .length($_) .  "bp;$count\n";
		print OUT $_ . "\n";
	}
	close (OUT);
}









sub name_gen
{
	my $new = 1;
	my %used = ();
	my $gen = $_[0];
	my %refs_tmp = ();
	my %gen_named = ();
	open (IN, $gen) or die;
	my $id = "";
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">")
		{
			my @data = split(" ",$_);
			$id = $data[1];
			next;
		}
		$refs_tmp{$id} .= $_;
	}
	close (IN);
	my %refs = ();
	my @alleles = keys %refs_tmp;
	foreach (@alleles)
	{
		my $seq = uc $refs_tmp{$_};
		$refs{$seq} .= ";$_";
	}
	undef %refs_tmp;
	my @refs_list = keys %refs;
	
	
	my %sequences = ();
	open (IN, "$output/$profile\.genomic.fas");
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {next;}
		if ($_ eq "") {next;}
		$sequences{uc $_}++;
	}
	close (IN);
	
	my @sequence_list = keys %sequences;
	
	foreach(@sequence_list)
	{
		my $seq = $_;
		my $id = "";
		$id = $refs{$seq};
		
		if ($id eq "") 
		{
			foreach (@refs_list)
			{
				my $ref = $_;
				if ($seq =~ /$ref/) 
				{
					$id .= $refs{$ref} . ".L";
				}
				if ($ref =~ /$seq/) 
				{
					$id .= $refs{$ref} . ".S";
				}				
			}
		}
		
		if ($id eq "") 
		{
			$id = ";newgenomic_$new";
			$new++;
		}
		
		$used{$id}++;

		if ($used{$id} >= 2)
		{
			$id .= "." . $used{$id};
		}
		$gen_named{$seq} = substr($id,1);
	}
	
	open (OUT, ">$output/$profile\.genomic.counted.fas");
	foreach(@sequence_list)
	{
		my $count = $sequences{$_};
		my $size = length($_);
		print OUT ">" . $gen_named{$_} . ";$size" . "bp;$count\n";
		print OUT $_ . "\n";
	}
	close (OUT);
}
print "All done!\n";
exit;


sub help 
{

	print "Warnings:\n";
	foreach (@warning)
	{
		print $_ . "\n";
	}
	print "\n";
	print "\n";
	print "$sc\n";
	print "-p profile (e.g., HLA-A, HLA-B)\n";
	print "-v vcf_file (compressed with bgzip and indexed with tabix)\n";
	print "-o output_folder\n";
	print "-r reference_fasta (e.g.,chr6.fasta)\n";
	print "-x vcfx_path (default, vcfx)\n";
	print "-a alternative_profile\n";
	print "\n";
	exit;
}

