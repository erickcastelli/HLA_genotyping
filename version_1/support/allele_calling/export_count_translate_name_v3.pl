#!/usr/bin/perl

use strict;
use Getopt::Std;
use File::Basename;


our ($opt_v,$opt_o,$opt_t, $opt_b, $opt_c, $opt_p,$opt_r,$opt_s,$opt_e,$opt_g,$opt_i);

getopts('o:v:t:c:p:b:r:s:e:g:i');

#goto saida;

if (! -e $opt_v) {print "\nVCF file not informed or not detected\n\n";help();}
if (! -e $opt_r) {print "\nReference file not informed or not detected\n\n";help();}
if ($opt_o eq "") {print "\nPlease indicate the output folder -o\n\n";help();}
if ($opt_s eq "") {print "\nPlease indicate the genomic starting point -s\n\n";help();}
if ($opt_e eq "") {print "\nPlease indicate the genomic ending point -e\n\n";help();}
if (! -e $opt_b) {print "\nPlease indicate the BED file for CDS -b\n\n";help();}
if ($opt_t eq "") {print "\nPlease indicate a tag -t\n\n";help();}

if (! -e $opt_o) {system("mkdir $opt_o");}

my $vcfx = `which vcfx`;
chomp $vcfx;
if ($vcfx eq "") {"vcfx must be installed and available in your path\n";exit;}

my $transeq = `which transeq`;
chomp $transeq;
if ($transeq eq "") {"transeq from Emboss must be installed and available in your path\n";exit;}


##### GENOMIC

#goto saida2;

my $cmd = "";
if ($opt_i)
{
	$cmd = "$vcfx fasta start=$opt_s end=$opt_e output=$opt_o/$opt_t.genomic.tmp input=$opt_v reference=$opt_r";
}

if (! $opt_i)
{
	$cmd = "$vcfx fasta start=$opt_s end=$opt_e output=$opt_o/$opt_t.genomic.fas input=$opt_v reference=$opt_r";
}
print "\n$cmd\n\n";
system ($cmd);

if ($opt_i)
{
	open (OUT, ">$opt_o/$opt_t.genomic.fas");
	open (IN, "$opt_o/$opt_t.genomic.tmp");
	while (<IN>)
	{
		print OUT $_;
		my $seq = <IN>;
		chomp $seq;
		$seq = uc $seq;
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/;
		print OUT $seq . "\n";
	}
	close (IN);
	close (OUT);
	unlink "$opt_o/$opt_t.genomic.tmp";
}

saida2:
my %gen = ();
my %gen_ref = ();

if (-e $opt_g)
{
	print "\nLoading Genomic reference file...\n\n";

	my $id = "";
	open (IN, $opt_g);
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {my @tmp = split(" ",$_); $id = $tmp[1];next;}
		$gen{$id} .= uc $_;
	}
	close (IN);
	
	my @ids = keys %gen;
	foreach (@ids)
	{
		$gen_ref{$gen{$_}} .= ",$_";
	}
	undef %gen;
}

print "\nCounting Genomic sequences...\n\n";

open (IN, "$opt_o/$opt_t.genomic.fas");
my %gen = ();
while (<IN>)
{
	chomp $_;
	my $id = $_;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$gen{$seq}++;
}
close IN;

my @seqs = keys %gen;
my @refs = sort keys %gen_ref;
open (OUT, ">$opt_o/$opt_t.genomic.counted.fas");
my $unknown = 1;
foreach (@seqs)
{
	my $seq = $_;
	my $count = $gen{$seq};
	my $id = "unknown_$unknown";
	my $newid = "";

	foreach (@refs)
	{
		if (($_ =~ /$seq/) || ($seq =~ /$_/))
		{
			$newid .= "," . substr($gen_ref{$_},1);
			next;
		}
		else {
			my $subseq = $seq;
			my $subtest = $_;
			my $upseq = substr($seq,0,20);
			my $downseq = substr($seq,(length $seq)-20,20);
			my $uptest = substr($_,0,20);
			my $downtest = substr($_,(length $_)-20,20);
			
			if ($subtest =~ /$upseq/) 
			{
				my $result = index($subtest, $upseq);
				$subtest = substr($subtest,$result);
			}
			if ($subseq =~ /$uptest/) 
			{
				my $result = index($subseq, $uptest);
				$subseq = substr($subseq,$result);
			}
							
			if ($subseq =~ /$downtest/) 
			{
				my $result = index($subseq, $downtest);
				$subseq = substr($subseq,0,$result + length $downtest);
			}				
			if ($subtest =~ /$downseq/) 
			{
				my $result = index($subtest, $downseq);
				$subtest = substr($subtest,0,$result + length $downseq);
			}				

			if (($subtest =~ /$subseq/) || ($subseq =~ /$subtest/))
			{
				$newid .= "," . substr($gen_ref{$_},1);
			}
		}
	}
	if ($newid ne "") {$newid = substr($newid,1);}
	if ($newid eq "") {$newid = $id;}
	if ($newid eq "unknown_$unknown") {$unknown++;}
	
	my @tmp = split(",",$newid);
	@tmp = sort @tmp;
	$newid = join(",",@tmp);
	
	print OUT  ">$newid;$count\n$seq\n";
}


undef %gen;
undef %gen_ref;





### TRANSCRIPT
my $cmd = "";
if (! $opt_i) {
	$cmd = "$vcfx transcript bed=$opt_b output=$opt_o/$opt_t.cds.fas input=$opt_v reference=$opt_r";
}
if ($opt_i) {
	$cmd = "$vcfx transcript bed=$opt_b output=$opt_o/$opt_t.cds.tmp input=$opt_v reference=$opt_r";
}

print "\n$cmd\n\n";
system ($cmd);


if ($opt_i)
{
	open (OUT, ">$opt_o/$opt_t.cds.fas");
	open (IN, "$opt_o/$opt_t.cds.tmp");
	while (<IN>)
	{
		print OUT $_;
		my $seq = <IN>;
		chomp $seq;
		$seq = uc $seq;
		$seq = reverse $seq;
		$seq =~ tr/ATCG/TAGC/;
		print OUT $seq . "\n";
	}
	close (IN);
	close (OUT);
	unlink "$opt_o/$opt_t.cds.tmp";
}




my %cds = ();
my %cds_ref = ();

if (-e $opt_c)
{
	print "\nLoading CDS reference file...\n\n";

	my $id = "";
	open (IN, $opt_c);
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {my @tmp = split(" ",$_); $id = $tmp[1];next;}
		$cds{$id} .= uc $_;
	}
	close (IN);
	
	my @ids = keys %cds;
	foreach (@ids)
	{
		$cds_ref{$cds{$_}} .= ",$_";
	}
	undef %cds;
}

print "\nCounting CDS sequences...\n\n";

open (IN, "$opt_o/$opt_t.cds.fas");
my %cds = ();
while (<IN>)
{
	chomp $_;
	my $id = $_;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$cds{$seq}++;
}
close IN;

my @seqs = keys %cds;
my @refs = keys %cds_ref;
open (OUT, ">$opt_o/$opt_t.cds.counted.fas");
my $unknown = 1;
foreach (@seqs)
{
	my $seq = $_;
	my $count = $cds{$seq};
	my $id = "unknown_$unknown";
	if ($cds_ref{$seq} ne "") {$id = substr($cds_ref{$seq},1);}
	if ($cds_ref{$seq} eq "") 
	{
		foreach (@refs)
		{
			if (($_ =~ /$seq/) || ($seq =~ /$_/))
			{
				$id = substr($cds_ref{$_},1) . "_compatible";
				last;
			}
		}
		if ($id eq "unknown_$unknown") {$unknown++;}
	}
	print OUT  ">$id;$count\n$seq\n";
}




### PROTEIN

print "\nTranslating...\n\n";
my $cmd = "$transeq --sequence $opt_o/$opt_t.cds.fas -out $opt_o/$opt_t.cds.protein.tmp";
system($cmd);

open (IN, "$opt_o/$opt_t.cds.protein.tmp");
my $id = "";
my %prot = ();
while (<IN>)
{
	if (substr($_,0,1) eq ">") {$_ =~ s/_1\n//; $id = $_;next;}
	chomp $_;
	$prot{$id} .= uc $_;
}
close (IN);

my @ids = sort keys %prot;
open (OUT, ">$opt_o/$opt_t.cds.protein.fas");
foreach (@ids)
{
	my @data = split("\\*",$prot{$_});
	print OUT "$_\n$data[0]\n";
}
unlink "$opt_o/$opt_t.cds.protein.tmp";
undef %prot;
close OUT;



my %prot = ();
my %prot_ref = ();

if (-e $opt_p)
{
	print "\nLoading Protein reference file...\n\n";

	my $id = "";
	open (IN, $opt_p);
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,1) eq ">") {my @tmp = split(" ",$_); $id = $tmp[1];next;}
		$prot{$id} .= uc $_;
	}
	close (IN);
	
	my @ids = keys %prot;
	foreach (@ids)
	{
		$prot_ref{$prot{$_}} .= ",$_";
	}
	undef %prot;
}




print "\nCounting Protein sequences...\n\n";

open (IN, "$opt_o/$opt_t.cds.protein.fas");
my %prot = ();
while (<IN>)
{
	chomp $_;
	my $id = $_;
	my $seq = <IN>;
	chomp $seq;
	$seq = uc $seq;
	$prot{$seq}++;
}
close IN;

my @seqs = keys %prot;
my @refs = keys %prot_ref;
my %used_small_count = ();

open (OUT, ">$opt_o/$opt_t.cds.protein.counted.fas");
my $unknown = 1;
foreach (@seqs)
{
	my $seq = $_;
	my $count = $prot{$seq};
	my $id = "unknown_$unknown";
	if ($prot_ref{$seq} ne "") {$id = substr($prot_ref{$seq},1);}
	if ($prot_ref{$seq} eq "") 
	{
		foreach (@refs)
		{
			if (($seq =~ /$_/) && (length $seq > length $_))
			{
				$used_small_count{$prot_ref{$_}}++;
				$id = "larger_than_ref_" . $used_small_count{$prot_ref{$_}} . "," . substr($prot_ref{$_},1);
				last;
			}
		}
		if ($id eq "unknown_$unknown") {$unknown++;}
	}
	print OUT  ">$id;$count\n$seq\n";
}
close OUT;



saida:

### DB


print "\nGenerating DB ...\n\n";


undef %prot;
undef %prot_ref;
undef %cds;
undef %cds_ref;
undef %gen;
undef %gen_ref;



my %gen = ();
open (IN, "$opt_o/$opt_t.genomic.fas");
while (<IN>)
{
	my $id = $_;
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$gen{$id} = $seq;
}
close (IN);

my %genomic_ref = ();
open (IN, "$opt_o/$opt_t.genomic.counted.fas");
while (<IN>)
{
	my $id = $_;
	my @id = split(";",$id);
	$id = $id[0];
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$gen_ref{$seq} = $id;
}
close (IN);



my %cds = ();
open (IN, "$opt_o/$opt_t.cds.fas");
while (<IN>)
{
	my $id = $_;
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$cds{$id} = $seq;
}
close (IN);

my %cds_ref = ();
open (IN, "$opt_o/$opt_t.cds.counted.fas");
while (<IN>)
{
	my $id = $_;
	my @id = split(";",$id);
	$id = $id[0];
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$cds_ref{$seq} = $id;
}
close (IN);


my %prot = ();
open (IN, "$opt_o/$opt_t.cds.protein.fas");
while (<IN>)
{
	my $id = $_;
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$prot{$id} = $seq;
}
close (IN);

my %prot_ref = ();
open (IN, "$opt_o/$opt_t.cds.protein.counted.fas");
while (<IN>)
{
	my $id = $_;
	my @id = split(";",$id);
	$id = $id[0];
	my $seq = uc <IN>;
	chomp $id;
	chomp $seq;
	$prot_ref{$seq} = $id;
}
close (IN);



open (OUT, ">$opt_o/$opt_t.db.txt");
print OUT "Sample\tVector\tGenomic\tCDS\tCDS_clean\tProt\tProt_clean\n";
my @ids = sort keys %cds;
foreach (@ids)
{
	my @data = split("_",substr($_,1));
	print OUT $data[0] . "\t" . $data[1] . "\t";

	my $seq = $gen{$_};
	print OUT substr($gen_ref{$seq},1);
	print OUT "\t";

	my $seq = $cds{$_};
	print OUT substr($cds_ref{$seq},1);
	print OUT "\t";
	
	my $listall = substr($cds_ref{$seq},1);
	my @list = split(",",$listall);
	my %itens = ();
	foreach (@list)
	{
		my @fields = split("\\:",$_);
		my $new = $fields[0];
		if ($fields[1] ne "") {$new .= ":" . $fields[1];}
		if ($fields[2] ne "") {$new .= ":" . $fields[2];}
		$itens{$new} = 1;
	}
	my @res = sort keys %itens;
	print OUT join(",",@res);
	undef @res;
	undef %itens;
	undef @list;	
		
	my $seq = $prot{$_};
	print OUT "\t" . substr($prot_ref{$seq},1);
	print OUT "\t";
	
	
	my $listall = substr($prot_ref{$seq},1);
	my $new = $listall;
	
	if ($listall !~ /larger/) {
		my @list = split(",",$listall);
		my %itens = ();
		foreach (@list)
		{
			my @fields = split("\\:",$_);
			$new = $fields[0];
			if ($fields[1] ne "") {$new .= ":" . $fields[1];}
			$itens{$new} = 1;
		}
		my @res = sort keys %itens;
		print OUT join(",",@res);
	}
	else {
		print OUT $new;
	}

	undef @res;
	undef %itens;
	undef @list;	
	
	
	
	print OUT "\n";
}
close (OUT);

sub help
{
	print "-v [vcf file]\n";
	print "-o [output folder]\n";
	print "-t [tag, file prefix]\n";
	print "-c [cds fasta file]\n";
	print "-p [protein fasta file]\n";
	print "-g [genomic fasta file]\n";
	print "-b [bed file for CDS]\n";
	print "-s [start position for genomic]\n";
	print "-e [end position for genomic]\n";
	print "-r [chr reference in fasta format]\n";
	print "-i [Invert sequences]\n";
	
	exit;
}