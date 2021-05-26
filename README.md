# HLA Genotyping, haplotyping, and allele calling from next-generation sequencing data
Tutorial for genotyping, haplotyping, and allele calling for HLA genes

Version 1.0 (May 25th, 2021)

Author: Erick C. Castelli (erick.castelli@unesp.br)

The advantage of this method is that it calls variants, haplotypes, and then HLA alleles from the phased VCF. Thus, you have many different levels of information.

This tutorial is suitable for whole-genome sequencing (WGS), whole-exome sequencing (WES), and amplicon sequencing. It has been tested with Illumina and Ion sequencing data. Please note that coverage is essential. We recommend a coverage of at least 30x for WGS, and 100x for WES and amplicons. Read size is also important, and you will get better results when dealing with a read size larger than 75 nucleotides.

This procedure needs sample size. The minimum sample size we have tested in 150 samples. For a single-sample allele call, please use another method (such as HLA-LA). However, if you have a large sample size, this procedure is highly accurate and will provide you with SNPs (in the hg38 reference genome), haplotypes, and allele calls, and allow you to detect new variants.

## How to cite this pipeline
This pipeline is described in:

Hla-C genetic diversity and evolutionary insights in two samples from Brazil and Benin. HLA. 2020 Oct;96(4):468-486. doi: 10.1111/tan.13996

Immunogenetics of resistance to SARS-CoV-2 infection in discordant couples. medRxiv, doi: 10.1101/2021.04.21.21255872

## Packages and software needed
- hla-mapper 4 (www.castelli-lab.net/apps/hla-mapper) 
- GATK 4 (https://gatk.broadinstitute.org/hc/en-us)
- GATK 3.8 or WhatsHap (https://whatshap.readthedocs.io/en/latest/)
- vcfx 2 (www.castelli-lab.net/apps/vcfx)
- phasex 0.8.2 (https://github.com/erickcastelli/phasex)
- samtools 1.12 (http://samtools.sourceforge.net)
- BWA 0.7.17 (https://sourceforge.net/projects/bio-bwa/files/)
- bcftools 1.12 (http://samtools.github.io/bcftools/)
- vcftools (http://vcftools.sourceforge.net)
- IGV (https://software.broadinstitute.org/software/igv/)

## STEP 1: Using hla-mapper to get unbiased rad mapping for HLA genes
This step is essential. You won't retrieve correct genotypes in HLA genes unless using a method tailored for these genes. We recommend the use of hla-mapper.

hla-mapper supports many genes in the MHC region. Please check its website for instructions (www.castelli-lab.net/apps/hla-mapper)

There are two possible inputs for hla-mapper, a BAM file (step 1a), or FASTQ files (step 1b).
 
### STEP 1A (IF YOU WANT ALSO TO GENOTYPE INTERGENIC VARIANTS): 
- Download a copy of the human reference genome (hg38) and prepare it for BWA. We recommend the use of a reference genome without the alternative contigs.
- Prepare it for BWA:
> bwa index reference_genome
- Using BWA MEM, map your reads against the reference genome, and sort it:
> bwa mem reference_genome R1.fastq R2.fastq > sample.sam
> 
> samtools sort sample.sam > sample.bam
> 
> samtools index sample.bam
> 
> hla-mapper dna bam=sample.bam db=hla_mapper_database sample=Sample_Name output=output_folder
- You need to repeat this last part for each sample.
- You need to indicate a different output folder for each sample


### STEP 1B (IF YOU WILL GENOTYPE ONLY THE HLA GENES): 
> hla-mapper dna r1=R1.fastq.gz r2=R2.fastq.gz db=hla_mapper_database sample=Sample_Name output=output_folder
- You need to repeat this last part for each sample.
- You need to indicate a different output folder for each sample

## STEP 2 - Check some of the hla-mapper BAM files using IGV 
Using IGV, please check some of the hla-mapper outputted BAM files (Sample_Name.adjusted.bam) using IGV. Make sure everything is OK. You can compare the original BAM (using BWA MEM) with the new one (using hla-mapper).

## STEP 3 - Variant call using GATK 4
We recommend GATK 4 HaplotypeCaller to call variants. The hla-mapper BAM file is already prepared for GATK.

For each sample, run GATK HaplotypeCaller in the GVCF mode, such as this example:

> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -R reference_genome -I Sample_name.adjusted.bam -O output_folder/Sample_name.MHC.g.vcf -L chr6:29700000-33150000 -ERC GVCF --max-num-haplotypes-in-population 256 --native-pair-hmm-threads thread_number --allow-non-unique-kmers-in-ref TRUE

```diff
- Attention: If you are interested only in one gene (e.g., HLA-A), you should adjust the interval accordingly.
```

You need to adjust the amount of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the hla-mapper output BAM (Sample_name.adjusted.bam), the output folder (output_folder), the sample name, and the number of threads (thread_number).

After processing all your samples, you need to concatenate all G.VCF files in a single one. There are two ways to do that, depending on the number of samples. The most common one is using GATK CombineGVCFs (https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs). The other is GenomicsDBImport (https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport). This tutorial does not cover this issue. Please follow the GATK instructions and combine all GVCFS in a single file.

Now, you can genotype your GVCF using GATK GenotypeGVCFs, as follows:

> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar GenotypeGVCFs -R reference_genome -O output_folder/MHC.vcf -L chr6:29700000-33150000 --variant output_folder/All_samples.MHC.g.vcf --dbsnp path_to_dbsnp_vcf

You need to ajust the amont of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the single GVCF file (output_folder/All_samples.MHC.g.vcf), the output folder (output_folder), and the path to dbsnp (path_to_dbsnp_vcf). dbSNP is optional.

```diff
- Attention: If you are interested only in one gene (e.g., HLA-A), you should adjust the interval accordingly.
```

## STEP 4 - Variant refinement
There are many ways to proceed with variant refinement, i.e., removing artifacts, including GATK VQRS (very good for full genomes and exomes) and vcfx (better for small datasets). For HLA genes, we recommend vcfx.

Recode the VCF file using vcftools. This is important for the next steps to correct some minor encoding errors introduced by GATK.

> vcftools --vcf VCF_FILE --recode --out VCF_FILE_RECODE

Using any application you want, **or the script provided in /support/unphase_genotypes.pl**, please change any "|" allele separator for "/" in the recoded VCF file. 

Use vcfx to filter out artifacts and variants with too many missing alleles, as follows. Please check the vcfx manual to understand what is going on here (www.castelli-lab.net/apps/vcfx)

> vcfx checkpl input=VCF_FILE_RECODE (this will create a .pl.vcf file next to the original VCF)
> 
> bcftools view --trim-alt-alleles --min-ac 1 VCF_FILE_RECODE.pl.vcf > VCF_FILE_RECODE.pl.trim.vcf
> 
> vcfx evidence input=VCF_FILE_RECODE.pl.trim.vcf (this will create a .pl.trim.evid.vcf)
> 
> vcfx filter input=VCF_FILE_RECODE.pl.trim.evid.vcf tag=PASS,WARN (this will create a .pl.trim.evid.filter.vcf)
>
> vcftools --vcf VCF_FILE_RECODE.pl.trim.evid.filter.vcf --recode --out VCF_FILE_RECODE.pl.trim.evid.filter.recoded.vcf

The last VCF file contains only the variants that have passed the vcfx checkpl/evidence workflow. For now on, we will refer to this VCF file as "VCF".



```diff
- Attention: In this step, you should manually check your VCF file and remove possible artifacts that may have passed the vcfx workflow.
+ Please note that the VCF generated up to this step is suitable for association studies and other purposes, but it consists of unphased genotypes with some missing alleles.
```

***Attention: We recommend performing the next steps for each gene separately. For that, you need to extract from the VCF the variants overlapping the gene you are interested, and perform the next steps. You can use vcftools for that, indicating the intervals with --from-bp and --to-bp. You may also try this pipiline with all variants together, but this decreases accuracy for detecting HLA alleles***

## STEP 6 - Normalize your multi-allelic VCF to biallelic VCF
> bcftools norm -m-any VCF > BIALLELIC.VCF

## STEP 7 - Calling phasing sets directly from the sequencing data
In this step, we will infer phase sets (the micro haplotypes) directly from the sequencing data. For that, there are two options: ReadBackedPhasing grom GATK 3.8, or WhatsHap. Both methods work well, but here we will address only the ReadBackedPhasing method.

To use ReadBackedPhasing and parallelize runs in different cores, please use the support script from phasex (https://github.com/erickcastelli/phasex). You should get familiar with phasex to understang how to do it, but it is a simple task using a Perl script. The script will split your VCF files, one for each sample, run ReadBackedPhasing in parallel, and join all files in a single VCF. 

The input for this step is the BIALLELIC.VCF file.

After that, the VCF file produced by the script (BIALLELIC.RBP.VCF) contains phase sets in the HP format. These phase sets will be considered in the upcoming haplotyping procedure.


## STEP 8 - Removing unphased singletons
Shapeit4, which is used combined with phasex to call haplotypes, or any probabilistic method for haplotyping, does not handle properly unphased singletons. Singletons are variants that have occurred in just one sample in heterozygosis. Thus, we need to remove them to proceed to the next step.

**To do that, we can use script /support/remove_unphased_singleton_from_normalized_rbp.pl**

This method will remove from your BIALLELIC.RBP.VCF file all unphased singletons, producing 3 files: one VCF without singletons (BIALLELIC.RBP.NOSINGLETON.VCF), a separate file with the unphased singletons (BIALLELIC.RBP.unphased_singletons.VCF), and a log file. Please keep all these files.


## STEP 9 - Calling haplotypes
We will use phasex to call haplotypes. Please check https://github.com/erickcastelli/phasex for instructions in how to do it.

For each gene and after removing the unphased singletons, run:

> phasex hp-ps vcf=BIALLELIC.RBP.NOSINGLETON.VCF output=BIALLELIC.RBP.NOSINGLETON.PS.VCF

The command above will convert the HP format (from ReadBackedPhasing) to the PS format, which is compatible with phasex.

To proceed to the haplotyping step, run:

> phasex phase-ps vcf=BIALLELIC.RBP.NOSINGLETON.PS.VCF select=0.51

The command above will perform a phasex run with the default parameters, i.e., 50 iterations, 50 replicates, using half the number of cores in the machine, using a threshold to fix a haplotype of 95%, and selecting all samples in which the same haplotype pair is present in at least 51% of the final runs.

The final VCF file (in the phasex output folder, file results.vcf) is a phased VCF with the haplotypes for samples that present a haplotype with a empirical P value higher than 51%. You may adjust that according to your needs. Please check the phasex manual for details. 

```diff
+ Please note that the phased VCF generated up to this step is suitable for association studies and other purposes. For example, you can call HLA alleles, or you can convert it to a phased ARLEQUIN file (using vcfx arlequin).
```

## STEP 10 - Introducing unphased singletons (optional)

This step is optional. This procedure will reintroduce the unphased singletons, that they will be recorded as unphased genotypes.

These unphased genotypes are not considered when we export complete sequences and call HLA alleles in the next step.

**To reintroduce unphased singletons, use the script /support/insert_singletons_back_to_phased_data.pl**


## STEP 11 - Calling complete sequences and HLA alleles
For this step, we will generate complete sequences for each HLA gene, and compare this sequences with known ones from the IPD-IMGT/HLA database.

**You need to use the script /support/allele_calling/export_count_translate_name_v3.pl**

> perl export_count_translate_name_v3.pl

-v [indicate the phased VCF produced by phasex or the one with the unphased singletons re-introduced]

-o [indicate where the outputs should be placed]

-t [indicate a prefix for the files, such as the name the gene. Do not use spaces.]

-c [indicate the fasta file with known exonic sequences (nuc) from IPD-IMGT/HLA]

-p [indicate the fasta file with known protein sequences (prot) from IPD-IMGT/HLA]

-g [indicate the fasta file with known genomic sequences (gen) from IPD-IMGT/HLA]

-b [indicate the BED file for this gene (e.g., HLA-G.exons.bed)]

-s [indicate the starting position for this gene, which is indicated in the gemomic BED file (e.g., HLA-G.genomic.bed]

-e [end position for genomic]
-r [chr reference in fasta format]
-i [Invert sequences]
