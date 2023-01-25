
# HLA genotyping, haplotyping, and allele calling from short-read next-generation sequencing data
Version 2.0 (Jan 25th, 2023)

Author: Erick C. Castelli (erick.castelli@unesp.br)


## Important notes:
This pipeline was designed to call SNPs and indels in genes from the MHC region, get the haplotypes, and call HLA alleles directly from the phased VCF data.

Data compatibility: this tutorial is compatible with whole-genome sequencing (WGS), whole-exome sequencing (WES), and amplicon sequencing. It was tested with short reads from Illumina (WGS and WES). It might work with Ion with some adjustments.

System compatibility: macOS and Linux. We have tested with MacOS 10.15 and Ubuntu 18.04. Other versions might be compatible.

Read depth: please note that read depth is essential. We recommend coverage of at least 30x for WGS and 50x for WES and amplicons. 

Read size: you will get better results when dealing with a read size larger than 75 nucleotides and paired-end sequencing, although the pipeline is also compatible with single-end sequencing data.

Sample size: The minimum sample size we have tested is 150 samples. You can proceed with a single-sample analysis up to step 3, but not further. Please use another method (such as HLA-LA) for a single-sample allele call. 


## Dependences
The following list contains all the software used in this pipeline and their indicated versions. Newer or older versions might also work, but we haven't tested them.

Attention: Please use bcftools 1.13. The pipeline will not work in newer versions. 

hla-mapper 4 (www.castelli-lab.net/apps/hla-mapper)

GATK, 4.2.0 or higher (https://gatk.broadinstitute.org/hc/en-us)

WhatsHap, 1.4 or higher (https://whatshap.readthedocs.io/en/latest/)

vcfx 2 (www.castelli-lab.net/apps/vcfx)

shapeit4 (https://odelaneau.github.io/shapeit4/)

samtools, 1.16 or higher (http://samtools.sourceforge.net)

BWA, 0.7.17 (https://sourceforge.net/projects/bio-bwa/files/)

bcftools, 1.13 (http://samtools.github.io/bcftools/)

IGV, any version (https://software.broadinstitute.org/software/igv/)

vcftools, any version (https://vcftools.sourceforge.net/) 

Emboss, 6.6 (https://emboss.sourceforge.net/download/) 

BGZIP and TABIX


## How to cite this pipeline:

You should cite hla-mapper:

Hla-mapper: an application to optimize the mapping of hla sequences produced by massively parallel sequencing procedures. Human Immunology 2018. doi: 10.1016/j.humimm.2018.06.010

You should cite all the applications described in previous section (samtools, bcftools, shapeit4, etc)

And this pipeline was used in these two studies:

MHC Variants Associated With Symptomatic Versus Asymptomatic SARS-CoV-2 Infection in Highly Exposed Individuals. Front. Immunol., 28 September 2021, https://doi.org/10.3389/fimmu.2021.742881

MUC22, HLA-A, and HLA-DOB variants and COVID-19 in resilient super-agers from Brazil. Front. Immunol., 25 October 2022, https://doi.org/10.3389/fimmu.2022.975918


## STEP 1: Using hla-mapper to get unbiased read alignment for HLA genes
This step is essential. You won't retrieve accurate genotypes in HLA genes unless you use an alignment method tailored for these genes.
 
hla-mapper supports many genes in the MHC region. Please check its website for instructions (www.castelli-lab.net/apps/hla-mapper)


There are two possible inputs for hla-mapper, a BAM file (step 1A) or FASTQ files (step 1B). Step 1A allows you to get genotypes in intergenic regions, while step 1B focuses only on the HLA genes. Step 1A is also suitable if you already have a BAM file with reads aligned to the hg38 reference genome.


### STEP 1A: 
If you already have a BAM file aligned to the human reference genome hg38, you may skip the first part and go directly to the hla-mapper part. 
Download a copy of the human reference genome (hg38) and prepare it for BWA. We recommend using a reference genome with the "chr" annotation for chromosomes, such as hg38DH used in the 1000Genomes project. All the scripts in this pipeline were designed for a reference genome with "chr6".


Prepare it for BWA:

> bwa index reference_genome


Using BWA MEM, map your reads against the reference genome, and sort it. You should include the sample name with the "-R" option, such as in the following example:

> bwa mem -R '@RG\tID:foo\tSM:foo' reference_genome R1.fastq R2.fastq | samtools sort - > sample.bam
> samtools index sample.bam

Get unbiased alignments with hla-mapper 4. The sample name should be the same as the one used in the previous BWA step, such as:

> hla-mapper dna threads=number_of_threads bam=sample.bam db=hla_mapper_database sample=foo output=output_folder

Repeat this last part for each sample. Please indicate a different output folder for each sample.


## STEP 1B:
Get unbiased alignments with hla-mapper 4:
> hla-mapper dna threads=number_of_threads r1=R1.fastq.gz r2=R2.fastq.gz db=hla_mapper_database sample=Sample_Name output=output_folder

Repeat this last part for each sample. Please indicate a different output folder for each sample.

Please note that hla-mapper can handle uncompressed and compressed FASTQ files.
```diff
- Attention: hla-mapper supports single-end sequencing data. Instead of r1= and r2=, use r0= to indicate your single-end fastq.
```

## STEP 2 - Check some of the BAM files using IGV
Using IGV, please check some of the hla-mapper BAM files produced by hla-mapper (Sample_Name.adjusted.bam or Sample_Name.adjusted.nodup.bam). 

Make sure everything is OK. For step 1A, you can compare the original BAM (before hla-mapper optimization) with the new ones.


## STEP 3 - Variant call using GATK 4
We recommend GATK 4 HaplotypeCaller to call variants when using Illumina data. 

Freebayes also works very well, but this pipeline is focused on GATK. 

The hla-mapper BAM file is already prepared for GATK and freebayes. For Ion data, freebayes should produce better results. 

For each sample, run GATK HaplotypeCaller in the GVCF mode, such as this example:
> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -R reference_genome -I Sample_name.adjusted.bam -O output_folder/Sample_name.MHC.g.vcf -L chr6:29700000-33150000 -ERC GVCF --max-num-haplotypes-in-population 256 --native-pair-hmm-threads thread_number --allow-non-unique-kmers-in-ref TRUE

```diff
- Attention: If you are interested only in one gene (e.g., HLA-A), you can adjust the interval accordingly. You need to adjust the amount of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the hla-mapper output BAM (Sample_name.adjusted.bam), the output folder (output_folder), the sample name, and the number of threads (thread_number).
```

Repeat this step for every sample.

After processing all your samples, you need to combine all G.VCF files into one. There are two ways to do that, depending on the number of samples. The most common one is using GATK CombineGVCFs (https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs). The other is GenomicsDBImport (https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport). This tutorial does not cover this issue. Please follow the GATK instructions and combine all GVCFS into one.

Now, you can genotype your GVCF using GATK GenotypeGVCFs. If you used CombineGVFs, one example is this:
> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar GenotypeGVCFs -R reference_genome -O output_folder/MHC.vcf -L chr6:29700000-33150000 --variant output_folder/All_samples.MHC.g.vcf --dbsnp path_to_dbsnp_vcf
```diff
- Attention: If you are interested only in one gene (e.g., HLA-A), you can adjust the interval accordingly. You need to adjust the amount of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the hla-mapper output BAM (Sample_name.adjusted.bam), the output folder (output_folder), the sample name, and the number of threads (thread_number). dbSNP is optional.
```

## STEP 4 - Variant refinement
There are many ways to proceed with variant refinement, i.e., removing artifacts. Here, we will combine GATK VQSR (better for large sample sizes, WGS and WES) and vcfx (better for small datasets).

Recode the VCF file using vcftools. This is important for the following steps to correct some minor encoding errors sometimes introduced by GATK.
> vcftools --vcf VCF_FILE --recode --out VCF_FILE_RECODE

Use sed to change any "|" allele separator for "/" in VCF file.
> sed 's/|/\//g' VCF_FILE_RECODE.VCF > VCF_FILE_RECODE_TREATED

Compress the file with BGZIP and index it with TABIX.
> bgzip VCF_FILE_RECODE_TREATED 
> tabix -p vcf VCF_FILE_RECODE_TREATED_GZ


Use GATK VQSR to filter out artifacts. Please follow the GATK 4 instructions. 
An example of this step is as follows. The MHC.select.vcf.gz file is provided in the supplementary files.

> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar VariantRecalibrator -R hg38.fasta -V VCF_FILE_RECODE_TREATED_GZ -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode BOTH -O vcf.recal --tranches-file vcf.tranches --resource:local,known=false,training=true,truth=true,prior=15.0 MHC.select.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 All_20180418.chr6.vcf.gz --resource:mills,known=false,training=true,truth=false,prior=12.0 /resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar ApplyVQSR -R hg38.fasta -V VCF_FILE_RECODE_TREATED_GZ -O VCF_FILE_RECODE_TREATED_VQSR --tranches-file vcf.tranches --recal-file vcf.recal --mode BOTH --truth-sensitivity-filter-level 99.0

```diff
- Attention: You need to adjust the amount of memory (in this case, 32Gb), the path for the reference genome (you may use the provided chr6.fasta).
```

Use the provided script "filter_after_VQSR.pl" and file "MHC_All.vcf" to filter out the artifacts.

We will call the new VCF file after the VQSR procedure as VCF.VQSR.vcf


Use vcfx to introduce missing alleles in unbalanced heterozygous sites and in homozygous sites in regions with very low read depth. Please check the vcfx manual to understand what is going on here (www.castelli-lab.net/apps/vcfx).

vcfx checkad input=VCF.VQSR.vcf (this will create a .ad.vcf file next to the original VCF)

Use bcftools to remove alleles that no longer exist, and vcftools to recode the file.
> bcftools view --trim-alt-alleles VCF.VQSR.ad.vcf > VCF.VQSR.ad.trim.vcf
> bcftools view --min-ac 1 VCF.VQSR.ad.trim.vcf > VCF.VQSR.ad.trim.minac.vcf
> vcftools --vcf VCF.VQSR.ad.trim.minac.vcf --recode --out VCF.VQSR.ad.trim.minac.rec.vcf


The last VCF file contains only the variants that have passed the VQSR/vcfx workflow. For now on, we will refer to this VCF file as "VCF".
+ Please note that the VCF generated up to this step is suitable for association studies and other purposes. Still, it consists of unphased genotypes with some missing alleles.


## STEP 5 - Calling phasing sets directly from the sequencing data
In this step, we will infer phase sets (the micro haplotypes) directly from the sequencing data using WhatsHap. We will use these phase sets in the upcoming haplotyping procedure with shapeit4.
We have two options here. 
The first option is to run WhatsHap as recommended, using a single core, such as this:
whatshap phase --indels -o whatshap.vcf --reference chr6.fasta --tag PS VCF.VQSR.ad.trim.minac.rec.vcf bam1 bam2 bam3 …
However, we recommend the next option because it uses better your computation resources and maximizes the WhatsHap phasing capacity by transforming some multi-allelic variants into bi-allelic ones.


The second (and recommended) option is to run WhatsHap in parallel using the provided script (parallelize_whatshap.pl). This script will split your VCF files into single-sample VCF files, call WhatsHap for each sample in parallel, and join all files in a single VCF.


The input for this step is the VCF file produced in step 4 (VCF.VQSR.ad.trim.minac.rec.vcf). The output is a VCF file with phase sets (whatshap.vcf)
An example of the script to parallelize whatshap:
perl parallelize_whatshap.pl -v VCF.VQSR.ad.trim.minac.rec.vcf -b /Users/lab/bams/ -o /Users/lab/whatshap_out/ -r chr6.fasta


You will find a whatshap.vcf file in the output folder.


- Attention: You should copy all the .adjusted.bam and .adjusted.bam.bai files from each hla-mapper output to the same location, and indicate this location using the -b option.


STEP 6 - Normalize your WhatsHap VCF to a biallelic VCF
bcftools norm -m-any whatshap.vcf > whatshap.biallelic.vcf
Use bgzip and tabix to compress and index the whatshap.biallelic.vcf file
attention: you should use bcftools 1.13


STEP 7 - Calling haplotypes
We will use shapeit4 to call haplotypes. Please check https://odelaneau.github.io/shapeit4/ for instructions on how to do it.
An example of this run is as follows:
shapeit4 --input whatshap.biallelic.vcf.gz --map chr6.b38.gmap.gz --region chr6 --output whatshap.biallelic.shapeit.vcf --thread 10 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --sequencing --use-PS 0.0001 
The final VCF file is a phased biallelic VCF.
attention: you may adjust the number of threads, the interactions scheme. Include the map file --map is optional. However, never forget to include --sequencing and --use-PS

STEP 8 - Convert the biallelic VCF to multi-allelic VCF
bcftools norm -m+any whatshap.biallelic.shapeit.vcf  > whatshap.biallelic.shapeit.multi.vcf 
attention: you should use bcftools 1.13
+ Please note that the VCF generated up to this step is suitable for association studies and other purposes. 


STEP 9 - Calling complete sequences and HLA alleles
For this step, we will generate complete sequences for each HLA gene, and compare them with known ones from the IPD-IMGT/HLA database.
You should use the provided script hla-mapper_call_alleles.pl
The first step is compressing and indexing the VCF file from step 8 using BGZIP and TABIX.
A typical run would be something like this:
perl hla-mapper_call_alleles.pl -v whatshap.biallelic.shapeit.multi.vcf.gz -p HLA-B -o /home/lab/allele_out/ -r chr6.fasta

-v [the VCF file, compressed and indexed with BGZIP and TABIX)
-p [the profile, such as HLA-B]
-o [where the outputs should be placed]
-r [the reference sequence for chr6]


You can check the available profiles in the /bed folder. You may also create a new profile by adding/replacing files following the pattern observed in the /bed folder. 
The /imgt_db/HLA folder contains the IPD-IMGT/HLA database sequences. You can update this data if you want.
The script will create the output folder and a folder for each profile. These are the files in the output folder, using HLA-B as the profile:
HLA-B.genomic.fas
A fasta file with two complete sequences for every individual, one per chromosome [h1 and h2]. You should ignore this file when dealing with exomes, or when your sequencing data does not contain intronic and UTR data.


HLA-B.genomic.counted.fas
A fasta file with one copy of each different sequence, their names according to the IPD-IMGT/HLA database, or an indication that they are new, followed by their size and how many times they were detected in your dataset. You should ignore this file when dealing with exomes or when your sequencing data does not contain intronic and UTR data.


HLA-B.vcf
This is a copy of the VCF file but containing only the HLA-B region.


HLA-B.genomic.groups.fas
If necessary, alleles with the same sequence will be indicated here because you are dealing with a smaller region.


HLA-B.exon.fas
A fasta file with two exonic sequences for every individual, one per chromosome [h1 and h2]. Exons and concatenated into a single sequence. 


HLA-B.cds.fas
A fasta file with two CDS sequences for every individual, one per chromosome [h1 and h2]. These sequences start at the first translated ATG and end at the stop codon. Exons and concatenated into a single sequence. 
 
HLA-B.cds.counted.fas
A fasta file with one copy of each different CDS sequence, their names according to the IPD-IMGT/HLA database, or an indication that they are new, followed by their size and how many times they were detected in your dataset. 


HLA-B.prot.fas
A fasta file with two protein sequences for every individual, one per chromosome [h1 and h2]. They are a translation of the ones in the HLA-B.cds.fas. 
 
HLA-B.prot.counted.fas
A fasta file with one copy of each different protein sequence, their names according to the IPD-IMGT/HLA database, or an indication that they are new, followed by their size and how many times they were detected in your dataset. 


HLA-B.db.txt
A tab separated file with the results for all samples. Please ignore the Genomic data if your data does not contain introns. There is also an indication of the size of each sequence.


STEP 10 - Checking the alleles
At this step, please check the .db.txt file and evaluate with this data is as you expected. There are too many new alleles, and you did not expect that? Maybe an artifact may be passed the VQSR/vcfx workflow and needs to be manually removed, or something correct didn't pass the filer and must be manually recovered. For instance, if there are many samples with the same new allele, you may inspect the BAM file of one of these samples and the VCF using IGV to detect if there is a problem. If this is the case, you should perform the necessary adjustments and go back to step 5.


Known issues
Sometimes there is capture bias for exomes and panels, and some of the chromosomes are not captured or are under-captured. This might lead to genotyping errors beyond this pipeline's scope. This is frequently observed for HLA-DQA1 and HLA-DQB1 and exomes. We have solved this problem with imputation using the data from step 4 (please refer to doi 10.3389/fimmu.2022.975918).

