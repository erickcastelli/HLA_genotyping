# HLA Genotypingm haplotyping, and allele calls from next-generation sequencing
Tutorial for genotyping, haplotyping, and allele calls for HLA genes

Version 1.0 (May 25th 2021)

## Packages and software needed

- hla-mapper 4 (www.castelli-lab.net/apps/hla-mapper) 
- GATK 4 (https://gatk.broadinstitute.org/hc/en-us)
- vcfx 2 (www.castelli-lab.net/apps/vcfx)
- phasex 0.8.2 (https://github.com/erickcastelli/phasex)
- samtools 1.12 (http://samtools.sourceforge.net)
- BWA 0.7.17 (https://sourceforge.net/projects/bio-bwa/files/)
- bcftools 1.12 (http://samtools.github.io/bcftools/)
- vcftools (http://vcftools.sourceforge.net)
- IGV (https://software.broadinstitute.org/software/igv/)

 
## STEP 1A (IF YOU WANT TO ALSO GENOTYPE INTERGENIC VARIANTS): 
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

## STEP 1B (IF YOU WILL GENOTYPE ONLY THE HLA GENES): 
> hla-mapper dna r1=R1.fastq.gz r2=R2.fastq.gz db=hla_mapper_database sample=Sample_Name output=output_folder
- You need to repeat this last part for each sample.
- You need to indicate a different output folder for each sample

## STEP 2 - Check some of the hla-mapper BAM files using IGV 
Using IGV, please check some of the hla-mapper outputed BAM files (Sample_Name.adjusted.bam) using IGV.

## STEP 3 - Variant call using GATK 4
We recommend GATK 4 HaplotypeCaller to call variants. The hla-mapper BAM file is already prepared for GATK.

For each sample, run GATK HaplotypeCaller in the GVCF mode, such as this example:

> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -R reference_genome -I Sample_name.adjusted.bam -O output_folder/Sample_name.MHC.g.vcf -L chr6:29700000-33150000 -ERC GVCF --max-num-haplotypes-in-population 256 --native-pair-hmm-threads thread_number --allow-non-unique-kmers-in-ref TRUE

You need to ajust the amont of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the hla-mapper output BAM (Sample_name.adjusted.bam), the output folder (output_folder), the sample name, and the number of threads (thread_number).

After processing all your samples, you need to concatenate all G.VCF files in a single one. There are two ways to do that, depending on the number of samples. The most commom one is using GATK CombineGVCFs (https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs). The other is GenomicsDBImport (https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport). This tutorial does not cover this issue. Please follow the GATK instructions and combine all GVCFS in a single file.

Now, you can genotype your GVCF using GATK GenotypeGVCFs, as follows:

> java -Xmx32g -jar gatk-package-4.2.0.0-local.jar GenotypeGVCFs -R reference_genome -O output_folder/MHC.vcf -L chr6:29700000-33150000 --variant output_folder/All_samples.MHC.g.vcf --dbsnp path_to_dbsnp_vcf

You need to ajust the amont of memory (in this case, 32Gb), the path for the reference genome (reference_genome), the path for the single GVCF file (output_folder/All_samples.MHC.g.vcf), the output folder (output_folder), and the path to dbsnp (path_to_dbsnp_vcf). dbSNP is optional.

## STEP 4 - Variant refinement
There are many ways to proceed with variant refinement, i.e., the removal of artifacts, including GATK VQRS and vcfx. For HLA genes, we recommend vcfx.

Recode the VCF file using vcftools. This is important for the next steps.

> vcftools --vcf VCF_FILE --recode --out OUTPUT_FOLDER/VCF_FILE

Using any application or script you want, please change any "|" allele separator for "/" in the recoded VCF file.

Use vcfx to filter out artifacts and variants with too many missing alleles, as follows. Please check the vcfx manual to understand out is going on here (www.castelli-lab.net/apps/vcfx)

> vcfx checkpl input=VCF_FILE (this will create a .pl.vcf file next to the original VCF)
> 
> bcftools view --trim-alt-alleles --min-ac 1 VCF_FILE.pl.vcf > VCF_FILE.pl.trim.vcf
> 
> vcfx evidence input=VCF_FILE.pl.trim.vcf (this will create a .pl.trim.evid.vcf)
> 
> vcfx filter input=VCF_FILE.pl.trim.evid.vcf tag=PASS,WARN (this will create a .pl.trim.evid.filter.vcf)

The last VCF file contains only the variants that have passed the vcfx checkpl/evidence algorithm.


## STEP 5 - Haplotype calls

