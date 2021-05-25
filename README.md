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
> hla-mapper dna bam=sample.bam db=hla_mapper_database sample=Sample_Name


## STEP 1B (IF YOU WILL GENOTYPE ONLY THE HLA GENES): 
> hla-mapper dna r1=R1.fastq.gz r2=R2.fastq.gz db=hla_mapper_database sample=Sample_Name


