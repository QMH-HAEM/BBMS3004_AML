#!/bin/sh
###################################################################################
# This rudimentary pipeline is for the teaching purpose in BBMS3004 Targeted
# Sequencing in Acute Myeloid Leukaemia. It follows the GATK Best
# Practice for somatic SNV and indel to call variants in DNA samples. As some
# elements of this pipeline has been abridged to suit the computer resources and
# running time in the sessions, this pipeline and its associated reference files
# should NOT be used for any other purposes, including clinical and research
# purposes, it is likely to yield incorrect results if used otherwise.
#
# This script is designed to run in the directory:
# /home/ubuntu/Documents/practical3
#
# Division of Haematology, Department of Pathology, Queen Mary Hospital
# Authors: Wing-Fai Tang, Ho-Wan Ip
#
# Created on 02/08/2019
# Modified on 21/01/2021
# Modified on 11/02/2025
###################################################################################

set -e
set -u
set -o pipefail

if [ "$#" -lt 1 ]
then
   echo "You have input too few arguments. Correct command usage: bash toyPipeline.sh sample_number"
   exit 1
fi

# Define sample name
SN=$1

# Perform Trimmomatic
cp -n ~/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa ./
java -jar ~/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${SN}_1.fastq.gz ${SN}_2.fastq.gz ${SN}_1_trimmed_paired.fq.gz ${SN}_1_trimmed_unpaired.fq.gz ${SN}_2_trimmed_paired.fq.gz ${SN}_2_trimmed_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40

# Perform BWA
bwa mem -t 4 -M ~/Ref/ucsc.hg19.fasta ${SN}_1_trimmed_paired.fq.gz ${SN}_2_trimmed_paired.fq.gz > ${SN}.sam

# Perform sorting on SAM file and convert to BAM
gatk SortSam -I ${SN}.sam -O ${SN}_sorted.bam --SORT_ORDER coordinate

# Collect statistics from BAM file
bamtools stats -insert -in ${SN}_sorted.bam > ${SN}_sorted.bamtools.stats

# Perform adjustment procedures
gatk AddOrReplaceReadGroups -I ${SN}_sorted.bam -O ${SN}_RG.bam --RGID SPACE --RGLB panel --RGPL ILLUMINA --RGPU unit1 --RGSM ${SN}
gatk MarkDuplicates -I ${SN}_RG.bam -O ${SN}_MD.bam  -M ./${SN}_MD.stats --CREATE_INDEX true
gatk BaseRecalibrator -R ~/Ref/ucsc.hg19.fasta -I ${SN}_MD.bam -L ~/Ref/myeloid-targets.interval_list -ip 50 --known-sites ~/Ref/dbsnp_138.hg19.vcf --known-sites ~/Ref/Mills_and_1000G_gold_standard.indels.hg19.vcf -O ${SN}_recal_data.table
gatk ApplyBQSR -R ~/Ref/ucsc.hg19.fasta -I ${SN}_MD.bam --bqsr-recal-file ${SN}_recal_data.table -O ${SN}_BR.bam

# Collect QC metrics
gatk CollectMultipleMetrics -I ${SN}_BR.bam -O ./${SN}_GATK
gatk CollectReadCounts -I ${SN}_BR.bam -L ~/Ref/myeloid-targets.interval_list --interval-merging-rule OVERLAPPING_ONLY --format TSV -O ./${SN}.counts.tsv
gatk CollectHsMetrics -I ${SN}_BR.bam -O ./${SN}_hs_metrics.txt -R ~/Ref/ucsc.hg19.fasta -BI ~/Ref/myeloid-probe-coords.interval_list -TI ~/Ref/myeloid-targets.interval_list

# Perform variant calling by Mutect2
gatk Mutect2 -R ~/Ref/ucsc.hg19.fasta -I ${SN}_BR.bam -tumor ${SN} -L ~/Ref/myeloid-targets.interval_list  -germline-resource ~/Ref/af-only-gnomad.myeloid.bedtools.vcf.gz --f1r2-tar-gz f1r2.tar.gz -O ${SN}_unfiltered.vcf
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
gatk GetPileupSummaries -I ${SN}_BR.bam -V ~/Ref/small_exac_common_myeloid.vcf.gz -L ~/Ref/small_exac_common_myeloid.vcf.gz  -O getpileupsummaries.table
gatk CalculateContamination -I getpileupsummaries.table  -O calculatecontamination.table
gatk FilterMutectCalls -R ~/Ref/ucsc.hg19.fasta -V ${SN}_unfiltered.vcf --contamination-table calculatecontamination.table --ob-priors read-orientation-model.tar.gz -O ${SN}_filtered.vcf

# Filter the variants with VAF less than 0.05
bcftools filter -i " FORMAT/AF > 0.05 " ${SN}_filtered.vcf -o ${SN}_filtered_0.05.vcf

# Perform variant annotation by ANNOVAR
perl ~/Programs/annovar/table_annovar.pl ${SN}_filtered_0.05.vcf ~/Programs/annovar/humandb/ -buildver hg19 -out ${SN}_filtered_annotate -remove -protocol refGene,cosmic86,clinvar_20170905,exac03nontcga,gnomad_exome -operation g,f,f,f,f -nastring . -vcfinput
