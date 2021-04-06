# BBMS3004 Molecular Diagnostics Laboratory (2021)
## Practical: Targeted Sequencing and Its Clinical Application in Acute Myeloid Leukaemia

In this session. we will build a toy bioinformatic pipeline to get a feel of what it is like to practice bioinformatics using Linux. Our platform will be [Azure Lab Services](https://labs.azure.com/), which is a Microsoft cloud computing platform.

This session does not assume any Linux command-line skills from you. We will try to learn a few important Linux commands on-the-fly while we learn how to build our toy bioinformatic pipeline.

### Before the Practical Session...

Everyone of you should have received an invitation email from Microsoft Azure to register for this lab. In the email, click the tab "Register for the lab>" and follow the instructions to complete the registration process. If you are successful, you will arrive at a Azure Lab Services webpage with the header saying "My virtual machines". On the webpage, there is a virtual machine "DNAseqAML" waiting for you. If you reach this, you are prepared for the practical session.

You will notice the virtual machine is "Stopped". At this stage, please do not start the machine yet, as you will consume the assigned time for the session.

### Login to Azure Lab Services and Get Things Started

Before we begin with our discussion, let's take the time to initialise our Azure Lab Services and set our passwords, as these will take some time before you can login to the system.

1. Login to [My Virtual Machines](https://labs.azure.com/virtualmachines) in Azure Lab Services.
2. Click the switch of the DNAseqAML virtual machine to start it. This will take a few minutes. (If you are logging in for the first time, the system will ask you to set a new password. Proceed accordingly.)
3. When the virtual machine shows "Running", click on the little monitor at the lower right hand corner of the virtual machine tab. Click "Connect via SSH". A link will be shown.
4. Login to the virtual machine by PuTTY (for Windows) or in the terminal (for Mac OS) according to the instruction in class.
5. When you finish using the virtual machine. Type "exit" into the Linux shell to quit the session AND click the switch of the DNAseqAML in the webpage to stop it.

After you have login to our virtual machine, please follow the instructions in class to explore the environment.

### Basic Linux Commands
This section summarises the important Linux commands that will be useful during the practical session.

**Command: Show present working directory**
```bash
pwd
```

**Command: List files and sub-directories in the current directory**
```bash
ls
```

**Command: Change directory**
```bash
cd subdirectory               # Change to sub-directory
cd ..                         # Change to parent directory
```

**Command: Make a new directory**
```bash
mkdir newdirectory
```

### GATK Best Practice for Somatic SNVs and Indels

We will try to build a pipeline according [GATK Best Practices for somatic SNVs and indels](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146). Other GATK Best Practices pipelines can be found in the menu of the website as well.

### Building the Toy Pipeline

While we are trying to build a useable bioinformatic pipeline, but as some elements of this pipeline has been abridged to suit the computer resources and running time in the teaching session, this pipeline should NOT be used for any other purposes, including clinical practice and research.

**Change directory to the correct folder and create directory for QC statistics output**
```bash
cd Documents
mkdir stats
```

**Perform adaptor and quality trimming using Trimmomatic**

Run [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```bash
cp ~/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa ./

java -jar ~/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 17H0510082_1.fastq.gz 17H0510082_2.fastq.gz 17H0510082_1_trimmed_paired.fq.gz 17H0510082_1_trimmed_unpaired.fq.gz 17H0510082_2_trimmed_paired.fq.gz 17H0510082_2_trimmed_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40
```
*Please try to run fastqc again on the 2 adaptor and quality trimmed paired fastq files.

**Perform sequence alignment by BWA**

Sequencing mapping by [BWA](http://bio-bwa.sourceforge.net/)

```bash
bwa mem -t 4 -M ~/Ref/ucsc.hg19.fasta 17H0510082_1_trimmed_paired.fq.gz 17H0510082_2_trimmed_paired.fq.gz > 17H0510082.sam
```

**Perform sorting on SAM file and convert to BAM**
```bash
gatk SortSam -I 17H0510082.sam -O 17H0510082_sorted.bam --SORT_ORDER coordinate
```
*Inspect content of SAM and BAM files

**Inspecting a BAM file**
```bash
samtools view 17H0510082_BR.bam | less -S
```

**Collect statistics from BAM file**
Determine mappable read number:
```bash
bamtools stats -insert -in 17H0510082_sorted.bam > 17H0510082_sorted.bamtools.stats
```

**Perform adjustment procedures**
```bash
gatk AddOrReplaceReadGroups -I 17H0510082_sorted.bam -O 17H0510082_RG.bam --RGID SPACE --RGLB panel --RGPL ILLUMINA --RGPU unit1 --RGSM 17H0510082
gatk MarkDuplicates -I 17H0510082_RG.bam -O 17H0510082_MD.bam  -M ./stats/17H0510082_MD.stats --CREATE_INDEX true
gatk BaseRecalibrator -R ~/Ref/ucsc.hg19.fasta -I 17H0510082_MD.bam -L ~/Ref/myeloid-targets.interval_list -ip 50 --known-sites ~/Ref/dbsnp_138.hg19.vcf --known-sites ~/Ref/Mills_and_1000G_gold_standard.indels.hg19.vcf -O 17H0510082_recal_data.table
gatk ApplyBQSR -R ~/Ref/ucsc.hg19.fasta -I 17H0510082_MD.bam --bqsr-recal-file 17H0510082_recal_data.table -O 17H0510082_BR.bam
```

**Collect QC metrics**
```bash
gatk CollectMultipleMetrics -I 17H0510082_BR.bam -O ./stats/17H0510082_GATK
gatk CollectReadCounts -I 17H0510082_BR.bam -L ~/Ref/myeloid-targets.interval_list --interval-merging-rule OVERLAPPING_ONLY --format TSV -O ./stats/17H0510082.counts.tsv
gatk CollectHsMetrics -I 17H0510082_BR.bam -O ./stats/17H0510082_hs_metrics.txt -R ~/Ref/ucsc.hg19.fasta -BI ~/Ref/myeloid-probe-coords.interval_list -TI ~/Ref/myeloid-targets.interval_list
```

**Perform variant calling by Mutect2**
```bash
gatk Mutect2 -R ~/Ref/ucsc.hg19.fasta -I 17H0510082_BR.bam -tumor 17H0510082 -L ~/Ref/myeloid-targets.interval_list  -germline-resource ~/Ref/af-only-gnomad.myeloid.bedtools.vcf.gz --f1r2-tar-gz f1r2.tar.gz -O 17H0510082_unfiltered.vcf
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
gatk GetPileupSummaries -I 17H0510082_BR.bam -V ~/Ref/small_exac_common_myeloid.vcf.gz -L ~/Ref/small_exac_common_myeloid.vcf.gz  -O getpileupsummaries.table
gatk CalculateContamination -I getpileupsummaries.table  -O calculatecontamination.table
gatk FilterMutectCalls -R ~/Ref/ucsc.hg19.fasta -V 17H0510082_unfiltered.vcf --contamination-table calculatecontamination.table --ob-priors read-orientation-model.tar.gz -O $17H0510082_filtered.vcf
```
*Inspect content of VCF file.

**Filter variants**
```bash 
bcftools filter -i " FORMAT/AF > 0.05 " 17H0510082_filtered.vcf -o 17H0510082_filtered_0.05.vcf
```

**Perform variant annotation by ANNOVAR**
```bash
perl ~/Programs/annovar/table_annovar.pl 17H0510082_filtered_0.05.vcf ~/Programs/annovar/humandb/ -buildver hg19 -out 17H0510082_filtered_annotate -remove -protocol refGene,cosmic86,clinvar_20170905,exac03nontcga,gnomad_exome -operation g,f,f,f,f -nastring . -vcfinput
```
*Inspect variants in IGV.

### A Wrapped Bioinformatic Pipeline

To run the above steps in a single command, you can run this [shell script](https://github.com/QMH-HAEM/clinical-bioinformatics-3/raw/master/toyPipeline.sh) (please right click and save as in the directory of "practical3"), which contains all the commands above.

**Running a shell script**
```bash
bash toyPipeline.sh 17H0510082
```
