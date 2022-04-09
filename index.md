# BBMS3004 Molecular Diagnostics Laboratory (2022)
## Practical: Targeted Sequencing and Its Clinical Application in Acute Myeloid Leukaemia

In this session. we will build a toy bioinformatic pipeline to get a feel of what it is like to practice bioinformatics using Linux. Our platform will be [Azure Lab Services](https://labs.azure.com/), which is a Microsoft cloud computing platform.

This session does not assume any Linux command-line skills from you. We will try to learn a few important Linux commands on-the-fly while we learn how to build our toy bioinformatic pipeline. But due to COVID-19 wave, this session will not be held in a face-to-face manner. This may increase the difficulties a little bit when you start using the platform, but we will organise a face-to-face Zoom question-and-answer session to address any problem that you may have.

Please be assured that this practical serves to give you some hands-on experience in dealing with NGS sequencing data. You will not be asked about command-line operation in your examination.

### Before the Practical Session...

Everyone of you should have received an invitation email from Microsoft Azure to register for this lab. In the email, click the tab "Register for the lab>" and follow the instructions to complete the registration process. You should use your email with the address ending in "connect.hku.hk" to register. If you are successful, you will arrive at a Azure Lab Services webpage with the header saying "My virtual machines". On the webpage, there is a virtual machine "NGS in AML" waiting for you. If you reach this, you are prepared for the practical session. In case you come across an error message that reads "You do not have permission to access this lab...", please email me at iphowan@hku.hk and let me know your email address that ends with "connect.hku.hk".

You will notice the virtual machine is "Stopped". Please start the machine only when you are ready to go through the session, as a running machine will consume your assigned quota of 4 hours (per student).

### Login to Azure Lab Services and Get Things Started

Before we begin with our discussion, let's take the time to initialise our Azure Lab Services and set our passwords, as these will take some time before you can login to the system.

1. Login to [My Virtual Machines](https://labs.azure.com/virtualmachines) in Azure Lab Services.
2. Click the switch of the "NGS in AML" virtual machine to start it. This will take a few minutes. (If you are logging in for the first time, the system will ask you to set a new password. Proceed accordingly.)
3. When the virtual machine shows "Running", click on the little monitor at the lower right hand corner of the virtual machine tab. Click "Connect via SSH". A link will be shown.
4. Login to the virtual machine by PuTTY (for Windows) or in the terminal (for Mac OS, search for "Terminal" in Spotlight). For PuTTY, follow the instructions in the PowerPoint slides. For Mac terminal, simply paste the command (starting with "ssh ...") generated in Azure.
5. When you finish using the virtual machine. Type "exit" into the Linux shell to quit the session AND click the switch of the "NGS in AML" in the webpage to stop it.

After you have login to our virtual machine, please follow the instructions in the video to explore the Linux environment.

### Basic Linux Commands
This section summarises the important Linux commands that will be useful during the practical session. In case you are very interested in learning about Linux command-line operation, you may want to go through an [online tutorial course](https://rnabio.org/module-00-setup/0000/08/01/Unix/).

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

**Command: Copy a file**
```bash
cp ../abc.txt ./              # Copy abc.txt from the parent directory to the current directory
```

**Command: Read the contents of a file**
```bash
less -S ucsc.hg19.fasta              # View the human reference genome hg19 (in fasta format)
zless -S 17H1220080_1.fastq.gz       # View the zipped fastq file of the AML patient (in gzipped fastq format)
```

### GATK Best Practice for Somatic SNVs and Indels

We will try to build a pipeline according [GATK Best Practices for Data Preprocessing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) and [GATK Best Practices for somatic SNVs and indels](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-). Other GATK Best Practices pipelines can be found in the menu of the website as well.

### Building the Toy Pipeline

While we are trying to build a useable bioinformatic pipeline, but as some elements of this pipeline has been abridged to suit the computer resources and running time in the teaching session, this pipeline should NOT be used for any other purposes, including clinical practice and research.

**Change directory to the correct folder and create directory for QC statistics output**
```bash
cd Documents
mkdir stats
```

Before the start of analysis, the usual step is to inspect the basic quality of the sequencing data. A common software to use is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Since the output of FastQC is in html format and cannot be visualised in our current server setup, this step will be demonstrated instead. To inspect the output of the FastQC step for our data, please refer to [FastQC output for 17H1220080_1.fastq.gz](https://htmlpreview.github.io/?https://github.com/QMH-HAEM/BBMS3004_AML/blob/main/17H1220080_1_fastqc.html) and [FastQC output for 17H1220080_2.fastq.gz](https://htmlpreview.github.io/?https://github.com/QMH-HAEM/BBMS3004_AML/blob/main/17H1220080_2_fastqc.html).

**Perform adaptor and quality trimming using Trimmomatic**

Run [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```bash
cp ~/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa ./

java -jar ~/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 17H1220080_1.fastq.gz 17H1220080_2.fastq.gz 17H1220080_1_trimmed_paired.fq.gz 17H1220080_1_trimmed_unpaired.fq.gz 17H1220080_2_trimmed_paired.fq.gz 17H1220080_2_trimmed_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40
```

**Perform sequence alignment by BWA**

Sequencing mapping by [BWA](http://bio-bwa.sourceforge.net/)

```bash
bwa mem -t 4 -M ~/Ref/ucsc.hg19.fasta 17H1220080_1_trimmed_paired.fq.gz 17H1220080_2_trimmed_paired.fq.gz > 17H1220080.sam
```

**Perform sorting on SAM file and convert to BAM**
```bash
gatk SortSam -I 17H1220080.sam -O 17H1220080_sorted.bam --SORT_ORDER coordinate
```
*Inspect content of SAM and BAM files

**Inspecting a BAM file**
```bash
samtools view 17H1220080_BR.bam | less -S
```

**Collect statistics from BAM file**
Determine mappable read number:
```bash
bamtools stats -insert -in 17H1220080_sorted.bam > 17H1220080_sorted.bamtools.stats
```

**Perform adjustment procedures**
```bash
gatk AddOrReplaceReadGroups -I 17H1220080_sorted.bam -O 17H1220080_RG.bam --RGID SPACE --RGLB panel --RGPL ILLUMINA --RGPU unit1 --RGSM 17H1220080
gatk MarkDuplicates -I 17H1220080_RG.bam -O 17H1220080_MD.bam  -M ./stats/17H1220080_MD.stats --CREATE_INDEX true
gatk BaseRecalibrator -R ~/Ref/ucsc.hg19.fasta -I 17H1220080_MD.bam -L ~/Ref/myeloid-targets.interval_list -ip 50 --known-sites ~/Ref/dbsnp_138.hg19.vcf --known-sites ~/Ref/Mills_and_1000G_gold_standard.indels.hg19.vcf -O 17H1220080_recal_data.table
gatk ApplyBQSR -R ~/Ref/ucsc.hg19.fasta -I 17H1220080_MD.bam --bqsr-recal-file 17H1220080_recal_data.table -O 17H1220080_BR.bam
```

**Collect QC metrics**
```bash
gatk CollectMultipleMetrics -I 17H1220080_BR.bam -O ./stats/17H1220080_GATK
gatk CollectReadCounts -I 17H1220080_BR.bam -L ~/Ref/myeloid-targets.interval_list --interval-merging-rule OVERLAPPING_ONLY --format TSV -O ./stats/17H1220080.counts.tsv
gatk CollectHsMetrics -I 17H1220080_BR.bam -O ./stats/17H1220080_hs_metrics.txt -R ~/Ref/ucsc.hg19.fasta -BI ~/Ref/myeloid-probe-coords.interval_list -TI ~/Ref/myeloid-targets.interval_list
```

**Perform variant calling by Mutect2**
```bash
gatk Mutect2 -R ~/Ref/ucsc.hg19.fasta -I 17H1220080_BR.bam -tumor 17H1220080 -L ~/Ref/myeloid-targets.interval_list  -germline-resource ~/Ref/af-only-gnomad.myeloid.bedtools.vcf.gz --f1r2-tar-gz f1r2.tar.gz -O 17H1220080_unfiltered.vcf
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
gatk GetPileupSummaries -I 17H1220080_BR.bam -V ~/Ref/small_exac_common_myeloid.vcf.gz -L ~/Ref/small_exac_common_myeloid.vcf.gz  -O getpileupsummaries.table
gatk CalculateContamination -I getpileupsummaries.table  -O calculatecontamination.table
gatk FilterMutectCalls -R ~/Ref/ucsc.hg19.fasta -V 17H1220080_unfiltered.vcf --contamination-table calculatecontamination.table --ob-priors read-orientation-model.tar.gz -O 17H1220080_filtered.vcf
```
*Inspect content of VCF file.

**Filter variants**
```bash 
bcftools filter -i " FORMAT/AF > 0.05 " 17H1220080_filtered.vcf -o 17H1220080_filtered_0.05.vcf
```

**Perform variant annotation by ANNOVAR**
```bash
perl ~/Programs/annovar/table_annovar.pl 17H1220080_filtered_0.05.vcf ~/Programs/annovar/humandb/ -buildver hg19 -out 17H1220080_filtered_annotate -remove -protocol refGene,cosmic86,clinvar_20170905,exac03nontcga,gnomad_exome -operation g,f,f,f,f -nastring . -vcfinput
```
*Inspect variants in IGV.

You may find it a bit difficult to view the VCF file in the command-line interface. Another option is to download the [annotated file](https://github.com/QMH-HAEM/BBMS3004_AML/raw/main/17H1220080_filtered_annotate.hg19_multianno.txt) (right click and select "Save link as...") to your computer and view it using your spreadsheet viewer.

### A Wrapped Bioinformatic Pipeline

To run the above steps in a single command, you can run the toyPipeline.sh shell script, which contains all the commands you have entered above, in the following manner.

**Running a shell script**
```bash
bash toyPipeline.sh 17H1220080
```

### Visualisation of Sequencing Data in IGV

To visualise the sequencing data, you can download the [BAM file](https://github.com/QMH-HAEM/BBMS3004_AML/raw/main/17H1220080_BR.bam) produced after base recalibration, along with its [index file](https://github.com/QMH-HAEM/BBMS3004_AML/raw/main/17H1220080_BR.bai), and load the BAM file into [IGV](https://software.broadinstitute.org/software/igv/download) for visualisation. Please make sure you have selected hg19 human reference genome for visualisation. You may input the chromosome number and genomic coordinate of a variant in IGV to visualise it, e.g. for the IDH1 R132 variant, you can enter chr2:209113112.

### After the session

Thank you for going through this session. If you have any inquiry regarding the course content or issues with running the server, please do not hesitate to email me. After the date of our session, the virtual machines will be kept open for your access until after the examination. Do feel free to explore the Linux platform until your assigned time has been used up.
