---
title: "Variant calling workflow"
author: "Mary Piper, Meeta Mistry"
date: "Thursday, March 31, 2016"
---

Approximate time: 120 minutes

## Learning Objectives:

* Understand the different steps involved in inferring expression
* Use a series of command line tools to execute an expression workflow
* Becoming familiar with data formats encountered during NGS analysis
* Automate a workflow by grouping a series of sequential commands into a script


## Setting up

To get started with this lesson, make sure you are in `dc_workshop/data`. We've ignored it so far, but there is also a reference transcriptome we will be using in the `reference` directory.

Your directory structure should now look like this:

<pre>
dc_workshop
├── data
    ├── ref_genome
        └── GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz
    ├── untrimmed_fastq
    	├── SRR2040575_brain_1.fastq.gz
	├── SRR2040575_brain_2.fastq.gz
	├── SRR2040576_brain_1.fastq.gz
	├── SRR2040576_brain_2.fastq.gz
	├── SRR2040577_heart_1.fastq.gz
	├── SRR2040577_heart_2.fastq.gz
	├── SRR2040578_heart_1.fastq.gz
	├── SRR2040578_heart_2.fastq.gz
        ├── ...
    └── trimmed_fastq
	├── SRR2040575_brain_1_1P_fastqc.html
	├── SRR2040575_brain_1_1P_fastqc.zip
	├── SRR2040575_brain_1_1P.fastq.gz
	├── SRR2040575_brain_1_1U_fastqc.html
	├── SRR2040575_brain_1_1U_fastqc.zip
	├── SRR2040575_brain_1_1U.fastq.gz
    	├── ...
 ├── results
    ├── ...
 └── docs
    ├── ...

</pre>

You will also need to create directories for the results that will be generated as part of the workflow: 
```bash
$ mkdir results/sam results/bam results/counts
```

> *NOTE: All of the tools that we will be using in this workflow have been pre-installed on our remote computer*

## Alignment to a reference genome

We have already trimmed our reads so now the next step is alignment of our quality reads to the reference genome.

![workflow_align](../img/expression_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to choose from and while there is no gold standard there are some tools that are better suited for particular NGS analyses. We will be using the [Burrows Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent sequences against a large reference genome. The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

### Index the reference genome
Our first step is to index the reference genome for use by BWA. *NOTE: This only has to be run once*. The only reason you would want to create a new index is if you are working with a different reference  genome or you are using a different tool for alignment.

```bash   
$ # This step helps with the speed of alignment
$ bwa index data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz     
```

Eventually we will loop over all of our files to run this workflow on all of our samples, but for now we're going to work with the quality filtered paired sequences only for just one sample in our dataset `SRR2040575_brain_1_1P.fastq.gz` and `SRR2040575_brain_1_2P.fastq.gz`:

```bash
$ ls -alh ~/dc_workshop/data/trimmed_fastq/SRR2040575_brain_1_[12]P.fastq.gz 
```

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an aligner. BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate, so we will be using that. 
    
Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using and try to understand the options.*

```bash
$ bwa mem data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz \
    data/trimmed_fastq/SRR2040575_brain_1_1P.fastq.gz data/trimmed_fastq/SRR2040575_brain_1_2P.fastq.gz \
    > results/sam/SRR2040575_brain.aligned.sam
```
### Convert the format of the alignment to SAM/BAM

Explore the information within your SAM file:

```bash
$ head results/sam/SRR2040575_brain.aligned.sam
```	
Now convert the SAM file to BAM format for use by downstream tools: 

```bash
$ samtools view -S -b results/sam/SRR2040575_brain.aligned.sam > results/bam/SRR2040575_brain.aligned.bam
```

#### SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam), is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not have time to go in detail of the features of the SAM format, the paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification. **The binary version of SAM is called a BAM file.**

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.

![sam_bam1](../img/sam_bam.png)

![sam_bam2](../img/sam_bam3.png)

## Alignment cleanup

![workflow_clean](../img/expression_workflow_cleanup.png)

Post-alignment processing of the alignment file includes:

1. Sorting the BAM file by coordinate
2. Indexing the BAM file
3. Calculating the mapping statistics
4. Calculating the mapping depth (i.e., expression) per gene

### Sort BAM file by coordinates

Sort the BAM file:

```bash
$ samtools sort results/bam/SRR2040575_brain.aligned.bam results/bam/SRR2040575_brain.aligned.sorted
```

*SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.*

### Indexing the BAM file

Index the BAM file so it is easier to work with:

```bash
$ samtools index results/bam/SRR2040575_brain.aligned.sorted.bam
```

### Gather BAM mapping statistics

```bash
$ samtools flagstat results/bam/SRR2040575_brain.aligned.sorted.bam > results/bam/SRR2040575_brain.aligned.sorted.stats.txt
$ cat results/bam/SRR2040575_brain.aligned.sorted.stats.txt
1546135 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
1192782 + 0 mapped (77.15%:-nan%)
1546135 + 0 paired in sequencing
773100 + 0 read1
773035 + 0 read2
1179079 + 0 properly paired (76.26%:-nan%)
1187365 + 0 with itself and mate mapped
5417 + 0 singletons (0.35%:-nan%)
8060 + 0 with mate mapped to a different chr
2261 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Calculate the mapping depth for each reference sequence

RNAseq essentially randomly sequences mRNAseq data within a tissue. Therefore, gene expression is essentially a count of the number of sequencing reads that map to a given reference sequence. We can easily calculate this, and let's put it in our `results` directory.

```bash
$ samtools idxstats results/bam/SRR2040575_brain.aligned.sorted.bam > results/counts/SRR2040575_brain.aligned.sorted.counts.txt
$ head results/counts/SRR2040575_brain.aligned.sorted.counts.txt
NM_001005484.1_OR4F5	918	0	0
NM_001005221.2_OR4F29	939	0	0
NM_001005277.1_OR4F16	939	0	0
NM_198317.2_KLHL17	2581	8	0
NM_001160184.1_PLEKHN1	2299	0	0
NM_032129.2_PLEKHN1	2404	0	0
NM_001291366.1_PERM1	3417	0	0
NM_001291367.1_PERM1	3064	0	0
NM_021170.3_HES4	962	3	1
NM_001142467.1_HES4	1053	2	0
```
The resulting table contains the transcript ID in the first column, the length of that sequence, the number of mapped reads, and the number of unmapped reads.

## Using pipling to reduce input/output

One awesome feature of *BWA* and *samtools* is that they allow piping, so we can easily string commands together to more efficiently perform our mapping. This reduces the amount of input/output of data to our hard drive, making the process significantly faster. Moreover, it helps save disk space, which can rapidly disappear with the large file sizes produced in genomic data. This is just like the piping commands that you learned in the introduction to Bash, making these software very modular and powerful.

Therefore, instead of mapping the `*.fastq` files to produce `*.sam` files, converting them to binary `*.bam` files, and then sorting to produce a final `*.bam` file, we will be able to just go from the `*.fastq` files to the sorted `*.bam` files using one command! With this in mind, let's just get rid of our directory for `*.sam` files.

```bash
$ # removed recursively
$ rm -r results/sam
```

Now we can run out single mapping command by stringing the commands from above together with some special syntax.

```bash
$ bwa mem data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz \
    data/trimmed_fastq/SRR2040575_brain_1_1P.fastq.gz data/trimmed_fastq/SRR2040575_brain_1_2P.fastq.gz | \
    samtools view -S -b - | \
    samtools sort - results/bam/SRR2040575_brain.aligned.sorted
```

Notice a couple things about this command:
1. The pipe (`|`) symbols between the commands, which basically tell the computer to use the output from the first command (called STDOUT) as the input for the second command (called STDIN), without writing it to the disk.
2. The dashes (`-`) where input files names would normally be, which tells the respective command to look for input data from STDIN. And no output file names are specified except for the last command, because they are not needed until that point.
    
## Putting everything together into a mapping script

That was a lot of work, yes? But you have more FASTQ files to go...

- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?

### Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...

#### Collect all required commands
Using your command `history`, create a script file that will repeat these commands
for you. Name your script `expression_pipeline.sh`.

```bash
# aligning reads to produce sorted BAM file (1 piped command)
bwa mem data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz \
    data/trimmed_fastq/SRR2040575_brain_1_1P.fastq.gz data/trimmed_fastq/SRR2040575_brain_1_2P.fastq.gz | \
    samtools view -S -b - | \
    samtools sort - results/bam/SRR2040575_brain.aligned.sorted
    
# indexing the BAM file
samtools index results/bam/SRR2040575_brain.aligned.sorted.bam

# calculating mapping stats
samtools flagstat results/bam/SRR2040575_brain.aligned.sorted.bam > results/bam/SRR2040575_brain.aligned.sorted.stats.txt

# generating count data
samtools idxstats results/bam/SRR2040575_brain.aligned.sorted.bam > results/counts/SRR2040575_brain.aligned.sorted.counts.txt
```

### Positional parameters (input files)
At the moment this script will only call variants for `SRR2040575_brain`, but we want to grant it more flexibility so that it can be run on any sample that we give it as input.

We can specify a filename as input using **positional parameters**. Positional parameters allow flexibility within a script.

The command-line arguments $1, $2, $3,...$9 are positional parameters, with $0 pointing to the actual command, program or shell script, and $1, $2, $3, ...$9 as the arguments to the command." This basically means that "Script Name" == $0, "First Parameter" == $1, "Second Parameter" == $2 and so on...

At the start of the script let's create a variable that captures an input parameter that must be supplied with the script name, when the user is running it. This input parameter will be the name of the forward (read 1) file we want to work on.

```bash
# Get input read 1 file 
fq1=$1
```

But, wait! What about the ouput files? Or the second read file? They are all named using the SRA ID from the original FASTQ file. We can use a unix command to extract the base name of the file (without the path and .fastq.gz extension) and ue this for naming all of the reverse read and output files. We'll assign it to the `$base` variable:

```bash
# Grab base of filename for future naming
# lets get rid of _1_1P.fastq.gz
base=`basename $fq1 _1_1P.fastq.gz`
```

We can use an analogous command to get the directory name containing our target file.

```
dir=`dirname $base`
```

Now we can put together a variable for the second read.

```bash
fq2=${base}_1_2P.fastq.gz
```

Now that we have the input read names and the directory they are found in, we are ready to go. Let's initialize a script and define these variables. Let's also define a 2nd variable for the reference transcriptome.

```bash
fq1=$1
base=`basename $fq1 _1_1P.fastq.gz`
dir=`dirname $base`
fq2=${base}_1_2P.fastq.gz
ref=$2
```

Now we can create the four commands we want to run and save them to variables as well.

```bash
map="bwa mem $ref \
    ${dir}/${base}_1_1P.fastq.gz ${dir}/${base}_1_2P.fastq.gz | \
    samtools view -S -b - | \
    samtools sort - results/bam/${base}.aligned.sorted"
index="samtools index results/bam/${base}.aligned.sorted.bam"
stats="samtools flagstat results/bam/${base}.aligned.sorted.bam > results/bam/${base}.aligned.sorted.stats.txt"
counts="samtools idxstats results/bam/${base}.aligned.sorted.bam > results/counts/${base}.aligned.sorted.counts.txt"
```

We can finish this script using the material above and calling the commands we've created. It is always good to also put documentation on our script at the top. Save this script in `docs` as `expression_pipeline.sh`.

```bash
cat docs/expression_pipeline.sh
```

```bash
#!/bin/bash

# USAGE: sh run_variant_call_on_file.sh <read 1 fastq file> <reference file>
# This script will take the location and name of a fastq file and perform the following steps on it in a new directory. 
    ## map the reads to the designated reference using bwa mem, creating a sorted BAM file 
    ## index the BAM file and calcualting mapping stats 
    ## calculate expression counts for each file

# make sure we are in the directory we want
cd ~/dc_workshop

# create input/output variables
fq1=$1
base=`basename $fq1 _1_1P.fastq.gz`
dir=`dirname $fq1`
fq2=${base}_1_2P.fastq.gz
ref=$2

# create command variables
map="bwa mem $ref \
    ${dir}/${base}_1_1P.fastq.gz ${dir}/${base}_1_2P.fastq.gz | \
    samtools view -S -b - | \
    samtools sort - results/bam/${base}.aligned.sorted"
index="samtools index results/bam/${base}.aligned.sorted.bam"
stats="samtools flagstat results/bam/${base}.aligned.sorted.bam > results/bam/${base}.aligned.sorted.stats.txt"
cat_stats="cat results/bam/*stats.txt > docs/alignment_stats.txt"
counts="samtools idxstats results/bam/${base}.aligned.sorted.bam > results/counts/${base}.aligned.sorted.counts.txt"

# print the $base of the sample to keep track
echo $base

# print commands in order and run them
# mapping to sorted BAM
echo $map
eval $map

# indexing BAM
echo $index
eval $index

# calculate BAM stats
echo $stats
eval $stats
echo $cat_stats
eval $cat_stats

# infer expression counts
echo $counts
eval $counts
```

Finally, we can run this script in a loop to perform the expression workflow on all of our samples.

```bash
$ pwd

$ for sample in data/trimmed_fastq/*1_1P.fastq.gz
do bash docs/expression_pipeline.sh $sample data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.nodups.fasta.gz
done
```

You should conclude by copying the resuling `*.counts.txt` files to your personal computer, saving them to a directory called `dc_workshop` on your Desktop.
