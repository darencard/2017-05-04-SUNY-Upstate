---
title: "Variant calling workflow"
author: "Mary Piper, Meeta Mistry"
date: "Thursday, March 31, 2016"
---

Approximate time: 120 minutes

## Learning Objectives:

* Understand the different steps involved in variant calling
* Use a series of command line tools to execute a variant calling workflow
* Becoming familiar with data formats encountered during variant calling


## Setting up

To get started with this lesson, make sure you are in `dc_workshop/data`. We've ignored it so far, but there is also a reference transcriptome we will be using in the `reference` directory.

Your directory structure should now look like this:

<pre>
dc_workshop
├── data
    ├── ref_genome
        └── GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.fasta.gz
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
$ mkdir results/sam results/bam
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
$ bwa index data/reference/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.fasta.gz     # This step helps with the speed of alignment
```

Eventually we will loop over all of our files to run this workflow on all of our samples, but for now we're going to work with the quality filtered paired sequences only for just one sample in our dataset `SRR2040575_brain_1_1P.fastq.gz` and `SRR2040575_brain_1_2P.fastq.gz`:

```bash
$ ls -alh ~/dc_workshop/data/trimmed_fastq/SRR2040575_brain_1_[12]P.fastq.gz 
```

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an aligner. BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate, so we will be using that. 
    
Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using and try to understand the options.*

```bash
$ bwa mem data/ref_genome/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.fasta.gz \
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

### Sort BAM file by coordinates

Sort the BAM file:

```bash
$ samtools sort results/bam/SRR2040575_brain.aligned.bam results/bam/SRR2040575_brain.aligned.sorted
```

*SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.*

## Using pipling to reduce input/output

One awesome feature of *BWA* and *samtools* is that they allow piping, so we can easily string commands together to more efficiently perform our mapping. This reduces the amount of input/output of data to our hard drive, making the process significantly faster. Moreover, it helps save disk space, which can rapidly disappear with the large file sizes produced in genomic data. This is just like the piping commands that you learned in the introduction to Bash, making these software very modular and powerful.

Therefore, instead of mapping the `*.fastq` files to produce `*.sam` files, converting them to binary `*.bam` files, and then sorting to produce a final `*.bam` file, we will be able to just go from the `*.fastq` files to the sorted `*.bam` files using one command! With this in mind, let's just get rid of our directory for `*.sam` files.

```bash
$ # removed recursively
$ rm -r results/sam
```

Now we can run out single mapping command by stringing the commands from above together with some special syntax.

```bash
$ bwa mem data/ref_genome/GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.fasta.gz \
    data/trimmed_fastq/SRR2040575_brain_1_1P.fastq.gz data/trimmed_fastq/SRR2040575_brain_1_2P.fastq.gz | \
    samtools view -S -b - | \
    samtools sort - results/bam/SRR2040575_brain.aligned.sorted
```

Notice a couple things about this command:
1. The pipe (`|`) symbols between the commands, which basically tell the computer to use the output from the first command (called STDOUT) as the input for the second command (called STDIN), without writing it to the disk.
2. The dashes (`-`) where input files names would normally be, which tells the respective command to look for input data from STDIN. And no output file names are specified except for the last command, because they are not needed until that point.
    
