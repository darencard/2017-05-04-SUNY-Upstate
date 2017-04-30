---
layout: lesson
root: .
title: Introduction to workshop dataset
minutes: 30
---

## Learning Objectives
* Understand how and why we choose this dataset
* Have a general idea of the experiment and its objectives
* Get a general understanding of the NCBI SRA
* Understand how to download a summary run table of the SRA
* Be able to open a run table using a spreadsheet editor

## Introduction to the dataset

#### Description of dataset

RNAseq has emerged as an important technology for understanding how gene expression modulates organism phenotypes. The most common form of RNAseq targets mature mRNAs and therefore provides information about gene expression that usually underlies protein formation (though not all mRNAs become proteins). There are other forms of RNAseq that focus on other RNA molecules, like small RNA sequencing, but we will focus on mRNAseq.

Despite having the same underlying genome, cells in a single organisms body can vary radically in their phenotypes, and thus their functions in a multicellular organism. In humans, for instance, there are very different types of cells in our distinct tissues and organs, and even within organs there are different cell types that perform unique functions. Differential regulation of mRNA is the primary mediator explaining how the same genome can create very different phenotypes. Therefore, if one profiles the expression of different organs, he or she would expect to see very different levels of gene expression across genes. This is the basis for this workshop, where we will be using mRNAseq data from 9 human samples, comprising 4 tissues, to understand the general expression patterns that vary between tissues and the specific genes that are most different in expression between tissues.

| SRA Run Number | Organ |
| -------------- | ----- |
| SRR2040575 | Brain |
| SRR2040576 | Brain |
| SRR2040577 | Heart |
| SRR2040578 | Heart |
| SRR2040579 | Liver |
| SRR2040580 | Liver |
| SRR2040581 | Testes |
| SRR2040582 | Testes |
| SRR2040583 | Testes |


We want to be able to quantify gene expression across these samples and determine how different expression profiles are between organs. Ultimately, we will answer the questions:

- What is the gene expression profile across each organ?
- How different is gene expression across genes between organs?

#### Accessing the original archived data
The sequencing dataset was attained from the [NCBI Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra), which is a large (>3 quadrillion basepairs as of 2014) repository for next-generation sequence data. Like many NCBI databases, it is complex and mastering its use is greater than the scope of this lesson. Very often, as in this study, there will be a direct link (perhaps in the supplemental information) to where on the SRA the dataset can be found. The link from the paper is: [http://www.ncbi.nlm.nih.gov/sra?term=SRP058740](http://www.ncbi.nlm.nih.gov/sra?term=SRP058740)

###### A. Locate the Run Accessory Table for the Dataset on the SRA

1. Access the dataset from the provided link: [http://www.ncbi.nlm.nih.gov/sra?term=SRP058740](http://www.ncbi.nlm.nih.gov/sra?term=SRP058740).  
You will be presented with a page for the overall SRA accession SRP058740 - this is a collection of all the experimental data
2. This dataset includes expression from humans, which we will be using, and other primates. At the top right you will see the `Results by taxon`, where you can click on `Homo sapiens` to pare down the dataset to only the human data of interest.
3. Click on the first entry ([GSM1695913](https://www.ncbi.nlm.nih.gov/sra/SRX1038913[accn])); this will take you to a page for an SRX (Sequence Read eXperiment). Take a few minutes to examine some of the descriptions on the page
4. Click on the ['All runs'](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP058740) link under where it says **Study**; this is a description of all of the NGS datasets related to the experiment.
5. Go to the top of the page and in the **Total** row you will see there are 27 runs, 234.83 Gb data, and 376.77 Gbases of data. Click the 'RunInfo Table' button.

We are not downloading any actual sequence data here! This is only a text file that fully describes the entire dataset

You should now have a file called **SraRunTable.txt**

###### B. Review the SraRunTable in a spreadsheet program

1. Using your choice of spreadsheet program open the **SraRunTable.txt** file. If prompted this is a tab-delimited file.

###### C. Download and prepare the dataset

1. We've prepared this dataset for you and you should download it to be prepared for the rest of the workshop. The dataset has been archived at [figshare](https://figshare.com/). Just visit the figshare page for this dataset at [https://figshare.com/s/ec43b22e2054b2557498](https://figshare.com/s/ec43b22e2054b2557498) and click `Download`.
2. The dataset is large so it may take some time to download The dataset is stored as a `.zip` archive, so be sure that your computer does not automatically unzip it. Once it is downloaded, you should transfer it to your remote cloud computer and place it in a newly created directory in $HOME called `data_analysis`. Once the `.zip` file is there, you can unzip it using `unzip`. Use `man unzip` to determine how to do this.
3. Once unzipped, you will see the file `GCF_000001405.36_GRCh38.p10_rna.known_refseq_proteincoding.fasta.gz`, which is the reference transcriptome we will be using. There are also directories named `raw`, which contains the raw sequencing data we will be starting with, and `backup`, which includes several pre-prepared files for various stages of our analysis pipeline. The data within `backup` should be ignored and is only provided for reference and to help keep learners on track if problems arise.

Note: The downloaded dataset is not the full dataset you explored on the SRA and we've taken steps to reduce the data file size (considerably) so that you can efficiently analyze these data during the workshop.

***
**Exercises**

Discuss with the person next to you:

1. Besides humans, what other organisms are included in the dataset?
2. Looking at some of the other metadata, is there anything about the way the data were collected that may confound our ability to generalize the patterns we resolve?
3. What was the sequencing platform used for this experiment?
4. What samples in the experiment contain [paired end](http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html) sequencing data?
5. What other kind of data is available?
6. Are you collecting this kind of information about your sequencing runs?
***
After answering the questions, you should avoid saving this file; we don't want to make any changes. If you were to save this file, make sure you save it as a plain **.txt** file.


## Where to learn more

#### About the Sequence Read Archive

* You can learn more about the SRA by reading the [SRA Documentation](http://www.ncbi.nlm.nih.gov/Traces/sra/)
* The best way to transfer a large SRA dataset is by using the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc)

#### References

Ruiz-Orera, J., Hernandez-Rodriguez, J., Chiva, C., Sabido, E., Kondova, I., Bontrop, R., Marques-Bonet, T., Alba, M.M. Origins of de novo genes in human and chimpanzee (2015) PLOS Genetics, 11 (12), e1005721.
[Paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005721)

Data on NCBI SRA: [http://www.ncbi.nlm.nih.gov/sra?term=SRP058740](http://www.ncbi.nlm.nih.gov/sra?term=SRP058740)
