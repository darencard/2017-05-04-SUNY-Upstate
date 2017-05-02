---
layout: lesson
title: Planning for DC projects
date: 2017-05-04
---

## Learning Objectives 
* Think about and understand the types of metadata a sequencing experiment will generate
* Make decisions about how (if) data will be stored, archived, shared, etc. 
* Anticipate strategies we will need to learn in the rest of the lesson set.
* Introduction to fastq file.


*planning, documenting, and organizing* are the key to good reproducible science. 

## Discussion Questions

Before we go any further here are some important questions to consider. If you are learning at a workshop, please discuss these questions with your neighbor, and your instructors will collect your answers (on minute cards or in the Etherpad).

**Metadata** 

1. What is your definition of metadata?
2. What kinds of metadata might a sequencing project generate?

**Working with sequence data**

1. What challenges do you think you'll face (or have already faced) in working with a large sequence dataset?
2. What is your strategy for saving and sharing your sequence files?
3. How can you be sure that your raw data have not been unintentionally corrupted?
4. Where/how will you (did you) analyze your data - what software, what computer(s)?

## Exercises

**A. Examining a fastq file**

Although it looks complicated (and maybe it is), its easy to understand the fastq format with a little decoding. Some rules about the format include…

Line	Description
1	Always begins with ‘@’ and then information about the read
2	The actual DNA sequence
3	Always begins with a ‘+’ and sometimes the same info in line 1
4	Has a string of characters which represent the quality scores; must have same number of characters as line 2

so for example in our data set, one complete read is: 

```
head -n4 ~/dc_sample_data/untrimmed_fastq/SRR098026.fastq
```

    @SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
    NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
    +SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
    !!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!

**B. Work with read data**

First, copy sample data we will be using for the workshop.

```
scp data.zip dcuser@ip_address
```

Then we need to uncompress the zipped data.

```
gunzip data.zip
```

Now we can list and check the content of the data directory

```
cd data
pwd
ls
cd untrimmed_fastq
ls
gunzip SRR2040575_brain_1.fastq.gz
```

Similarly to how we examine the fastq file in the sample dataset, we now checkout the fastq file SRR2040575_brain_1.fastq

```
head SRR2040575_brain_1.fastq
```

You will see:

    @SRR2040575.1 1/1
    GGTCATCGTGTGATTGGCTCAGGCCAGCTTCCTTCTGCCTAGGCCTCCGGCGGAAGATGAGCTTAGTGCCGCTGCGCAGGAAACTCACTTTGTGATCCTT
    +
    BC@FFFFFHHCCFHIIJIJIIJJGHHGIIJGHIIGJIGIJIGIIIJIJJJIJFFCDBDDCDCCDCDCD@CDDDDD>BBDDBDDACACCDDDDDDDD@@CC
    @SRR2040575.2 2/1
    GTCCGCTTCTGTTAATTTATTCTTTTTAATAAACTGGGTGTTGAGATACCTATATAAGCAGTCCATATAGTCTGCACCCTTGCTGTATTCTTCCCAGTAC
    +
    BB@FFFFFHHHHHIGHJJIJJJJJJJJJGJJGJJJJIJHHIJIIHGIJJIIIIIIIIJJIJJJHJGIIIJIJJIJFHGHFFFFFEC@CEEEEDDCDDDCD
    
Note this file might look slightly different from the sample file. But both are in valid FASTQ format. 

**C. More on the data **

**Excercise**

1. Explore the data directory. What do we have?
2. What do the trimmed/untrimmed mean in the directory names?
3. What does the '.gz' extension on the filenames indicate?
4. What is the total file size - what challenges in downloading and sharing these data might exist?

## Summary 

Before analysis of data has begun, there are already many potential areas for errors and omissions. Keeping organized and keeping a critical eye can help catch mistakes. 

One of Data Carpentry's goals is to help you achieve *competency* in working with bioinformatics. This means that you can accomplish routine tasks, under normal conditions, in an acceptable amount of time. While an expert might be able to get to a solution on instinct alone - taking your time, using Google, and asking for help are all valid ways of solving your problems. As you complete the lessons you'll be able to use all of those methods more efficiently.

## Where to go from here?

What are the minimum metadata standards for your experiment/datatype - 

[Biosharing.org's listing of minimum information standards](https://biosharing.org/standards/?selected_facets=isMIBBI:true&selected_facets=domains_exact:DNA%20sequence%20data)

More reading about '[FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).


