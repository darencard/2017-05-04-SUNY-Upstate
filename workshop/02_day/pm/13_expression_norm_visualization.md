---
layout: topic
title: Intro to data frames
author: Data Carpentry contributors
minutes: 30
---

------------

> ## Learning Objectives
>
> * Understand *tidy* data and how to create it.
> * Learn how to use `dplyr` to subset and summarize data frames.
> * Normalize expression and infer differential expression.
> * Use `ggplot2` to create plots commonly used for expression analysis.

------------

# *Tidy* data

So far, we've been working with our dataset in the form of a matrix, where samples and metadata
are columns and observations (i.e., transcripts) are rows. This is a logical way to structure data and
is intutive to us. However, this data structure isn't as easy to manipulate (take our word for it for now)
and tasks like subsetting, summarizing, and plotting, while still possible, are less efficient. The best way
to organize data is using a structure we call *tidy*, where variables are stored as columna and observations
are stored as rows. The best way to understand this is to make our dataset *tidy* (it currently isn't).

The *tidy* paradigm has emerged and expanded in recent years, and as a result many of the commonly-used tools
that create and manipulate *tidy* data have been combined together into a single reference, called the 
*tidyverse*. Let's start by installing this resource (called a `package` in the R world). Installing in R is 
quick and easy, which is why many people like it so much, and can be accomplished with a single command. Once
we've installed a package, we can load it into our *Environment* so it is available for us to use.

```{r, eval=TRUE,  purl=FALSE}
# install.packages("tidyverse")
library(tidyverse)
```

If your data isn't currently loaded, please do so using the same command as before.

```{r, eval=TRUE,  purl=FALSE}
raw_exp <- read.table("full_expression_counts.raw.txt", header=TRUE, sep="\t")
```

Now let's look at the structure of our current dataset.

```{r, eval=TRUE,  purl=FALSE}
head(raw_exp)
```

Note the characteristics of the current `data.frame`. For one thing is is pretty wide, meaning it has several
columns, especially when you consider it is only really recording one variable, expression. 'Wide' data is
usually a pretty good indicator that a dataset isn't tidy.

We're going to use a command from `tidyr` a subset of the `tidyverse` package we just downloaded. We're going to 
manipulate our dataset using `dplyr` another subset of `tidyverse`, which allows us to chain commands together
like we did when we used pipes (`|`) in Bash. Instead of a pipe, the syntax is `ddplyr` is `%>%`, and as we would
expect, our output is used as input (normally the first argument to a function) for the next command. Let's 
execute the following command to make our dataset tidy.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts <- raw_exp %>% 
  gather(sample, expression, -transcript, -length)
```

Since there are a lot of new concepts here, let's break this command down:
1. We are setting whatever the output of the command is to a variable called `tidy_counts`
2. We use the `%>%` pipe operator to pass the output of the first command as the input
   of the second command.
   a. `raw_exp` just prints the untidy `data.frame`
   b. `gather` does the tidying. This command looks across the columns in the `data.frame`
       and gathers them such that the header line describes some variable of an observation
       and each of the rows is an individual observation. We provide the argument `sample`
       to name the variables encompassed by the header (i.e., the samples) and `expression` to
       name the variable corresponding to the expression values. The `-transcript` and `-length`
       tell `gather` to ignore those two fields from the original dataset and treat them as separate
       variables describing each observation. If this doesn't make much sense, hopefully the output
       will make it more obvious.
       
```{r, eval=TRUE,  purl=FALSE}
head(tidy_counts)
str(tidy_counts)
```

Notice how we now have only 4 columns. `transcript` and `length` have stayed the same, since we told
gather to exclude them. `transcript` is the unit of measurement here, and `length` is a variable of the
transcript length. The newly formated columns are `sample` and `expression`. `sample` now represents the 
variable reflecting the samples we are analyzing. You'll see that it repeats as the same value for each
transcript before repeating when a new sample is encountered. Finally, for each combination of variables
in the other 3 columns, you have the `expression` variable that provides the expression under that exact
circumstance.

This is *tidy* data and you'll notice that we only have variables as columns and observations as rows.
As you have probably noticed, tidy data tends to be in a 'long' format, as the data that was 
side-by-side has essentially been stacked. Great, our data is now tidy!

Let's demo the package `dplyr` some more. You've already seen the pipe operator `%>%` but `dplyr` is
also designed for manipulating tidy datasets. `dplyr` makes it easy to subset.

```{r, eval=TRUE,  purl=FALSE}
# extract the first sample and show the head
tidy_counts %>%
    subset(sample=="SRR2040575_brain") %>%
    head()
```

Now you might be starting to see why tidy data is so nice. We can also summarize the data in many ways.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
    subset(sample=="SRR2040575_brain") %>%
    summarize(mean=mean(expression))
```

This command is pretty easy, it simply outputs the mean expression of the target sample. But what if
we want the mean for each sample? We could write a `for` loop, but `dplyr` already has that covered.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
    group_by(sample) %>%
    summarize(mean=mean(expression))
```

In this command, the `group_by` command essentially bins the samples, and `summarize` operates on each bin.

> ### Challenge
>
> Based on some of the commands you've used, can you answer the following questions?
>
> * How many genes have an expression level greater than 10 in each sample?
> * What is the mean and standard deviation of expression in each sample?

```{r, echo=FALSE, purl=TRUE}

## Challenge
## Based on some of the commands you've learned, can you answer the following questions?
## * How many genes have an expression level greater than 10 in each sample?
## * What is the mean and standard deviation of expression in each sample?

```

<!---
```{r, echo=FALSE, purl=FALSE}
## Answers
##
## * tidy_counts %>% subset(expression>10) %>% group_by(sample) %>% count()
## * tidy_counts %>% group_by(sample) %>% summarize(mean=mean(expression), stdev=sd(expression))
```
--->

# Plotting with `ggplot`

Okay, now let's do some plotting, one of the most powerful functions in R. Today we are going to do 
our plotting with the package `ggplot2`, which makes plot construction intuitive and produces pretty
nice looking default plots. `ggplot2` is part of the `tidyverse`.

To begin and cover the basics, we'll do a simple scatterplot. We only really have 2 variables we can
compare, transcript length and expression, so let's do that.

```{r, eval=TRUE,  purl=FALSE}
ggplot(data=tidy_counts) +
  geom_point(aes(x=length, y=expression))
```

This command produces a scatter plot with transcript length on the x-axis and expression on the y-axis.
The first line calls the ggplot function and provides it with the `data.frame` we want to use. The second
line, separated by the `+`, layers a `geom` on top. `geom` is a `ggplot` construct and is referring to the
type of plot to create, in this case a `point` plot or scatterplot. Within `geom_point` we provide the
argument `aes`, which is short for aesthetic. Within `aes` we are mapping a given data variable to a visual
property. In this case, we are mapping `length` to the visual property `x`, which is the x-axis. `aes` can
also be put in the initial `ggplot` command, as these mappings are conserved throughout a plot.

With `ggplot` being in the `tidyverse`, we can also use it within `%>%` pipes.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  ggplot() +
  geom_point(aes(x=length, y=expression))
```
You've probably noticed that we are plotting all observations, not taking any information about what sample
they came from into account. It is easy to do this in `ggplot`.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  ggplot() +
    geom_point(aes(x=length, y=expression, color=sample))
```

We can also easy separate each sample into subplots.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  ggplot() +
    geom_point(aes(x=length, y=expression)) + 
    facet_wrap(~sample)
```

This is great, but overall we are noticing that transcript length vs. expression isn't all that interesting.
We're most interested in comparing expression between samples. Let's create a boxplot instead that shows
the distribution of expression for each sample.

```{r, eval=TRUE,  purl=FALSE}
ggplot(tidy_counts) +
  geom_boxplot(aes(x=sample, y=expression))
```

> ### Challenge
>
> The boxplots leave a lot to be desired. The aesthetics are easy to fix. The bigger problem
> is that our plot is dominated by the outliers. Construct a command that allows us to better
> look at the distribution of expression by only looking at transcripts with an expression of
> less than 100.
>
> What do you conclude based on this plot?

```{r, echo=FALSE, purl=TRUE}

## Challenge
## The boxplots leave a lot to be desired. The aesthetics are easy to fix. The bigger problem
## is that our plot is dominated by the outliers. Construct a command that allows us to better
## look at the distribution of expression by only looking at transcripts with an expression of
## less than 100.
##
## What do you conclude based on this plot?

```

<!---
```{r, echo=FALSE, purl=FALSE}
## Answers
##
## * tidy_counts %>% filter(expression < 100) %>% ggplot() + geom_boxplot(aes(x=sample, y=expression))
## * expression is uneven, and is likely a product of unequal sequence depth
```
--->

# Expression normalization, analysis, and visualization

The boxplot you just generated hopefully made you realize that expression levels are often variable
across samples, and this is due to uneven amounts of sequencing and other factors. Therefore, comparing
expression levels based on the raw expression data is like comparing apples to oranges. We get around this
problem using *normalization*, where we essentially try to adjust each sample such that their overall expression
profiles are more similar, and thus the samples are more comparable.

There are many different statistical methods for normalizing data and several statistical packages (many in R) that
do it. We'll use one popular one called `DESeq2`. We won't talk about the specifics of the normalization method 
and if you find yourself doing RNAseq work, you should research the normalization method that is best for you.

`DESeq2` is available from `Bioconductor`, so the installation process is a little different, but still very easy.
You can find more information about installation and the `DESeq2` package at [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

```{r, eval=TRUE,  purl=FALSE}
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
```

We're going to be using our original, untidy dataset for `DESeq2`, as it fits the format the package requires much better.
Transcript length has no bearing on normalization, so we will first get rid of that column. Let's also import another
`data.frame` that provides metadata on our samples, namely on the organs that samples came from.

```{r, eval=TRUE,  purl=FALSE}
raw_exp$length <- NULL
metadata <- read.table("full_expression_counts.metadata.txt", header=TRUE, sep="\t")
```

`DESeq2` is very easy to use. We begin by loading our data into a `DESeq2` dataset (`dds`), which is a special object
unique to `DESeq2` that stores the data and other information in one nice structure. Objects in `DESeq2` and other R packages
are the reason why we consider R and *object-oriented* programming language. We'll use a command to load our data, specifying
our `metadata` dataset to provide experimental context and telling `DESeq2` that we want to use the assignments in the `organ`
field of `metadata` as the data levels. There are 4 data levels here, the individual organs, but you could also do a
treatment vs. control experiment.

```{r, eval=TRUE,  purl=FALSE}
dds <- DESeqDataSetFromMatrix(countData=raw_exp, 
                              colData=metadata, 
                              design=~organ, 
                              tidy=TRUE)
```

Once we have our data loaded, it is very easy to have `DESeq2` perform expression normalization. This is
what the following command does. It also does other analysis, namely inferring differential expression across
the data levels in our dataset. Again, there are several ways of inferring differential expression, with the
two most common being linear models and Fisher's exact tests. We're bypass the details of this with `DESeq2`,
but you should definitely research this if you plan to do any differential expression analyses.

```{r, eval=TRUE,  purl=FALSE}
dds <- DESeq(dds)
```

The results of all of these analyses are automatically stored in the same `DESeq2` object. We can get the 
results, including information about the level of expression change and the statistical difference.

```{r, eval=TRUE,  purl=FALSE}
res <- results(dds, tidy=TRUE)
```

You'll see that the output is just a simple `data.frame`. We can also use `dplyr` to figure out how
many genes have significantly different expression levels across our samples.

```{r, eval=TRUE,  purl=FALSE}
head(res)
res %>% 
  filter(padj < 0.05) %>%
  summarise(count=n())
```

The other important data are the normalized expression levels, which are useful for expression visualization.

```{r, eval=TRUE,  purl=FALSE}
norm <- counts(dds, normalized=TRUE)
```

The normalized expression dataset unfortunately isn't tidy, so let's make it tidy before moving forward.

```{r, eval=TRUE,  purl=FALSE}
tidy_norm <- norm %>% 
  as.data.frame() %>%
  cbind(transcripts = rownames(norm)) %>%
  gather(sample, expression, -transcripts)
```

> ### Challenge
>
> Create a boxplot of the normalized expression per sample. What can you conclude by comparing the pre-normalization
> and post-normalization datasets?

```{r, echo=FALSE, purl=TRUE}

## Challenge
## Create a boxplot of the normalized expression per sample. What can you conclude by comparing the pre-normalization
## and post-normalization datasets? You may have to restrict the expression values to really see the core of the
## distribution.

```

<!---
```{r, echo=FALSE, purl=FALSE}
## Answers
##
## * tidy_norm <- norm %>% filter(expression < 100) %>% ggplot() + geom_boxplot(aes(x=sample, y=expression))
## * Normalization made the expression distributions much more similar and comparable.
```
--->

Okay, finally to the best part, producing high-quality figures that you can use to show off your work. One
common analysis you see in gene expression studies is a PCA, which allows you to visually cluster samples
based on expression profile. We can perform an analogous analysis using a similar metric right in `DESeq2`.
The function we will use is `vst` and you can view the help page to get more information.

```{r, eval=TRUE,  purl=FALSE}
vsdata <- vst(dds)
```

Once we have performed this analysis, `DESeq2` also has a nice function (`plotPCA`) to plot the result. This 
function is what we call a *wrapper*, because it uses `ggplot2` for plotting the results. Let's produce this 
plot and save it to a variable `pca`. Since `pca` is now a `ggplo2` object, we can customize it. In this case
we will add a theme.

```{r, eval=TRUE,  purl=FALSE}
pca <- plotPCA(vsdata, intgroup="organ")
pca + theme_bw()
```

```{r, eval=TRUE,  purl=FALSE}
scripts <- arrange(res, padj)$row[1:50]
```

```{r, eval=TRUE,  purl=FALSE}
tidy_norm <- norm %>% 
  t() %>%
  scale() %>%
  t() %>% 
  as.data.frame() %>%
  cbind(transcripts = rownames(norm)) %>%
  filter(transcripts %in% scripts) %>%
  gather(sample, expression, -transcripts)
```

```{r, eval=TRUE,  purl=FALSE}
heatmap <- ggplot(tidy_norm) +
  geom_tile(aes(x=sample, y=transcripts, fill=expression))

heatmap
```

```{r, eval=TRUE,  purl=FALSE}
# install.packages("viridis")
library(viridis)

heatmap +
  scale_fill_viridis()
```

```{r, eval=TRUE,  purl=FALSE}
volcano <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj),
                           text=row)) +
  geom_point(aes(color=factor(res$padj <= 0.01 & abs(res$log2FoldChange) >= 5))) + 
  theme_bw()

volcano
```

```{r, eval=TRUE,  purl=FALSE}
# install.packages("plotly")
library(plotly)

ggplotly(volcano)
```
