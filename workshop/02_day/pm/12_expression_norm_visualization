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
expect, our output is used as input (normally the first argument to a function) for the next command.

```{r, eval=TRUE,  purl=FALSE}
tidy_counts <- raw_exp %>% 
  gather(sample, expression, -transcript, -length) %>%
  as.data.frame()
```

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  ggplot() +
    geom_point(aes(x=length, y=expression))
```

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  ggplot() +
    geom_point(aes(x=length, y=expression, color=sample))
```

```{r, eval=TRUE,  purl=FALSE}
ggplot(tidy_counts) +
  geom_boxplot(aes(x=sample, y=expression))
```

```{r, eval=TRUE,  purl=FALSE}
tidy_counts %>%
  filter(expression < 100) %>%
  ggplot() +
    geom_boxplot(aes(x=sample, y=expression))
```

```{r, eval=TRUE,  purl=FALSE}
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
```

```{r, eval=TRUE,  purl=FALSE}
raw_counts$length <- NULL
metadata <- read.table("~/Desktop/dc/full_expression_counts.metadata.txt", header=TRUE, sep="\t")
```

```{r, eval=TRUE,  purl=FALSE}
dds <- DESeqDataSetFromMatrix(countData=raw_counts, 
                              colData=metadata, 
                              design=~organ, 
                              tidy=TRUE)
```

```{r, eval=TRUE,  purl=FALSE}
dds <- DESeq(dds)
```

```{r, eval=TRUE,  purl=FALSE}
res <- results(dds, tidy=TRUE)
```

```{r, eval=TRUE,  purl=FALSE}
res %>% 
  filter(padj < 0.05) %>%
  summarise(count=n())
```

```{r, eval=TRUE,  purl=FALSE}
norm <- counts(dds, normalized=TRUE)
```

```{r, eval=TRUE,  purl=FALSE}
tidy_norm <- norm %>% 
  as.data.frame() %>%
  cbind(transcripts = rownames(norm)) %>%
  gather(sample, expression, -transcripts)
```

```{r, eval=TRUE,  purl=FALSE}
tidy_norm %>% 
  filter(expression < 100) %>%
  ggplot() +
    geom_boxplot(aes(x=sample, y=expression))
```

```{r, eval=TRUE,  purl=FALSE}
vsdata <- vst(dds)
```

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
