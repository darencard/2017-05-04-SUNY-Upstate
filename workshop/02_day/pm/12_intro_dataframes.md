---
layout: topic
title: Starting with data
author: Data Carpentry contributors
minutes: 20
---

------------

> ## Learning Objectives
>
> * Load external tabular data from a text file into R.
> * Describe what an R data frame is.
> * Summarize the contents of a data frame in R.
> * Perform basic subsetting of a data frame in R.
> * Understand categorical data in R using factors.

------------

The data file required for this lesson can be [downloaded directly here](https://raw.githubusercontent.com/darencard/2017-05-04-SUNY-Upstate/gh-pages/workshop/data/full_expression_counts.raw.txt) or [viewed in Github](../../data/full_expression_counts.raw.txt).

We are going to use the R function `read.table()` to load the data file into memory (as a `data.frame`):

```{r, eval=TRUE,  purl=FALSE}
raw_exp <- read.table("full_expression_counts.raw.txt", header=TRUE, sep="\t")
```

This statement doesn't produce any output because assignment doesn't display
anything. If we want to check that our data has been loaded, we can print the
variable's value: `exp_data`

Alternatively, wrapping an assignment in parentheses will perform the assignment
and display it at the same time.

```{r, eval = TRUE, purl = FALSE}
(raw_exp <- read.table("full_expression_counts.raw.txt", header=TRUE, sep="\t"))
```

Wow... that was a lot of output. At least it means the data loaded properly. Let's check the top (the first 6 lines) of this `data.frame` using the function `head()`:

```{r, results='show', purl=FALSE}
head(raw_exp)
```

We've just done two very useful things.
1. We've read our data in to R, so now we can work with it in R
2. We've created a data frame (with the `read.table` command) the 
standard way R works with data. 

# What are data frames?

`data.frame` is the _de facto_ data structure for most tabular data and what we
use for statistics and plotting.

A `data.frame` is a collection of vectors of identical lengths. Each vector
represents a column, and each vector can be of a different data type (e.g.,
characters, integers, factors). The `str()` function is useful to inspect the
data types of the columns.

A `data.frame` can be created by the functions `read.csv()` or `read.table()`, in
other words, when importing spreadsheets from your hard drive (or the web).

By default, `data.frame` converts (= coerces) columns that contain characters
(i.e., text) into the `factor` data type. Depending on what you want to do with
the data, you may want to keep these columns as `character`. To do so,
`read.csv()` and `read.table()` have an argument called `stringsAsFactors` which
can be set to `FALSE`:

Let's now check the __str__ucture of this `data.frame` in more details with the
function `str()`:

```{r, purl=FALSE}
str(raw_exp)
```

# Inspecting `data.frame` objects

We already saw how the functions `head()` and `str()` can be useful to check the
content and the structure of a `data.frame`. Here is a non-exhaustive list of
functions to get a sense of the content/structure of the data.

* Size:
    * `dim()` - returns a vector with the number of rows in the first element, and
    the number of columns as the second element (the __dim__ensions of the object)
    * `nrow()` - returns the number of rows
    * `ncol()` - returns the number of columns
* Content:
    * `head()` - shows the first 6 rows
    * `tail()` - shows the last 6 rows
* Names:
    * `names()` - returns the column names (synonym of `colnames()` for `data.frame`
	objects)
   * `rownames()` - returns the row names
* Summary:
   * `str()` - structure of the object and information about the class, length and
	content of  each column
   * `summary()` - summary statistics for each column

Note: most of these functions are "generic", they can be used on other types of
objects besides `data.frame`.

> ### Challenge
>
> Based on some of the commands you've used, can you answer the following questions?
>
> * What is the class of the object `raw_exp`?
> * How many rows and how many columns are in this object?
> * How many transcripts have been evaluated?
> * What is the average transcript length?

```{r, echo=FALSE, purl=TRUE}

## Challenge
## Based on some of the commands you've learned, can you answer the following questions?
## * What is the class of the object `raw_exp`?
## * How many rows and how many columns are in this object?
## * How many transcripts have been evaluated?
## * What is the average transcript length?

```

<!---
```{r, echo=FALSE, purl=FALSE}
## Answers
##
## * class: data frame
## * how many rows: 38681,  how many columns: 11
## * how many transcripts: 38681
## * average length: 3392
```
--->


## Indexing and subsetting data frames

Our survey data frame has rows and columns (it has 2 dimensions), if we want to
extract some specific data from it, we need to specify the "coordinates" we
want from it. Row numbers come first, followed by column numbers. However, note
that different ways of specifying these coordinates lead to results with
different classes.

```{r, purl=FALSE}
raw_exp[1]      # first column in the data frame (as a data.frame)
raw_exp[, 1]    # first column in the data frame (as a vector)
raw_exp[1, 1]   # first element in the first column of the data frame (as a vector)
raw_exp[1, 6]   # first element in the 6th column (as a vector)
raw_exp[1:3, 7] # first three elements in the 7th column (as a vector)
raw_exp[3, ]    # the 3rd element for all columns (as a data.frame)
head_raw_exp <- raw_exp[1:6, ] # equivalent to head(surveys)
```

`:` is a special function that creates numeric vectors of integers in increasing
or decreasing order, test `1:10` and `10:1` for instance.

You can also exclude certain parts of a data frame using the "`-`" sign:

```{r, purl=FALSE}
raw_exp[,-1]          # The whole data frame, except the first column
raw_exp[-c(7:34786),] # Equivalent to head(surveys)
```

As well as using numeric values to subset a `data.frame` (or `matrix`), columns
can be called by name, using one of the four following notations:

```{r, eval = FALSE, purl=FALSE}
raw_exp["transcript"]       # Result is a data.frame
raw_exp[, "transcript"]     # Result is a vector
raw_exp[["transcript"]]     # Result is a vector
raw_exp$transcript          # Result is a vector
```

For our purposes, the last three notations are equivalent. RStudio knows about
the columns in your data frame, so you can take advantage of the autocompletion
feature to get the full and correct column name.

> ### Challenge
>
> 1. Create a `data.frame` (`raw_exp_200`) containing only the observations from
>    row 200 of the `raw_exp` dataset.
>
> 2. Notice how `nrow()` gave you the number of rows in a `data.frame`?
>
>      * Use that number to pull out just that last row in the data frame.
>      * Compare that with what you see as the last row using `tail()` to make
>        sure it's meeting expectations.
>      * Pull out that last row using `nrow()` instead of the row number.
>      * Create a new data frame object (`raw_exp_last`) from that last row.
>
> 3. Use `nrow()` to extract the row that is in the middle of the data
>    frame. Store the content of this row in an object named `raw_exp_middle`.
>
> 4. Combine `nrow()` with the `-` notation above to reproduce the behavior of
>    `head(surveys)` keeping just the first through 6th rows of the surveys
>    dataset.
>
> 5. Pick one of the samples in this `data.frame` and calculate the mean expression
>    of only that column.


```{r, echo=FALSE, purl=TRUE}
### Challenges:
###
### 1. Create a `data.frame` (`raw_exp_200`) containing only the
###    observations from row 200 of the `raw_exp` dataset.
###
### 2. Notice how `nrow()` gave you the number of rows in a `data.frame`?
###
###      * Use that number to pull out just that last row in the data frame
###      * Compare that with what you see as the last row using `tail()` to make
###        sure it's meeting expectations.
###      * Pull out that last row using `nrow()` instead of the row number
###      * Create a new data frame object (`raw_exp_last`) from that last row
###
### 3. Use `nrow()` to extract the row that is in the middle of the
###    data frame. Store the content of this row in an object named
###    `surveys_middle`.
###
### 4. Combine `nrow()` with the `-` notation above to reproduce the behavior of
###    `head(surveys)` keeping just the first through 6th rows of the surveys
###    dataset.
### 5. Pick one of the samples in this `data.frame` and calculate the mean expression
###    of only that column.

```

<!---
```{r, purl=FALSE}
## Answers
raw_exp_200 <- raw_exp[200, ]
raw_exp_last <- raw_exp[nrow(raw_exp), ]
raw_exp_middle <- raw_exp[nrow(raw_exp)/2, ]
raw_exp_head <- raw_exp[-c(7:nrow(raw_exp)),]
mean(raw_exp[,3])
```
--->


## Factors

Factors are used to represent categorical data. Factors can be ordered or
unordered and are an important class for statistical analysis and for plotting.

Factors are stored as integers, and have labels associated with these unique
integers. While factors look (and often behave) like character vectors, they are
actually integers under the hood, and you need to be careful when treating them
like strings.

In the data frame we just imported, let's do 
```{r, purl=TRUE}
str(raw_exp)
```

We can see the names of the multiple columns. And, we see that 
some say things like `Factor w/ 38681`

When we read in a file, any column that contains text is automatically
assumed to be a factor. Once created, factors can only contain a pre-defined set values, known as
*levels*. By default, R always sorts *levels* in alphabetical order, but one can rearrange the if
desired.

### Using `stringsAsFactors=FALSE`

By default, when building or importing a data frame, the columns that contain
characters (i.e., text) are coerced (=converted) into the `factor` data
type. Depending on what you want to do with the data, you may want to keep these
columns as `character`. To do so, `read.csv()` and `read.table()` have an
argument called `stringsAsFactors` which can be set to `FALSE`.

In most cases, it's preferable to set `stringsAsFactors = FALSE` when importing
your data, and converting as a factor only the columns that require this data
type.


