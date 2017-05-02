---
title: Introduction to R
author: Data Carpentry contributors
minutes: 45
---


```{r, echo=FALSE, purl=FALSE, message = FALSE}
source("setup.R")
```

------------

> ### Learning Objectives
>
> * Define the following terms as they relate to R: object, assign, call,
>   function, arguments, options.
> * Create objects and and assign values to them.
> * Use comments to inform script.
> * Do simple arithmetic operations in R using values and objects.
> * Call functions and use arguments to change their default options.
> * Inspect the content of vectors and manipulate their content.

------------


## Creating objects in R

```{r, echo=FALSE, purl=TRUE}
### Creating objects in R
```

You can get output from R simply by typing math in the console:

```{r, purl=FALSE}
3 + 5
12 / 7
```

However, to do useful and interesting things, we need to assign _values_ to
_objects_. To create an object, we need to give it a name followed by the
assignment operator `<-`, and the value we want to give it:

```{r, purl=FALSE}
weight_kg <- 55
```

`<-` is the assignment operator. It assigns values on the right to objects on
the left. So, after executing `x <- 3`, the value of `x` is `3`. The arrow can
be read as 3 **goes into** `x`.  For historical reasons, you can also use `=`
for assignments, but not in every context. Because of
the
[slight](http://blog.revolutionanalytics.com/2008/12/use-equals-or-arrow-for-assignment.html) [differences](https://web.archive.org/web/20130610005305/https://stat.ethz.ch/pipermail/r-help/2009-March/191462.html) in
syntax, it is good practice to always use `<-` for assignments.

In RStudio, typing <kbd>Alt</kbd> + <kbd>-</kbd> (push <kbd>Alt</kbd> at the
same time as the <kbd>-</kbd> key) will write ` <- ` in a single keystroke.


Objects can be given any name such as `x`, `current_temperature`, or
`subject_id`. You want your object names to be explicit and not too long. They
cannot start with a number (`2x` is not valid, but `x2` is). R is case sensitive
(e.g., `weight_kg` is different from `Weight_kg`). There are some names that
cannot be used because they are the names of fundamental functions in R (e.g.,
`if`, `else`, `for`,
see
[here](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Reserved.html)
for a complete list). In general, even if it's allowed, it's best to not use
other function names (e.g., `c`, `T`, `mean`, `data`, `df`, `weights`). If in doubt,
check the help to see if the name is already in use. It's also best to avoid
dots (`.`) within a variable name as in `my.dataset`. There are many functions
in R with dots in their names for historical reasons, but because dots have a
special meaning in R (for methods) and other programming languages, it's best to
avoid them. It is also recommended to use nouns for variable names, and verbs
for function names. It's important to be consistent in the styling of your code
(where you put spaces, how you name variables, etc.). Using a consistent coding
style makes your code clearer to read for your future self and your
collaborators. In R, two popular style guides
are [Hadley Wickham's](http://adv-r.had.co.nz/Style.html)
and [Google's](https://google.github.io/styleguide/Rguide.xml)
(this [one](http://jef.works/R-style-guide/) is also comprehensive).

When assigning a value to an object, R does not print anything. You can force R to print the value by using parentheses or by typing the object name:

```{r, purl=FALSE}
weight_kg <- 55    # doesn't print anything
(weight_kg <- 55)  # but putting parenthesis around the call prints the value of `weight_kg`
weight_kg          # and so does typing the name of the object
```

Now that R has `weight_kg` in memory, we can do arithmetic with it. For
instance, we may want to convert this weight into pounds (weight in pounds is 2.2 times the weight in kg):

```{r, purl=FALSE}
2.2 * weight_kg
```

We can also change a variable's value by assigning it a new one:

```{r, purl=FALSE}
weight_kg <- 57.5
2.2 * weight_kg
```

This means that assigning a value to one variable does not change the values of
other variables.  For example, let's store the animal's weight in pounds in a new
variable, `weight_lb`:

```{r, purl=FALSE}
weight_lb <- 2.2 * weight_kg
```

and then change `weight_kg` to 100.

```{r, purl=FALSE}
weight_kg <- 100
```

What do you think is the current content of the object `weight_lb`? 126.5 or 220?

### Comments

The comment character in R is `#`, anything to the right of a `#` in a script
will be ignored by R. It is useful to leave notes, and explanations in your
scripts.
RStudio makes it easy to comment or uncomment a paragraph: after selecting the
lines you  want to comment, press at the same time on your keyboard
<kbd>Ctrl</kbd> + <kbd>Shift</kbd> + <kbd>C</kbd>.

> ### Challenge
>
> What are the values after each statement in the following?
>
> ```{r, purl=FALSE}
> mass <- 47.5            # mass?
> age  <- 122             # age?
> mass <- mass * 2.0      # mass?
> age  <- age - 20        # age?
> mass_index <- mass/age  # mass_index?
> ```

```{r, echo=FALSE, purl=TRUE}
### Challenge
##
## What are the values after each statement in the following?
##
## mass <- 47.5            # mass?
## age  <- 122             # age?
## mass <- mass * 2.0      # mass?
## age  <- age - 20        # age?
## mass_index <- mass/age  # mass_index?
```

### Functions and their arguments

Functions are "canned scripts" that automate more complicated sets of commands
including operations assignments, etc. Many functions are predefined, or can be
made available by importing R *packages* (more on that later). A function
usually gets one or more inputs called *arguments*. Functions often (but not
always) return a *value*. A typical example would be the function `sqrt()`. The
input (the argument) must be a number, and the return value (in fact, the
output) is the square root of that number. Executing a function ('running it')
is called *calling* the function. An example of a function call is:

```{r, eval=FALSE, purl=FALSE}
b <- sqrt(a)
```

Here, the value of `a` is given to the `sqrt()` function, the `sqrt()` function
calculates the square root, and returns the value which is then assigned to
variable `b`. This function is very simple, because it takes just one argument.

The return 'value' of a function need not be numerical (like that of `sqrt()`),
and it also does not need to be a single item: it can be a set of things, or
even a dataset. We'll see that when we read data files into R.

Arguments can be anything, not only numbers or filenames, but also other
objects. Exactly what each argument means differs per function, and must be
looked up in the documentation (see below). Some functions take arguments which
may either be specified by the user, or, if left out, take on a *default* value:
these are called *options*. Options are typically used to alter the way the
function operates, such as whether it ignores 'bad values', or what symbol to
use in a plot.  However, if you want something specific, you can specify a value
of your choice which will be used instead of the default.

Let's try a function that can take multiple arguments: `round()`.

```{r, results='show', purl=FALSE}
round(3.14159)
```

Here, we've called `round()` with just one argument, `3.14159`, and it has
returned the value `3`.  That's because the default is to round to the nearest
whole number. If we want more digits we can see how to do that by getting
information about the `round` function.  We can use `args(round)` or look at the
help for this function using `?round`.

```{r, results='show', purl=FALSE}
args(round)
```

```{r, eval=FALSE, purl=FALSE}
?round
```

We see that if we want a different number of digits, we can
type `digits=2` or however many we want.

```{r, results='show', purl=FALSE}
round(3.14159, digits=2)
```

If you provide the arguments in the exact same order as they are defined you
don't have to name them:

```{r, results='show', purl=FALSE}
round(3.14159, 2)
```

And if you do name the arguments, you can switch their order:

```{r, results='show', purl=FALSE}
round(digits=2, x=3.14159)
```

It's good practice to put the non-optional arguments (like the number you're
rounding) first in your function call, and to specify the names of all optional
arguments.  If you don't, someone reading your code might have to look up the
definition of a function with unfamiliar arguments to understand what you're
doing.

### Objects vs. variables

What are known as `objects` in `R` are known as `variables` in many other
programming languages.  Depending on the context, `object` and `variable` can
have drastically different meanings.  However, in this lesson, the two words are
used synonymously.  For more information see:
https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Objects

## Vectors and data types

```{r, echo=FALSE, purl=TRUE}
### Vectors and data types
```

A vector is the most common and basic data type in R, and is pretty much
the workhorse of R. A vector is composed by a series of values, which can be   either numbers or characters. We can assign a series of values to a vector using the `c()` function. For example we can create a vector of animal weights and assign it to a new object `weight_g`:

```{r, purl=FALSE}
weight_g <- c(50, 60, 65, 82)
weight_g
```

A vector can also contain characters:

```{r, purl=FALSE}
animals <- c("mouse", "rat", "dog")
animals
```

The quotes around "mouse", "rat", etc. are essential here. Without the quotes R
will assume there are objects called `mouse`, `rat` and `dog`. As these objects don't exist in R's memory, there will be an error message.

There are many functions that allow you to inspect the content of a
vector. `length()` tells you how many elements are in a particular vector:

```{r, purl=FALSE}
length(weight_g)
length(animals)
```

An important feature of a vector, is that all of the elements are the same type of data.
The function `class()` indicates the class (the type of element) of an object:

```{r, purl=FALSE}
class(weight_g)
class(animals)
```

The function `str()` provides an overview of the structure of an object and its
elements. It is a useful function when working with large and complex
objects:

```{r, purl=FALSE}
str(weight_g)
str(animals)
```

You can use the `c()` function to add other elements to your vector:
```{r, purl=FALSE}
weight_g <- c(weight_g, 90) # add to the end of the vector
weight_g <- c(30, weight_g) # add to the beginning of the vector
weight_g
```

In the first line, we take the original vector `weight_g`,
add the value `90` to the end of it, and save the result back into
`weight_g`. Then we add the value `30` to the beginning, again saving the result
back into `weight_g`.

We can do this over and over again to grow a vector, or assemble a dataset.
As we program, this may be useful to add results that we are collecting or
calculating.

We just saw 2 of the 6 main **atomic vector** types (or **data types**) that R
uses: `"character"` and `"numeric"`. These are the basic building blocks that
all R objects are built from. The other 4 are:

* `"logical"` for `TRUE` and `FALSE` (the boolean data type)
* `"integer"` for integer numbers (e.g., `2L`, the `L` indicates to R that it's an integer)
* `"complex"` to represent complex numbers with real and imaginary parts (e.g.,
  `1+4i`) and that's all we're going to say about them
* `"raw"` that we won't discuss further

Vectors are one of the many **data structures** that R uses. Other important
ones are lists (`list`), matrices (`matrix`), data frames (`data.frame`),
factors (`factor`) and arrays (`array`).


> ### Challenge
>
>
> * We’ve seen that atomic vectors can be of type character,
>   numeric, integer, and logical. But what happens if we try to mix these types in
>   a single vector?
> <!-- * _Answer_: R implicitly converts them to all be the same type -->
>
> * What will happen in each of these examples? (hint: use `class()`
>   to check the data type of your objects):
>
>     ```r
>     num_char <- c(1, 2, 3, 'a')
>     num_logical <- c(1, 2, 3, TRUE)
>     char_logical <- c('a', 'b', 'c', TRUE)
>     tricky <- c(1, 2, 3, '4')
>     ```
>
> * Why do you think it happens?
> <!-- * _Answer_: Vectors can be of only one data type. R tries to convert (coerce)
>   the content of this vector to find a "common denominator". -->
>
> * You've probably noticed that objects of different types get
>   converted into a single, shared type within a vector. In R, we
>   call converting objects from one class into another class
>   _coercion_. These conversions happen according to a hierarchy,
>   whereby some types get preferentially coerced into other
>   types. Can you draw a diagram that represents the hierarchy of how
>   these data types are coerced?
> <!-- * _Answer_: `logical -> numeric -> character <-- logical` -->

```{r, echo=FALSE, eval=FALSE, purl=TRUE}
## We’ve seen that atomic vectors can be of type character, numeric, integer, and
## logical. But what happens if we try to mix these types in a single
## vector?

## What will happen in each of these examples? (hint: use `class()` to
## check the data type of your object)
num_char <- c(1, 2, 3, "a")

num_logical <- c(1, 2, 3, TRUE)

char_logical <- c("a", "b", "c", TRUE)

tricky <- c(1, 2, 3, "4")

## Why do you think it happens?

## You've probably noticed that objects of different types get
## converted into a single, shared type within a vector. In R, we call
## converting objects from one class into another class
## _coercion_. These conversions happen according to a hierarchy,
## whereby some types get preferentially coerced into other types. Can
## you draw a diagram that represents the hierarchy of how these data
## types are coerced?
```