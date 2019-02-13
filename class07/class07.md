Class 7 R Functions and Packages
================
Luke Wang
1/30/2019

Revisiting Functions
====================

Revisiting the first example function we created last week - rescale()

Load in the rescale function with **source()**

``` r
source("http://tinyurl.com/rescale-R")
```

Test the rescale() function with some examples

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
# rescale(c(1:5, "string"))
```

Handling Errors and Unexpected Inputs
-------------------------------------

You can use **warning()** or **stop()** to warn users when using a function warning() will continue to execute the function stop() will stop the action of the function

the *call.* argument can display the line of code that produced the error Here is the example of the error handling of rescale2() function when encountering a non numeric value

``` r
# rescale2(c(1:5,"string"))
```

Using the **is.numeric** function

``` r
is.numeric(1:5)
```

    ## [1] TRUE

Practive writing a new function
-------------------------------

``` r
# Lets define an example x and y
x <- c( 1, 2, NA, 3, NA) 
y <- c(NA,3,NA,3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
# Comparing two vecotrs, or and
is.na(x) | is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE  TRUE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
# This will show you the position of the missing value in a vector
which(is.na(x) & is.na(y))
```

    ## [1] 3

Now we take our wokring snippet and turn it in to a function

``` r
both_na <- function (x,y) {
  # Check for NA elements in both input vectors
  sum(is.na(x) & is.na(y))
}
```

Running the *both\_na* function

``` r
both_na(x,y)
```

    ## [1] 1

Let's try to break the *both\_na* function

``` r
x <-  c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both_na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
# This will produce an error and an answer. The answer will be based on recycling of the value
```

Installing a bioconductor package
=================================

``` r
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("GenomicFeatures")
```
