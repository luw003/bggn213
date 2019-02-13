Class 6 Writing Functions in R
================
Luke Wang
1/25/2019

File Reading!
-------------

Using **data.table()** and related functions to bring data in to R properly.

Insert a code chunk. You can bring small files into R directly through an URL

``` r
read.table("https://bioboot.github.io/bggn213_S18/class-material/test1.txt", header= TRUE, sep=",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Assign URL to file1

``` r
file1 <- "https://bioboot.github.io/bggn213_S18/class-material/test1.txt"
data1 <- read.csv(file1)
data1
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

Reading table with custom seperator ($)

``` r
file2 <- "https://bioboot.github.io/bggn213_S18/class-material/test2.txt"
data2 <- read.table(file2, header =TRUE, sep="$")
data2
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

If you get a weird data format, try read.table

``` r
file3 <- "https://bioboot.github.io/bggn213_S18/class-material/test3.txt"
data3 <- read.table(file3)
data3
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

Writing Function
================

Writing a function If in the argument, a variable has a =, then it means it has a default value

Basic Addition Function
-----------------------

``` r
add <- function(x, y=1) {
# Sum the input x and y 
  x+y
}
```

Let's use the **add()** function

``` r
add(1)
```

    ## [1] 2

``` r
add(1,5)
```

    ## [1] 6

``` r
add(c(1,2,3,4))
```

    ## [1] 2 3 4 5

``` r
add(c(1,2,3,4),c(1,2,3,4))
```

    ## [1] 2 4 6 8

What happens with the unexpected?

``` r
# add(1,2,3)

# Unused argument
```

``` r
# add(x=1,y="b")
# Non numeric argument to binary operator
```

Writing More Complex Function
-----------------------------

### Rules of Writing a Function

Simplify & Reduce Calculation Duplication

### Rescale Function

You need a *name*, *arguments*, and *body*

``` r
rescale <- function(x) {
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

Testing with Small Vector

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

### Exceptions in Functions

``` r
rescale (c(1,2,NA,3,10))
```

    ## [1] NA NA NA NA NA

Diagnosing Unexpected Result

``` r
z <- c(1,2,NA,3,10)
# range(x)
```

Dealing with NA

``` r
# Range() has a arugment for na.rm
rescale2 <- function(x) {
  rng <- range(x, na.rm=TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2 (c(1,2,NA,3,10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

Some more complex function

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
   if(na.rm) {
     rng <-range(x, na.rm=na.rm)
   } else {
     rng <-range(x)
   }
   print("Hello")
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   return(answer)
   print("is it me you are looking for?")
   if(plot) {
      plot(answer, typ="b", lwd=4)
     print("Please don't ever sing again")
}
   print("I can see it in ...")
   
   # return() will end the function and return whatever variable in the return
}
```

``` r
rescale3(1:10, plot=T)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Anallysis of Protein Drug Interaction with Bio3D
------------------------------------------------

### Using Bio3D

``` r
# Load bio3D
library(bio3d)
s1 <- read.pdb("4AKE")  # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE")  # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y")  # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-20-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-20-3.png)

``` r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class06_files/figure-markdown_github/unnamed-chunk-20-4.png) Examing what **read.pdb** does

``` r
library(bio3d)
pdb <- read.pdb("1hbs")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hbs")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 9104,  XYZs#: 27312  Chains#: 8  (values: A B C D E F G H)
    ## 
    ##      Protein Atoms#: 8760  (residues/Calpha atoms#: 1148)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 344  (residues: 8)
    ##      Non-protein/nucleic resid values: [ HEM (8) ]
    ## 
    ##    Protein sequence:
    ##       VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK
    ##       KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA
    ##       VHASLDKFLASVSTVLTSKYRVHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQ
    ##       RFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGT...<cut>...HKYH
    ## 
    ## + attr: atom, xyz, seqres, helix, calpha,
    ##         remark, call
