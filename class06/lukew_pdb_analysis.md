PDB Structure Analysis
================
Luke Wang
1/28/2019

Introduction
------------

Here we will create a function **pdb\_analysis()** to analyze .pdb files using the Bio3D package.

The function in its current state will take any vector of PDB accession numbers, and analyze all the C-alpha atoms on the "A" chain of the protein. The results will be returned in individual plots of each residue plotted against its B-factor as well as secondary structure information on the top and bottom margin of the plot.

Inputs for pdb\_analysis()
--------------------------

There are only two arguments for the function.

First is the object that contains pdb accession numbers. This argument has no default value.

Second is the *plot* argument. This argument will set whether to plot the b factor values. The default is set at TRUE. If the argument is set to false, there will be no output from this function.

Output for pdb\_analysis()
--------------------------

If the arugment *plot* is set to true, then the function will generate plots of B factors of CA atoms on the A chain for each protein that is being analyzed. If argument *plot* is set to false, there is no output from the function.

Original Code
-------------

The objective of this excercise is to reduce the complexity of the original code listed below

``` r
## Original Code
# s1 <- read.pdb("4AKE")  # kinase with drug
# s2 <- read.pdb("1AKE")  # kinase no drug
# s3 <- read.pdb("1E4Y")  # kinase with drug
# s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
# s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
# s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
# s1.b <- s1.chainA$atom$b
# s2.b <- s2.chainA$atom$b
# s3.b <- s3.chainA$atom$b
# plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
# plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
# plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

New Function!
-------------

Here we will create a function based off of the original code listed above that will reduce the complexity of this analysis.

``` r
# Load bio3D packages for functions to analyze PDB files
library(bio3d)

# Establishing the pdb_analysis function with 2 arguments
pdb_analysis <- function(protein, plot=TRUE){
  for (x in protein){
    protein1 <- read.pdb(x)
    pro1.chainA <- trim.pdb(protein1, chain ="A", elety = "CA")
    pro1.b <- pro1.chainA$atom$b
    if (plot == TRUE){
      plotb3(pro1.b, sse=pro1.chainA, typ="l", ylab="Bfactor")
    }
  }
}
```

Some Examples
-------------

``` r
# Example with a vector of 3 proteins, with default plot setting
pdb_analysis(c("4AKE", "1AKE", "1E4Y"))
```

    ##   Note: Accessing on-line PDB file

![](lukew_pdb_analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

![](lukew_pdb_analysis_files/figure-markdown_github/unnamed-chunk-3-2.png)

    ##   Note: Accessing on-line PDB file

![](lukew_pdb_analysis_files/figure-markdown_github/unnamed-chunk-3-3.png)

``` r
# Example with a string of single protein, with default plot setting
pdb_analysis("4AKE")
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 5v/mdyrfcnd70qdx50pkx1w5hzm0000gn/T//Rtmp4WzuGi/4AKE.pdb exists. Skipping
    ## download

![](lukew_pdb_analysis_files/figure-markdown_github/unnamed-chunk-3-4.png)

``` r
# Example with a vector of 3 proteins but no plots
pdb_analysis(c("4AKE", "1AKE", "1E4Y"), plot = FALSE)
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 5v/mdyrfcnd70qdx50pkx1w5hzm0000gn/T//Rtmp4WzuGi/4AKE.pdb exists. Skipping
    ## download

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 5v/mdyrfcnd70qdx50pkx1w5hzm0000gn/T//Rtmp4WzuGi/1AKE.pdb exists. Skipping
    ## download

    ##    PDB has ALT records, taking A only, rm.alt=TRUE
    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 5v/mdyrfcnd70qdx50pkx1w5hzm0000gn/T//Rtmp4WzuGi/1E4Y.pdb exists. Skipping
    ## download

``` r
# No useful output is produced from the previous example. The protein that was analyzed will not be passed outside of the function for furthur use.
```
