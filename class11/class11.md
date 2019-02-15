Class 10 Structrual Bioinformatics
================

PDB Statistics
--------------

Analyzing the PDB statistics about percentage of structures solved by X-Ray and Electron Microscopy

``` r
stats <- read.csv("Data Export Summary.csv", row.names=1) # read csv file and set first column as row name
total <- sum(stats$Total) # get total number of structures
xray <- round(stats$Total[1]/total*100,2) # calculate the percentage of x ray solved structures
em <- round(stats$Total[3]/total*100,2) # calculate the percentage of em solved structures
paste(xray,"% of total structures are solved by X-Ray")
```

    ## [1] "89.43 % of total structures are solved by X-Ray"

``` r
paste(em,"% of total structures are solved by Electron Microscopy")
```

    ## [1] "1.89 % of total structures are solved by Electron Microscopy"

``` r
# Here is an alternative way to get all the percentages
per.by.method <- round(stats$Total/sum(stats$Total)*100,1)
names(per.by.method) <- rownames(stats) # assign row names to each of the vector
paste(per.by.method,"% of total structures are ",rownames(stats), sep = '')
```

    ## [1] "89.4% of total structures are X-Ray"             
    ## [2] "8.4% of total structures are NMR"                
    ## [3] "1.9% of total structures are Electron Microscopy"
    ## [4] "0.2% of total structures are Other"              
    ## [5] "0.1% of total structures are Multi Method"

Let's look at what proportion of structures are actually protein

``` r
per.by.protein <- round(sum(stats$Proteins)/sum(stats$Total)*100,1)
paste(per.by.protein,"% of total structures are proteins", sep = '')
```

    ## [1] "92.8% of total structures are proteins"

Getting started with Bio3D and Bio3D View
-----------------------------------------

``` r
library(bio3d)
library(bio3d.view)
```

Loading the 1HSG pdb file

``` r
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

Using a developmental tool - bio3d-view to view 1HSG in R

``` r
view(pdb,"overview",col="sse")
```

    ## Computing connectivity from coordinates...

Atom Selection with bio3d
-------------------------

Selecting all alpha carbons

``` r
# Selecting all the alpha carbon atoms and return their indicies

calpha.inds <- atom.select(pdb,"calpha")
calpha.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
    ## 
    ##    Atom Indices#: 198  ($atom)
    ##    XYZ  Indices#: 594  ($xyz)
    ## 
    ## + attr: atom, xyz, call

Selecting just the chain A

``` r
a.inds <- atom.select(pdb,chain="A")
# Selecting C-alphs of chain A
ca.inds <- atom.select(pdb,"calpha",chain="A")
```

Extract the protein only portion of this PDB structure and write it out to a new PDB file

``` r
pro.inds <- atom.select(pdb,"protein")
pro <- trim.pdb(pdb,inds=pro.inds)
write.pdb(pdb=pro,file="1hsg-no-mk1.pdb")
```

Extract the ligand (i.e drug) and write it to a new PDB file

``` r
ligand.inds <- atom.select(pdb,"ligand")
mk1 <- trim.pdb(pdb,inds=ligand.inds)
write.pdb(pdb=mk1, file="mk1.pdb")
```
