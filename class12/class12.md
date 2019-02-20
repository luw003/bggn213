Class 12 Drug Discovery
================
Luke Wang
2/20/2019

Computer Aided Drug Discovery
=============================

In silico docking of drugs to HIV-1 Protease
--------------------------------------------

We are first going to obtain and inspect the input structrue

``` r
library(bio3d)
pdb.code <- "1hsg" # change the pdb code to anything you like for future use
pdb_hiv <- get.pdb(pdb.code) # download 1hsg pdb file
```

    ## Warning in get.pdb(pdb.code): ./1hsg.pdb exists. Skipping download

``` r
file.name <- read.pdb(pdb_hiv) # using file name to make later use easier
```

Now we have the pdb file we are going to remove the ligand that is currently in the structure and produce two seperate files. One of just the protein and one for the ligand.

We are not hard coding any of hte name

``` r
prot <- trim.pdb(file.name,"protein") # trim the original pdb file to just protein
lig <- trim.pdb(file.name,"ligand") # trim to ligand

prot.filename <- paste0(pdb.code,"_protein.pdb") # storing the name of the protein code with protein.pdb suffix
lig.filename <- paste0(pdb.code,"_ligand.pdb")

write.pdb(prot, file=prot.filename)
write.pdb(lig, file=lig.filename)
```

Docking Ligands to HIV-Protease
-------------------------------

We will now read the docking results

``` r
res <- read.pdb("all.pdbqt", multi = TRUE) # read the pdbqt file generated with auto docking VINA

write.pdb(res, "results.pdb")
```

We will procceed to analyze the results quantitatively by calculating the Root Mean Square Distance

``` r
ori <- read.pdb("ligand.pdbqt")

rmsd(ori,res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

As we can see that RMSD and ranking by the docking is quite different. The RMSD doesnt necesary account for the energy.
