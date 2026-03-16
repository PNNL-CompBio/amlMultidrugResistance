This directory houses all of the code for analysis of PTRC Exp24/Exp27 proteomics data. There are several versions of data generated through analysis with DIA-NN or Spectronaut, based on considerations such as correcting for FAIMS CV. The first sets of analyses aim to compare the different versions of the data- how similar are they in terms of protein IDs, how do the resulting statistics from an identical analysis of each dataset compare, etc. Additionally, we test how filtering out patient samples and how designing the contrast matrix for limma analysis affect the results within any given dataset. The metrics from these analyses guided our decision for choosing a final dataset to move forward with, which is then analyzed in downstream notebooks.



## previous

This directory contains all of the scripts that were used in prior analysis of Exp24 + Exp27 data. Keeping this all for future reference, but nothing in here is needed for the analysis performed in the directories listed below.



## 01\_initial\_data\_processing

The file "01\_Exp24Exp27\_data\_processing.Rmd" performs the initial steps of ingesting results from DIA-NN software, cleaning the data, mapping from UniProt accession to gene names, evaluating and filtering the data based on missingness, median-centering the data within samples, and writing out sample-level data tables for downstream use.



## 02\_dataset\_comparisons

Comparing results from some basic analysis (CD14 vs CD34) run exactly the same way in the 4 different datasets. Additionally, checking within a single dataset at a time, what differences do we get when a) removing some patient samples that had low protein coverage (as was done previously), and b) comparing all CD14 samples vs all CD34 samples (and controlling for sort type (bead vs. flow)) versus only comparing bead-sorted CD14 vs bead-sorted CD34 populations (as was done previously. For this I create a factor called "cellsort" that combines cell type and sort type (i.e. CD14\_bead, CD14\_flow, CD34\_bead, CD34\_flow) that can be used in the linear model to specify which groups to compare.



## 03\_final\_stats

Having completed initial comparisons across the datasets, opted to move forward with one version of the data for final analyses. This notebook pulls msnsets that were generated and saved in the previous folder (02\_dataset\_comparisons), but this can be updated to source from synapse once things are uploaded there. Here we again look at cell type differences, and begin looking into protein and pathway-level differences from cell type and drug sensitivity comparisons.

