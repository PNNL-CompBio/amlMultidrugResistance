# Patient outcome disparities in AML

This work focuses on understanding why black patients with AML fare worse than white patients with 
the same molecular profiles. Much of this work is premised on a publication and samples 
from [Stiff et al.](https://www.nature.com/articles/s41588-024-01929-x). 

## Experiment 26 Data
From the ~100 or so patient samples we measured:
- proteomics
- phosphoproteomics
- metabolomics
- lipidomics
- acetylomics
That data is all stored in the [Experiment 26 directory](https://www.synapse.org/Synapse:syn62750464) 
on Synapse and currently under embargo. Additionally we plan to integrate data from the Stiff et al paper above.

## Project goals
The goal of this project is to understand why known genetic biomarkers of AML prognosis 
fail to predict response in black patients, and possibly see if we can find better
predictors.

## Data processing and batch correction
The first step of this analysis is to ensure that we can merge two cohorts together: the Beat AML patient cohort together with the
primarily black cohort from Ohio. 

### Lipids
Lipid batch correction information can be found in the 
[lipid directory](src/r/lipids) with the processed data uploaded to synapse
(where is final file?)

### Metabolites
- Metabolites: [src/r/metab](src/r/metab)

### Proteins and phosphosites

### Acetyl sites

### Figures
The results of the batch correction should be documented in the initial
manuscript figures under [analysis/Fig1](analysis/Fig1). 

