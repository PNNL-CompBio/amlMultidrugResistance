# Metabolite batch correction

Here we combine two datasets to create a single metabolomic measurement dataset
for the manuscript.

## Starting files
We require specific files for the analysis, documented below.

### Metadata
To combine the two cohorts we use the metadata file assembled for this purpose.
The metadata is found [here](https://www.synapse.org/Synapse:syn69692583). 

### Beat AML metabolite analysis

Beat AML metabolite data can be found [on synapse](https://www.synapse.org/Synapse:syn53678273)
for processing. 

### Black patient metabolite analysis

The black patient files can be found *** WHERE??*** of the
Experiment 26 directory. 

*** Javi please specific initial files to start with ***

All datafiles are pulled directly from synapse - I did not process any raw datafiles for this analysis.

## Batch correction steps

*** Javi what order do I run the scripts in? ***

Both scripts do the same thing, however one uses a metadatafile that you prepared (i.e. bap_metab_integration_wmetadatafile.R). Running only one of the scripts is necessary.

## Final combined data

*** Javi what file is the combined metabolomics data? ***

The final combined dataset isn't written out anywhere, but it is created in lines 1142-1290 in bap_metab_integration_wmetadatafile.R (or lines 1092-1240 in bap_metab_integration.R)
