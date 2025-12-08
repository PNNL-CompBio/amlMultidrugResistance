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

The black patient files can be found 
[on synapse](https://www.synapse.org/Synapse:syn68710369) in the
Experiment 26 directory. 

## Batch correction steps

The batch correction is handled in a single script 
[bap_metab_integration_wmetadatafile.R]()

## Final combined data

There are four csv files generated for the each metabolite mode data for a total of eight csv files:

- hilic_normalized_e_data_lipid.csv (normalized log2 abundance values)
- hilic_normalized_combat_edata_lipid.csv (normalized log2 abundance values that have been batch corrected via ComBat)
- hilic_f_data_lipid.csv (sample information)
- hilic_e_meta_lipid.csv (biomolecule information)
- rp_normalized_e_data_lipid.csv (normalized log2 abundance values)
- rp_normalized_combat_edata_lipid.csv (normalized log2 abundance values that have been batch corrected via ComBat)
- rp_f_data_lipid.csv (sample information)
- rp_e_meta_lipid.csv (biomolecule information)

These files can be found [in the metabolites folder](https://www.synapse.org/Synapse:syn71850016) of the
Experiment 26 directory. 

