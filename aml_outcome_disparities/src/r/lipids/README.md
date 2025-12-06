# Lipid batch correction

Here we combine two datasets to create a single lipidomic measurement dataset
for the manuscript

## Starting files
We require specific files for the analysis, documented below.

### Metadata
To combine the two cohorts we use the metadata file assembled for this purpose.
The metadata is found [here](https://www.synapse.org/Synapse:syn69692583). 

### Beat AML lipid analysis

Beat AML lipid data can be found [on synapse](https://www.synapse.org/Synapse:syn52121001)
for processing. 

### Black patient lipid analysis

The black patient files can be found
[in the lipids folder](https://www.synapse.org/Synapse:syn71210014) of the
Experiment 26 directory. 

beat_lipids_proc.R: BEAT_AML_aligned_lipids_for_stats 2-2025-08-14.xlsx
data_processing_code/ptrc_lipids_proc.R: PTRC_lipids_POS_for_stats_aligned_with_BEAT.xlsx, 
PTRC_lipids_NEG_aligned_with_BEAT.xlsx, PTRC_run_order.xlsx

## Batch correction steps

For lipids, beat_lipids_proc.R should be run first to process the beat aml lipidomics data.
After processing that data, the file bap_lipid_integration_wmetadatafile.R should be ran.
In this file, data are preprocessed, normalized, integrated, and batch corrected.

## Final combined data

There are four csv files generated for the lipid data:

- normalized_e_data_lipid.csv (normalized log2 abundance values)
- normalized_combat_edata_lipid.csv (normalized log2 abundance values that have been batch corrected via ComBat)
- f_data_lipid.csv (sample information)
- e_meta_lipid.csv (biomolecule information)

These files can be found [in the lipids folder](https://www.synapse.org/Synapse:syn71210014) of the
Experiment 26 directory. 

