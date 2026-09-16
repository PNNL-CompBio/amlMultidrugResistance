# Proteomics batch correction

This directory contains scripts to batch correct phosphoproteomic and global proteomic measurements collected across 
cohorts for effective data integration. These scripts use PNNL's version of ComBat that most effectively handles
missing measurements.

## Starting files

All required files are downloaded via Synapse. These scripts expect that a .txt file, `auth_token.txt`, is saved to the
root directory of this project at `amlMultidrugResistance/aml_outcome_disparities`.

### Phosphoproteomics

Data from the BeatAML 210 cohort is collected from [here](https://www.synapse.org/Synapse:syn25714936); data from 
the racial outcomes pilot study is found [here](https://www.synapse.org/Synapse:syn69075545). 

### Global proteomics

Data from the BeatAML 210 cohort is collected from [here](https://www.synapse.org/Synapse:syn25714254); data from 
the racial outcomes pilot study is found [here](https://www.synapse.org/Synapse:syn69075554). 

## Batch correction steps

Batch correction is performed via the `correct_proteomics` function in`proteomic_batch_correction.R` that leverages 
PNNL's ComBat with missing data compatibility. Running this script directly will natively process both the 
phosphoproteomic and global proteomic measurements.

## Final combined data

When executed, `proteomic_batch_correction.R` produces four .csv files: `ba_global_corrected.csv`, 
`pilot_global_corrected.csv`, `ba_phospho_corrected.csv`, and `pilot_phospho_corrected.csv`. Each contains 
batch-corrected and pre-processed phosphoproteomic and global measurements for the BeatAMl 210 (ba) and racial outcomes
pilot (pilot). These files are stored in `aml_outcome_disparities/src/python/data/` for further python-based analyses.
