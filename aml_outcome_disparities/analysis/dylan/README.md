# Info
## Setup
- Follow instructions [here](../../src/python/README.md) to set up Python environment
    - must be able to import and use functions from [data import module](../../src/python/pilot/data_import.py)
- intermediate data files and generated figures are saved in the (untracked) `_cache` and `figures` directories, respectively
    - create these if they do not exist

# Analysis
## 1. Create a Combined Data table to support analyses
- [combined_data_table.ipynb](./combined_data_table.ipynb)
    - metadata:
        - basic sample-level metadata:
            - Age
            - Sex
            - Race
            - Study
        - Mutation states (important genes highlighted in Eisfeld paper)
            - For all genes, the "Not measured" label is converted to "WT", proportions pre-conversion included for reference
            - IDH1 and IDH2 to be combined into single IDH1+IDH2 mutant flag
            - All mutation states columns will be converted to boolean, True indicates mutant
            - filter to only include samples with race labeled as "White" or "Black"
    - input:
        - data tables produced by the functions from [data import module](../../src/python/pilot/data_import.py)
    - output: 
        - `_cache/meta.arrow`: Metadata table
            - columns: metadata fields
            - rows: samples
        - `_cache/combined.arrow`: Combined -omics data table
            - columns: 
                - indices: "Block", "Feature"
                - samples
            - rows: 
                - features

## 2. Analyze combinations of different mutations, stratified by race
- [mutations.ipynb](./mutations.ipynb)
    - input: `_cache/meta.arrow`
    - output: `_figures/mutations_upset_race-stratified.png`

## 3. Look at additional metadata comparisons
- [metadata_comparisons.ipynb](./metadata_comparisons.ipynb)
    - input: `_cache/meta.arrow`
    - output: 
        - Age:
            - `_figures/age_comparison_all.png`
            - `_figures/age_comparison_NPM1.png`
            - `_figures/age_comparison_NRAS.png`
            - `_figures/age_comparison_WT.png`
        - Sex 
            - `_figures/sex_comparison_all.png`
            - `_figures/sex_comparison_black.png`
            - `_figures/sex_comparison_white.png`

## 4. Create subsets of the -omics data for specific mutation state comparisons
- [subset_omics_data.ipynb](./subset_omics_data.ipynb)
    - input: 
        - `_cache/meta.arrow`
        - `_cache/combined.arrow`
    - output:
        - Black patients: NPM1 mutant vs. WT
            - `_cache/B_NPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/B_NPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_B_NPM1-WT.png`: plot of feature counts before and after subsetting
        - White patients: NPM1 mutant vs. WT
            - `_cache/W_NPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/W_NPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_W_NPM1-WT.png`: plot of feature counts before and after subsetting
        - Black patients: NRAS mutant vs. WT
            - `_cache/B_NRAS-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/B_NRAS-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_B_NRAS-WT.png`: plot of feature counts before and after subsetting
        - White patients: NRAS mutant vs. WT
            - `_cache/W_NRAS-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/W_NRAS-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_W_NRAS-WT.png`: plot of feature counts before and after subsetting
    - during subsetting, the omics data has some filtering/normalization applied:
        - filter columns to only include the selected samples for the subset (from WT and mutant sample groups)
        - filter rows separately for each sample group to remove rows with missingness fraction greater than a specified threshold
        - impute remaining missing values across all subset samples with row median
        - apply row-wise z-score normalization across subset samples 

## 5. DIABLO

### Prepare inputs
- [prep_diablo_data.ipynb](./prep_diablo_data.ipynb)
- [create_diablo_designs.ipynb](./create_diablo_designs.ipynb)

### Run parameter tuning
- run [MomDiablo.R](./MomDiablo.R) with "tune" option
- shell script for running parameter tuning for all subsets: [run_diablo_tune.zsh](./run_diablo_tune.zsh)


