# Info
## Setup
- Follow instructions [here](../../src/python/README.md) to set up Python environment
    - must be able to import and use functions from [data import module](../../src/python/pilot/data_import.py)
- intermediate data files and generated figures are saved in the (untracked) `_cache` and `_figures` directories, respectively
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
- mutations of interest from [Stiff, et al. Nat. Genet. (2024)](https://doi.org/10.1038/s41588-024-01929-x)
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
- [Singh, et al. Bioinformatics. (2019)](https://doi.org/10.1093/bioinformatics/bty1054)
- In DIABLO parlance a "block" is a single -omics dataset, a component of a multiomics dataset.

### Prepare inputs
- [prep_diablo_data.ipynb](./prep_diablo_data.ipynb)
    - prepares input CSVs (one per block) for DIABLO analysis
- [create_diablo_designs.ipynb](./create_diablo_designs.ipynb)
    - writes design matrices to CSV files (input for DIABLO analysis):
        - null design (all 0s)
            - model not constrained by pre-defined correlations between blocks 
        - design based on correlation between component 0 from single-block PLS-DA
            - correlating compounent 0 projections from PLS-DA on single blocks is meant to capture the extent to which the different blocks capture variance related to phenotype (Mutant vs. WT)
            - this imparts strong constraints on model optimization, usually leading to poorer classification performance but alsp hopefully a more coherent signature spanning the blocks.

### Run parameter tuning
- run [MomDiablo.R](./MomDiablo.R) with "tune" option
- shell script for running parameter tuning for all subsets: [run_diablo_tune.zsh](./run_diablo_tune.zsh)
- [create_diablo_keepX.ipynb](./create_diablo_keepX.ipynb)
    - Based on parameter tuning results, creates CSVs with the number of features to keep per block for each component (input for DIABLO analysis)

### Run final DIABLO analysis
- run [MomDiablo.R](./MomDiablo.R) with "final" option
- shell script for running final analysis for all subsets: [run_diablo_final.zsh](./run_diablo_final.zsh)
    - distance metric and number of components parameters are selected based on results from the tuning

