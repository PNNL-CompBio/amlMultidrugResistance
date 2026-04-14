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
            - `_figures/age_comparison_FLT3-ITD.png`
            - `_figures/age_comparison_NPM1xFLT3-ITD.png`
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
        - summary of omics data coverage for samples: `_figures/mutations_omics_coverage_sample_counts.png`
        - Black patients: NPM1 mutant vs. WT
            - `_cache/B_NPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/B_NPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_B_NPM1-WT.png`: plot of feature counts before and after subsetting
        - White patients: NPM1 mutant vs. WT
            - `_cache/W_NPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/W_NPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_W_NPM1-WT.png`: plot of feature counts before and after subsetting
        - Black patients: FLT3-ITD mutant vs. WT
            - `_cache/B_FLT3-ITD-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/B_FLT3-ITD-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_B_FLT3-ITD-WT.png`: plot of feature counts before and after subsetting
        - White patients: FLT3-ITD mutant vs. WT
            - `_cache/W_FLT3-ITD-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/W_FLT3-ITD-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_W_FLT3-ITD-WT.png`: plot of feature counts before and after subsetting
        - Black patients: FLT3-ITDxNPM1 mutant vs. WT
            - `_cache/B_FLT3-ITDxNPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/B_FLT3-ITDxNPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_B_FLT3-ITDxNPM1-WT.png`: plot of feature counts before and after subsetting
        - White patients: FLT3-ITDxNPM1 mutant vs. WT
            - `_cache/W_FLT3-ITDxNPM1-WT_cols.pkl`: tuple of lists (WT columns, mutant columns)
            - `_cache/W_FLT3-ITDxNPM1-WT.arrow`: combined -omics data table for selected subset
            - `_figures/block-feature-counts_W_FLT3-ITDxNPM1-WT.png`: plot of feature counts before and after subsetting
    - during subsetting, the omics data has some filtering/normalization applied:
        - filter columns to only include the selected samples for the subset (with two sample groups)
        - filter rows to remove those that do not have at least 5 non-null values in WT and mutant groups
        - apply row-wise z-score normalization across subset samples that are not null
            - require at least 5 samples within WT and mutant groups to z-score normalize

## 5. Single -ome differential expression
- [differential_expression.ipynb](./differential_expression.ipynb)
    - perform basic differential expression analysis (Welch's t-test) to identify top features within each block for all mutation state comparisons

## 6. Associate lipids with proteins
- [associate_lipids_and_proteins.ipynb](./associate_lipids_and_proteins.ipynb)
    - Uses coexpression (Spearman correlation) of lipids and proteins (Proteomics and Transcriptomics) to find proteins that are associated with significantly differentially expressed lipids
    - This generates lists of proteins that can be used for pathway enrichment
