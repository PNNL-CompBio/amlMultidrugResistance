# Setup
- Follow instructions [here](../../src/python/README.md) to set up Python environment
    - must be able to import and use functions from [data import module](../../src/python/pilot/data_import.py)
- intermediate data files and generated figures are saved in the (untracked) `_cache` and `_figures` directories, respectively

# Analysis
## 1. Create a Combined Data table to support analyses
- [1_combined_data_table.ipynb](./1_combined_data_table.ipynb)
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

## 2. Analyze sample counts for combinations of different mutations, stratified by race
- [2_mutations.ipynb](./2_mutations.ipynb)
    - input: `meta.parquet` (synapse: `syn########`)
    - output: `_figures/mutations_upset_race-stratified.png`
- mutations of interest (NPM1, FLT3-ITD) from [Stiff, et al. Nat. Genet. (2024)](https://doi.org/10.1038/s41588-024-01929-x)

## 3. 