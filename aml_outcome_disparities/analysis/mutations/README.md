# Setup
- Follow instructions [here](../../src/python/README.md) to set up Python environment
    - must be able to import and use functions from [data import module](../../src/python/pilot/data_import.py)
- intermediate data files and generated figures are saved in the (untracked) local `_cache` and `_figures` directories, respectively

# Analysis

## 1. Create a Combined Data table to support analyses

### [1_combined_data_table.ipynb](./1_combined_data_table.ipynb)
Assembles metadata and combined -omics data tables from individual data tables. Each individual -omics data table is subjected to null-filtering/filling and normalization as as needed, then all data tables are collected into a single -omics data table to support downstream analyses. 

- metadata:
    - basic sample-level metadata:
        - Age
        - Sex
        - Race
        - Study
    - Mutation states
        - important genes highlighted in Eisfeld paper ([Stiff, et al. Nat. Genet. (2024)](https://doi.org/10.1038/s41588-024-01929-x)): NPM1, FLT3-ITD
            - these are also the most abundant across full cohort
        - For all genes, the "Not measured" label is converted to "WT", proportions pre-conversion included for reference
        - All mutation states columns will be converted to boolean, True indicates mutant
        - filter to only include samples with race labeled as "White" or "Black"
- inputs:
    - data tables produced by the functions from [data import module](../../src/python/pilot/data_import.py)
        - `import_acetyl`
        - `import_global`
        - `import_lipids`
        - `import_meta`
        - `import_metabolites`
        - `import_phospho`
        - `import_rna`
- outputs: 
    - generates cached local files and uploads/updates them to synapse
    - Metadata table
        - columns: metadata fields
        - rows: samples
        - files: 
            - `_cache/metadata.csv` → (synapse: `syn74775940`)
            - `_cache/metadata.parquet` → (synapse: `syn74775994`)
    - Combined -omics data table
        - columns: 
            - indices: "OmicsDataType", "Feature"
            - samples
        - rows: 
            - features
        - files: 
            - `_cache/combined_omics_data.csv` → (synapse: `syn74775999`)
            - `_cache/combined_omics_data.parquet` → (synapse: `syn74775996`)


## 2. Analyze sample counts for combinations of different mutations, stratified by race

### [2_mutations.ipynb](./2_mutations.ipynb)
Generates a figure summarizing the counts of samples by race and mutation status.

- input: 
    - metadata table parquet (synapse: `syn74775940`)
- output: 
    - `_figures/mutations_upset_race-stratified.png`


## 3. Summarize -omics data coverage, stratified by race and mutation

### [3_omics_data_coverage.ipynb](./2_mutations.ipynb)
Generates a figure summarizing the coverage of different -omics data types by race and mutation status.

- inputs: 
    - metadata table, parquet (synapse: `syn74775940`)
    - combined omics data table, parquet (synapse: `syn74775996`)
- output: 
    - `_figures/mutations_omics_coverage_sample_counts.png`


## 4. Use `limma` to perform differential expression across mutation/race comparisons using different -omics data types

### [4a_set_up_limma.ipynb](./4a_set_up_limma.ipynb)
Defines all pairwise comparisons and creates the input parquets for differential expression analysis using `limma`.

- input: 
    - metadata table parquet (synapse: `syn74775940`)
- outputs:
    - `_cache/limma/design.parquet`
    - `_cache/limma/contrast.parquet`

### [4b_run_limma.R](./4b_run_limma.R)
Runs `limma` using combined -omics data table and pre-defined contrasts/design tables

- inputs:
    - _all inputs are read from local cached files_
    - `_cache/combined_omics_data.parquet`
    - `_cache/limma/design.parquet`
    - `_cache/limma/contrast.parquet`
- output:
    - `_cache/limma/results.parquet`

### [4c_analyze-limma-results.ipynb](./4c_analyze-limma-results.ipynb)
Generates summary figures for differentially expressed features from `limma` results, optionally upload/update results table to synapse.

- inputs: 
    - `_cache/limma/results.parquet` 
        - optionally, upload/update → (synapse: `syn78408378`)
- outputs:
    - `_figures/limma_results_sig_up_down.png`
    - `_figures/limma_results_sig_up_down_padj.png`