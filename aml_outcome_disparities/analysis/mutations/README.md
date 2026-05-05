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

## 2. Analyze combinations of different mutations, stratified by race
- mutations of interest from [Stiff, et al. Nat. Genet. (2024)](https://doi.org/10.1038/s41588-024-01929-x)
- [2_mutations.ipynb](./2_mutations.ipynb)
    - input: `_cache/meta.arrow`
    - output: `_figures/mutations_upset_race-stratified.png`

## 3. Look at additional metadata comparisons
- [3_metadata_comparisons.ipynb](./3_metadata_comparisons.ipynb)
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
- [4_subset_omics_data.ipynb](./4_subset_omics_data.ipynb)
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
- [5_differential_expression.ipynb](./5_differential_expression.ipynb)
    - perform basic differential expression analysis (Welch's t-test) to identify top features within each block for all mutation state comparisons
    - inputs (for all mutation comparison subsets):
        - `_cache/*_cols.pkl`: tuple of lists (WT columns, mutant columns) for each subset
        - `_cache/*.arrow`: combined -omics data table for each subset
    - outputs (for all mutation comparison subsets x -omes combinations): 
        - significantly up/down regulated feature lists: `_cache/diffex/*.txt`
        - significantly up/down regulated features with sample data: `_cache/diffex/*.csv`
        - volcano plots: `_figures/volcano/*.png`

## 6. Associate lipids with proteins
- [6_associate_lipids_and_proteins.ipynb](./6_associate_lipids_and_proteins.ipynb)
    - Uses coexpression (Spearman correlation) of lipids and proteins (Proteomics and Transcriptomics) to find proteins that are associated with significantly differentially expressed lipids to use as inputs for pathway enrichment
    - inputs: 
        - `_cache/*_cols.pkl`: tuple of lists (WT columns, mutant columns) for each subset
        - `_cache/*.arrow`: combined -omics data table for each subset
    - outputs: 
        - significantly up/down regulated lipid-associated protein lists: `_cache/diffex/*_Lipid-Protein-Association_*.txt`
        - significantly up/down regulated lipid-associated proteins with sample data: `_cache/diffex/*_Lipid-Protein-Association_*.csv`

## 7. Run pathway enrichment analysis (ORA)
- [7a_run_pathway_enrichment.Rmd](./7a_run_pathway_enrichment.Rmd)
    - uses KEGG, Reactome, and MSigDB Hallmark databases to perform pathway enrichment (ORA) to highlight pathways from significantly differentially expressed proteins
    - knit from command-line: `Rscript -e 'rmarkdown::render("7a_run_pathway_enrichment.Rmd")'`
    - inputs: 
        - lists of differentially expressed proteins (from Proteomics, Transcriptomics, or Lipid-Protein associations): `_cache/diffex/*.txt`
    - outputs: 
        - knitted PDF from Rmd: `7a_run_pathway_enrichment.pdf`
        - dot plot figures showing enriched pathways: `_figures/pathway_enrichment/*/*.png`
            - separate subdirectories for the different pathway databases
        - output tables with enriched pathways and metadata: `_cache/diffex/pathway_enrichment/*/*.csv`
            - separate subdirectories for the different pathway databases
- [7b_pathway_enrichment_summary.ipynb](./7b_pathway_enrichment_summary.ipynb)
    - Generate summary plots from pathway enrichment analysis across conditions and -omes
    - inputs:
        - output tables with enriched pathways and metadata: `_cache/diffex/pathway_enrichment/*/*.csv`