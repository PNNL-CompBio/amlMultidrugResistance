# NRAS
Investigates changes (proteomic or transcriptomic) upon NRAS knockdown since
NRAS mutations arise with acquired resistance to gilteritinib. Includes PTRC2
experiments 20, 21, and 22.

## Initial processing: for each experiment in its respective subfolder (exp20, exp21, or exp22)
### Create study design tables
0-create_study_design_tables.Rmd

### Initial TMT global proteomics processing
1-process_global_data.Rmd

### Initial TMT phospho proteomics processing
2-process_phospho_data.Rmd

### Prep KSTAR input
2.5-process_KSTAR_input.Rmd
Not currently used but available in case of future need.

### Normalization and batch correction
3-normalize_and_batch_correction.Rmd

### Upload crosstabs to Synapse
4-push_to_synapse.Rmd

## Analysis
### Differential expression and GSEA
4-panSEA_global,phospho_20250xxx.R
Uses functions available in the helperFunctions subfolder.

### Deconvolution (figure 4)
2025-01-20_v2_MonoSigComparison_cellFrac.R

### Other figure generation
Figure 3: exp20/Figure3_20250528.R
Figure 6: exp21/Figure6_20250430.R
Other figures for experiment 22: exp22/suppFigExp22_20250430.R
