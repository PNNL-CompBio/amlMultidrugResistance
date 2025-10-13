# Sorted proteomics
Investigates proteomic differences between CD14+ and CD34+ AML cells since
CD14+ AML tend to be resistant to venetoclax. Includes PTRC2 experiments 
24 (DIA & TMT proteomics) and 27 (DIA proteomics).

## Initial processing: DIA (experiments 24 and 27)
### Experiment 24 (exp24 folder)
final_processing_DIA_sorted_cells.Rmd

### Experiment 27 (exp27 folder)
final_processing_DIA_sorted_cells.Rmd

## Initial processing: TMT (experiment 24)
### Create study design tables
TMT/0-create_study_design_tables.Rmd

### Initial TMT global proteomics processing
TMT/1-process_global_data.Rmd

### Initial TMT phospho proteomics processing
TMT/2-process_phospho_data.Rmd

### Prep KSTAR input
TMT/2.5-process_KSTAR_input.Rmd
Not currently used but available in case of future need.

### Normalization and batch correction
- TMT/3-normalize_and_batch_correction.Rmd
- Also tried adding DIA phospho data in 
TMT/3-normalize_and_batch_correction-wDIAphospho_20250102.Rmd
but not currently in use.

### Upload crosstabs to Synapse
TMT/4-push_to_synapse.Rmd

## Analysis
### Differential expression and GSEA
4-panSEA_global,phospho_20250703.R

### Differential expression of published scRNAseq data
241202_diffexp_limma.R

### Correlation-based analyses
corrWithBulkAUC.R
Not currently in use.

### Signature comparison
2025-01-20_v2_MonoSigComparison_cellFrac.R

### Signature refinement for venetoclax AUC correlation
2025-02-20_drugSensPredictionAccuracy.R

### Other figure generation
- Figure 3: fig3.R
- Figure S1: figS1.R
- Other figures: figures_2025-01-20.R
