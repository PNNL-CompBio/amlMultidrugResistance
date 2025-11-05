# Sorted proteomics
Investigates proteomic differences between CD14+ and CD34+ AML cells since
CD14+ AML tend to be resistant to venetoclax. Includes PTRC2 experiments 
24 (DIA & TMT proteomics) and 27 (DIA proteomics).

## Install dependencies
installDependencies.R

## Initial processing: [DIA] (DIA/) (experiments 24 and 27)
### [Experiment 24] (DIA/exp24/)
DIA/exp24/final_processing_DIA_sorted_cells.Rmd

### [Experiment 27] (DIA/exp27/)
DIA/exp27/final_processing_DIA_sorted_cells.Rmd

## Initial processing: [TMT] (TMT/) (experiment 24)
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
Figures 1C, 4D; Table S1, S3: 5-panSEA_global,phospho_20250703.R

### Differential expression of published scRNAseq data
241202_diffexp_limma.R

### Correlation-based analyses
corrWithBulkAUC.R
Not currently in use.

### Signature comparison
Figures 2, 5B-C: Table S3: 2025-01-20_v2_MonoSigComparison_cellFrac.R

### Signature refinement for venetoclax AUC correlation
Figures 5D, S2: Table S5: 2025-02-20_drugSensPredictionAccuracy.R

### Other figure generation
- Figure 3: fig3.R
- Figure S1: figS1.R
- Figures 1D, 4A-C: figures_2025-01-20.R
