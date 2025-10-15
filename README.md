# AML Multidrug resistance studies
This repository contains the analysis and processing for multiple
projects focusing on multidrug resistance for the CPTAC funded
OHSU/PNNL PTRC. 



All data for this project resides on the Synapse at
[http://synapse.org/ptrc2]. To run the code for each project, you will
need to navigate to the site nad request access. 


## Directory structure

The repository is currently laid out as follows, with folders divided roughly according to
dataset or experiment that was analyzed

### Cell line work

In this project we developed a number of drug resistant cell lines (MOLM14) that were
cultured in combination with Gilteritinib, Venetoclax, Decitabine, or some combination 
there of. These data can be found here.

- [./prelim_analysis](./prelim_analysis): Here we have exploratory work
  from early projects and data, mainly for historic purposes.
- [./combination_cell_line_proteomics](./combination_cell_line_proteomics):
  Here we have data analysis from MOLM14 cells treated with
  combiantions of venetoclax, gilteritinib and HMA
- [./venetoclax_gilteritinib_resistance](./venetoclax_gilteritinib_resistance):
  This is analysis of MOLM14 cell lines that are resistant to ven +
  gilt in combination. Currently going into a manuscript.
- [./NRAS_ASO](./NRAS_ASO): This is analysis for a study focusing on
  the impacs of NRAS ASO on gilteritinib resistance. 

### Ex vivo cultures
We have also measured proteins expressed in ex vivo culture from patient samples. 
This analysis can be found below.

- [./drug_treated_samples](./drug_treated_samples): analysis of patient samples cultured for
months of treatment with single agents or combinations of drugs. 

### Patient samples
We also measured proteomics measurements of patient samples, all analysis for these is in the 
following directories.
- [./cd14_cd34_proteomics](./cd14_cd34_proteomics): This is analysis
  comparing expression of AML patient samples sorted by CD14+ and
  CD34+ markers. 
- [./aml_outcome_disparities](./aml_outcome_disparities): This
  analysis focuses on the identification of multiomic mechanisms of
  outcome disparities in Black AML patients. 
- [./ven_hma_samples](./ven_hma_samples): Here we have new analysis
  of data from patient samples before treatemnt with ven+HMA. 
