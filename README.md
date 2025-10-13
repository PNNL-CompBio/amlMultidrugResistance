# AML Multidrug resistance studies
This repository contains the analysis and processing for multiple projects focusing on multidrug resistance for the CPTAC funded OHSU/PNNL PTRC. 

## Basic data processing
We focus on three drugs: gilteritinib (gilt), decitabine (dec), and
venetoclax (ven).  We used single resistance models of gilt and dec,
plus dual resistance models of gilt/ven, gilt/dec, and finally
gilt/ven/dec. Enrichment analysis was performed with enrichR through
MAGINE. Current data processing and notebooks to do general analysis
are in the [data](./data) and [notebooks](./notebooks) directories respectively.

### Data location
All data are located on Synapse at
[http://synapse.org/ptrc2](http://synapse.org/ptrc2) and can only be
accessed with permission at the moment. 


## Manuscript-specific analysis code
Our current plan is to contribute to three individual manuscripts,
each one describing a new set of combination cell lines. 

### Gilteritinib/Venetoclax combination resistance
Here we explore the phenotypes that emerge during early and late
resistance to Gilteritinib and Venetoclax in combination. Current
analyses are located in [venGiltResistance](./venGiltResistance).


### Gilteritinib/Decitabine combination resistances
Here we explore the phenotype that emerges upon treatment with
Decitabine and Gilteritinib. Analyses are located in
[decGiltResistance](./decGiltResistance). 

### Triple resistance
The last manuscript explores resistance to Gilteritinib, Decitabine,
and Venetoclax in combination. Analyses are located in
[tripleResistance](./tripleResistnace). 

## Other manuscripts
Other manuscripts related to resistance to gilteritinib or venetoclax:

### NRAS
Investigates changes (proteomic or transcriptomic) upon NRAS knockdown since
NRAS mutations arise with acquired resistance to gilteritinib. Proteomic 
analyses are located in [NRAS](./NRAS).

### Sorted proteomics
Investigates proteomic differences between CD14+ and CD34+ AML cells since
CD14+ AML tend to be resistant to venetoclax. Analyses are located in 
[sortedProteomics](./sortedProteomics).