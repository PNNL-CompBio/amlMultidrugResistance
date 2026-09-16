Folder for DIA processing scripts and notes


Steps to DIA-NN library free processing:
1 - Raw files were converted to mzml format in MSConvert using the peak picking algorithm. 
2 - mzml files and relevant fasta (named "Human2023.fasta") were given as input to DIA-NN, which was run in library free mode (see "report.log.txt" for specific parameters)
3 - protein level output (pg_matrix.tsv) was used for subsequent processing in "final_processing_DIA_sorted_cells.Rmd", which contains SamO and Camilo's processing steps
