# Patient outcome disparities in AML: Python

This src directory contains Python-based analysis scripts supporting the AML outcome disparities project.

## Project Organization

├── src/python/ \
│   ├── data/                       # Utility and intermediate storage \
│   │   └── van_Galen_genes.csv     # Genes for van Galen cell types \
│   └── pilot/                      # Python module \
│       ├── figures/                # Preliminary data exploration figures \
│       │   └── figure_setup.py     # Figure plotting utilities (to ensure consistent formatting) \
│       ├── data_import.py          # Data importation utilities \
│       ├── decomposition.py        # Data-driven decomposition tools \
│       ├── gene_analysis.py        # Transcriptomic analysis tools \
│       └── utils.py                # Miscellaneous supporting utilities \
├── pyproject.toml                  # uv dependency file \
├── README.md                       # README for Python module \
└── uv.lock                         # Locked uv dependency file

## Python Installation & Package Management

This project uses [uv](https://docs.astral.sh/uv/) for Python dependency management. You can install uv following 
the [instructions here](https://docs.astral.sh/uv/#installation).

Once installed, you can build a Python environment by navigating to this directory 
(`aml_outcome_disparities/src/python`) and running:

```
uv lock
uv sync
```

This will download Python and any required dependencies as well as set up a virtualenv (under `.venv/`) that will 
automatically handle connections between Python and any required dependencies.

## Data Prerequisites

Before running any of the Python-based figure scripts, you'll need to add a couple of files to the `data/` 
directory.

1. Batch-corrected phosphoproteomic and global proteomic measurements for the BeatAML 210 and Pilot study cohorts. 
   These are automatically generated and saved to `data/` by executing the proteomic batch correction script
   `src/r/proteomics/proteomic_batch_correction.R`. This batch correction is done via R to leverage PNNL's improved 
   ComBat method that handles missing values.
2. Gene lengths to convert RNA-seq counts to transcripts-per-million, saved as `data/biomart_gene_lengths.txt.gz`. 
   These can be retrieved from [BioMart](https://useast.ensembl.org/info/docs/index.html), and require the columns 
   `Transcript length (including UTRs and CDS)` and `Gene stable ID`.
    * TODO: Move this to an automated .xml-based import for reproducibility.

Data on Synapse is downloaded directly through Synapse's Python client. To use, you'll need to include a `.txt` file 
with a personal access token that you can generate [here](
https://accounts.synapse.org/authenticated/personalaccesstokens). Once generated, save the token to `auth_token.txt` in
the `src/python/` directory.

## Code Execution

Files within the `pilot/` directory contain useful tools for analysis and data importation that are *not* intended 
to be executed directly. Rather, you can execute the figure files in the top-level 
`aml_outcome_disparities/analysis/` directory that will call functions from these scripts instead to generate 
manuscript figures.

To execute figure scripts, you can use the following command to ensure the script is run within the uv virtual 
environment:

```
uv run /path/to/figure.py
```

Each Python-based figure script has been configured to add scripts from the `pilot/` directory to the current 
path, so the above `uv run` command can be executed from any location within this repository. Alternatively, you can 
view the output of each figure script within the notebooks in the `aml_outcome_disparities/analysis/` directory that 
will automatically run and display the outputs of each figure script.
