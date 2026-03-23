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

Before running any of the Python-based figure scripts, you'll need to generate an authorization token for Synapse. 
All biological and clinical meta data are stored via Synapse for reproducibility and transparency, and the scripts 
within this repository will automatically download files from Synapse as needed for analyses.

You can generate a personal access token [here](https://accounts.synapse.org/authenticated/personalaccesstokens). 
Once generated, save the token to `src/python/auth_token.txt`. This will allow the data import scripts to access any 
files you have permission to access on Synapse.

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
