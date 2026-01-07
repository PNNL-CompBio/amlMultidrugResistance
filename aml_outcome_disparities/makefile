SHELL := /bin/bash
FIGURE_DIR = analysis
PYTHON_DIR = src/python
R_DIR = src/r
PYTHON_PATH_CMD = "import sys; sys.path.append($(PYTHON_DIR))"
INTEGRATED_PREREQS = $(PYTHON_DIR)/data/ba_global_corrected.csv \
					$(PYTHON_DIR)/data/ba_phospho_corrected.csv \
					$(PYTHON_DIR)/data/pilot_global_corrected.csv \
					$(PYTHON_DIR)/data/pilot_phospho_corrected.csv

.venv: uv.lock
	uv sync

lint: .venv
	black src/python/

uv.lock: pyproject.toml
	uv lock

figure1%: .venv $(INTEGRATED_PREREQS)
	@mkdir -p ./analysis/fig1
	@PYTHONPATH=$(PYTHONPATH):$(PYTHON_DIR) python $(PYTHON_DIR)/pilot/figures/figure1$*.py
	@mv ./figure1$*.svg ./analysis/fig1/

$(INTEGRATED_PREREQS):
	Rscript src/r/proteomics/proteomic_batch_correction.R
