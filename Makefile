# Reproduce everything from scratch:  make all
# Set VCG_DATA_DIR to keep the ~9 GB of raw data outside the repository.
PYTHON ?= python
export PYTHONPATH := src

.PHONY: all install data prepare evaluate summarise methods compare notebook test clean

all: data prepare evaluate summarise methods compare

install:
	$(PYTHON) -m pip install -r requirements.txt

data:
	$(PYTHON) data/download_ntp.py $(DOWNLOAD_FLAGS)

prepare:
	$(PYTHON) -m vcg.cli prepare

evaluate:
	$(PYTHON) -m vcg.cli evaluate

summarise:
	$(PYTHON) -m vcg.cli summarise

methods:
	$(PYTHON) -m vcg.cli methods

compare:
	$(PYTHON) -m vcg.cli compare

notebook:
	$(PYTHON) -m jupyter nbconvert --to notebook --execute --inplace notebooks/walkthrough.ipynb

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/tables/*.csv results/tables/*.parquet results/figures/*.png results/*.json
