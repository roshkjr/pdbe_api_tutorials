# PDBe API Tutorials

This repository contains marimo notebooks that demonstrate how to fetch and analyze
live data from the PDBe REST API.

## Install

Use `uv` to sync the project environment:

```bash
uv sync
```

## Run notebooks

Open the notebook in marimo edit mode:

```bash
uv run marimo edit notebooks/<python_script.py>
```

Run it as a read-only marimo app:

```bash
uv run marimo run notebooks/<python_script.py>
```

## Notebooks

|Notebook|Description|Open in molab|
|---------|-----------|-------------|
|egfr_ligand_batch_analysis.py|Analysis of ligand bound structures of EGFR in the PDB|    [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/roshkjr/pdbe_api_tutorials/blob/main/notebooks/egfr_ligand_batch_analysis.py)|

