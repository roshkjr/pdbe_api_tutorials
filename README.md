# PDBe API Tutorials

This repository contains marimo notebooks that demonstrate how to fetch and analyze
live data from the PDBe REST API.

## Install

Use `uv` to sync the project environment:

```bash
uv sync
```

## Run the EGFR notebook

Open the notebook in marimo edit mode:

```bash
uv run marimo edit notebooks/egfr_ligand_batch_analysis.py
```

Run it as a read-only marimo app:

```bash
uv run marimo run notebooks/egfr_ligand_batch_analysis.py
```

## Notes

- The notebook uses the live PDBe API and requires network access.
- The workflow starts from `TARGET_UNIPROT_ID = "P00533"` and focuses on EGFR-bound ligands.
- The analysis highlights a `GET` request to fetch EGFR ligand data and a batch `POST`
  request to fetch compound summaries.
