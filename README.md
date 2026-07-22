# Hi-C Sequencing Depth Requirements

### Code accompanying the paper

**"Systematic assessment of sequencing depth requirements for Hi-C-derived metrics"**

*A. Granovsky, K. Polovnikov*

---

## Overview

Systematic subsampling of twelve published Hi-C libraries ([Akgol Oksuz et al., 2021](https://doi.org/10.1038/s41592-021-01248-7); [Rao et al., 2014](https://doi.org/10.1016/j.cell.2017.09.026)) from ~360 M down to 1 M read pairs, with four Hi-C-derived metrics recomputed at every depth:

- **P(s)-derived loop density and loop size** — inferred via the log-derivative "dip" and "peak", respectively.
- **Compartment eigenvectors (E1).**
- **Insulation track and boundary calls.**

For each metric we compare the scatter between independent subsamples at a given depth against a biological floor — the variability between independent biological replicates at full depth. The minimum useful depth is the point at which subsampling scatter falls to that floor: below it the estimate is limited by sequencing, above it by biology.

Loop density and loop size are inferred using the model of [Polovnikov & Starkov (2026)](https://doi.org/10.1073/pnas.2534385123).

---

## Data

Raw `.fastq` files: [GSE163666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163666) (Oksuz et al.) and [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) (Rao et al.).

`.mcool` files were produced with the [distiller-nf](https://github.com/open2c/distiller-nf) pipeline (config in `distiller-nf/`).

Intermediate files are included in `data/`, so the analysis notebooks can be re-run without reprocessing raw `.mcool` files.

---

## Environment

A [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container definition is provided in `hic-tools.def`. It bundles all dependencies needed for the analysis (Python, cooltools, the Open2C ecosystem, etc.). To build and use:

```bash
singularity build hic-tools.sif hic-tools.def
singularity exec hic-tools.sif jupyter lab
```

---

## Repository Structure

```
├── data/
│   ├── compartments/            # E1 trackss per depth
│   ├── datasets/                # Subsampled .mcool metadata and sample tables
│   ├── insulation/              # Insulation profiles and boundary calls per depth
│   └── slopes/                  # P(s) curves and log-derivatives per depth
│
├── distiller-nf/
│   ├── cluster.config           # config for distiller-nf
│   └── params_file.yml     	 # pipeline parameters
│
├── figures/					 # Figure output folder
│
├── notebooks/
│   ├── compute_slopes.ipynb             # P(s) curves and log-derivatives from raw .mcool files
│   ├── loop_size_and_density.ipynb      # Loop size / density inference and stability vs. depth 
│   ├── compute_compartments.ipynb       # Compartment eigenvectors across depths
│   ├── compute_insulation.ipynb         # Insulation profiles and boundary calls across depths
│   ├── comp_insul_overview.ipynb        # Compartment & insulation metric families; biological floors 
│   └── other/                           # functions for compartment/insulation computing and other
│
├── raw/                                 # For raw .mcool files 
│
├── tables/								 # Table output folder
│
└── hic-tools.def                        # Singularity container definition
```

---

## Notebooks: Figures and Tables

| Notebook | Figures | Tables | Description |
|----------|---------|--------|-------------|
| `compute_slopes.ipynb` | — | — | Processes raw `.mcool` files into P(s) curves and log-derivatives |
| `loop_size_and_density.ipynb` | 2, S1 | S1, S2, S3 | Infers loop size and loop density |
| `compute_compartments.ipynb` | — | – | Compartment E1 across depths and correlation with full-depth references |
| `compute_insulation.ipynb` | — | – | Insulation profiles and boundary calls across depths |
| `comp_insul_overview.ipynb` | 3, 4, 5, S2, S3 | S4, S5 | Visualization of compartment and insulation metric families; floors |
