# Single-Cell Profiling Revealed the Unique Phenotypic Signatures of Malignant NK and T Cells and Immune Dynamics in Extranodal NK/T-Cell Lymphoma Before and After PD-1 Blockade

## Overview
Extranodal NK/T-cell lymphoma (ENKTCL) is an Epstein-Barr virus-associated malignancy with heterogeneous responses to immune checkpoint blockade. We performed longitudinal single-cell RNA sequencing with paired T-cell and B-cell receptor profiling in 11 patients treated with PD-1-based chemoimmunotherapy. Integrated transcriptomic and copy number analyses identified chromosome 6p21 amplification (chr6p21Amp) and deletion (chr6p21Del) as a major genomic axis stratifying malignant NK/T cells. Chr6p21Amp malignant NK cells corresponded to LMP-1 and PD-L1 positive classical states, whereas chr6p21Del cells exhibited reduced MHC class I gene expression and adopted stem-like programs. PD-1 therapy induces immune remodeling within the microenvironment, characterized by an expansion of CXCL13⁺ helper T cells, establishing a structured, long-term immunoactivated state. PD-1 therapy rapidly depleted malignant NK cells, whereas malignant T cells persisted with dynamic immune remodeling. Immunosuppression was sustained through PD-L1/PD-1 interactions and DPP4-associated chemokine signaling. These findings establish a mechanistic link between malignant genomic features and immune microenvironment remodeling following PD-1 blockade, guiding the rational design of next-generation immunotherapies for NKTCL. 

## User Guide

### Environmental Setup
This repository contains Jupyter notebooks (.ipynb) with all analysis code. Please ensure the following software versions are installed:
- Python: 3.11.11
- R: 4.3.3

**Python Dependencies**
Install all required Python packages using the provided `requirements.txt`:
```bash
pip install -r requirements.txt
```
**Data Download**
Processed data files are not included in the repository due to file size limits. Please use the provided shell script to download the data:
```bash
bash download.sh
```
This script will automatically download all processed data files into the `data/` directory. If you intend to run analyses on the external validation datasets, please also execute:
```bash
bash download-external.sh
```
### Running the Analysis

The Jupyter notebooks are organized into two categories:

1.  **Final Figures (Clean & Well-Documented)**
    - **File:** `Figure1-6.ipynb`
    - **Purpose:** This notebook is comprehensive and well-documented. It contains the finalized code, analysis, and plotting commands required to generate the main figures in the manuscript. For users wishing to understand the final methodology or reproduce the publication results, this is the primary file to review and execute.

2.  **Data Processing (Intermediate & Messy)**
    - **Files:** `01_NKTCL_basic.ipynb` to `06_NKTCL_monocle.ipynb`
    - **Purpose:** These notebooks contain the step-by-step raw data processing pipelines, quality control steps, and intermediate statistical analyses.
    - **Note:** The code in these notebooks is functional but may appear "messy" (i.e., contains verbose comments, debugging lines, and exploratory tests). They represent the initial data wrangling and analysis exploration phases. While they are provided for full transparency, they are less streamlined than the final figure notebook.