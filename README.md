# Integrated Kidney Human Cell Atlas

This repository contains the integration and analysis notebooks and results for the Human Kidney Cell Atlas (HKCA).

## Introduction

Despite numerous single-cell and single-nucleus transcriptomics studies, differences in cohorts, technologies, analytical pipelines, and annotation frameworks have limited the ability to define consensus kidney cell states, identify disease-associated populations and interpret kidney disease genetic susceptibility. Here, we assembled 18 human kidney single-cell and single-nucleus RNA-sequencing datasets spanning 232 donors and five major disease contexts into a uniformly processed and computationally integrated Human Kidney Cell Atlas (HKCA), comprising over one million high-quality cells (816,895) and nuclei (215,308). The HKCA resolves 63 cell types and 120 harmonized cell states, including rare epithelial and stromal populations associated with kidney diseases . Integration with spatial transcriptomics, intercellular communication networks, and human genetic association data further defined the anatomical context and disease relevance of these populations.


<img width="2291" height="1368" alt="Asset 2" src="https://github.com/user-attachments/assets/a8b49782-0c6b-48c8-8859-6ff248a19b65" />

## Data portal

`index.html` is the landing page for the atlas (abstract, data downloads, interactive embeddings):
https://zhou119010475.github.io/Kidney_integration/

## Interactive 3D UMAPs

Hosted via GitHub Pages (raw links in this repo won't render interactively — GitHub serves `.html` files as plain text):

- [HKCAsc_3DUmap](https://zhou119010475.github.io/Kidney_integration/HKCAsc_3DUmap.html) — single-cell 3D UMAP
- [HKCAsn_3DUmap](https://zhou119010475.github.io/Kidney_integration/HKCAsn_3DUmap.html) — single-nucleus 3D UMAP

## Models for automated annotation (Celltypist or scArches)
Find the trained models of the HKCA to use with your own dataset at https://zenodo.org/records/21468481

## Try the CellorNuc function to evaluate dissociation biases for single cell or nucleus RNAseq datasets:
  Example notebook: [notebooks/CellorNuc.ipynb](https://github.com/Peng-He-Lab/Kidney_integration/blob/main/notebooks/CellorNuc.ipynb)
