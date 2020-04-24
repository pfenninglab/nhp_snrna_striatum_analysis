# nhp_snrna_striatum_analysis

The code and jupyternotebooks to go along with the Medium Spiny Neuron Diversity in the Primate Nucleus Accumbens publication

## Requirements 

* R>3.6.0

* python>3.5.1

* STAR>2.7.2b

The rest of the requirements can be installed by running:

1. pip install pip_requirements.txt
2. Rscript installation.R

## Computational Workflow

Alignment of the rhesus macaque straitum nuclei data was performed using [STAR version 2.7.2b](https://github.com/alexdobin/STAR) on the rhemac10 genome.
Scripts for the alignmet are located in the alignment scripts folder.

The analysis was composed of the following steps and can be viewed in the analysis jupyer notebooks.

| Computational Step   | Software Packages  | Associated Notebook  |
| ------------- |:-------------:| -----:|
| Empty Droplet Removal   | [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) | [remove_noise_droplets.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/remove_noise_droplets.ipynb) |
| Doublet Removal    | [SCDS]( https://bioconductor.org/packages/release/bioc/html/scds.html)   |   [remove_noise_droplets.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/remove_noise_droplets.ipynb) |
| Combining Brain Regions  | ----  |   [combine_brain_regions.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/combine_brain_regions.ipynb) |
| Poor Quality Nuclie Removal | [Scanpy](https://scanpy.readthedocs.io/en/stable/)     |   [preclustering_and_qc.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/preclustering_and_qc.ipynb) |
| Preclustering | [Scanpy](https://scanpy.readthedocs.io/en/stable/)    |    [preclustering_and_qc.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/preclustering_and_qc.ipynb) |
| Normalization | [Scran](https://bioconductor.org/packages/release/bioc/html/scran.html)     |    [scran_normalization.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/scran_normalization.ipynb) |
| Full Striatum Clustering and Marker analyis|  [Scanpy](https://scanpy.readthedocs.io/en/stable/)    |  [full_nuclei_analysis.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/full_nuclei_analysis.ipynb)   |
| Medium Spiny Neuron Clustering and Marker analyis|  [Scanpy](https://scanpy.readthedocs.io/en/stable/)    |  [msn_analysis.ipynb](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/blob/master/analysis/msn_analysis.ipynb)   |

## Statistical testing

We are also providing our own code for the following statistical tests 

| Statistical Test  | Purpose | Script  |
| ------------- |:-------------:| -----:|
|Cluster-cell type enrichment test | Computes the marker gene enrichment between a single nuclei cluster and set of cell types in a gene marker database | [cluster_profile.py](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/tree/master/src/cluster_profile.py) |
|Cluster split test | Computes a gene-wise split test of whether a cluster split splits a gene expression distribution based on its bimodality | [bimodal_split_test.py](https://github.com/pfenninglab/nhp_snrna_striatum_analysis/tree/master/src/bimodal_split_test.py) |




Data will be available via SRA and GEO soon!




