# Re-evaluating the evidence for fecal microbiota transplantation "super-donors" in inflammatory bowel diseases

Author: Scott Olesen <solesen@openbiome.org>

## Files

- `README.md`: This file
- `donor-effects-env.txt`: Conda environment file
- `install-packages.R`: R script that installs packages needed for analyses
- `Snakefile`: A [Snakemake](https://snakemake.readthedocs.io/) file with instructions for how to repeat the analyses
- `16S-common.snakefile`: A Snakemake file with instructions for 16S rRNA data processing
- `16S-single.snakefile`: 16S processing for single-ended reads (Goyal *et al*.)
- `16S-paired.snakefile`: 16S processing for paired-ended reads (Jacob *et al*. and Kump *et al*.)
- `costello/`: Re-analysis of Costello *et al*.
- `goyal/`: Re-analysis of Goyal *et al*.
- `jacob/`: Re-analysis of Jacob *et al*.
- `kump/`: Re-analysis of Kump *et al*.
- `moayyedi/`: Re-analysis of Moayyedi *et al*.
- `nishida/`: Re-analysis of Nishida *et al*.
- `nusbaum/`: Re-analysis of Nusbaum *et al*.
- `paramsothy/`: Re-analysis of Paramsothy *et al*.
- `pools/`: Analysis of the 3 "pool stool" studies
- `rossen/`: Re-analysis of Rossen *et al*.
- `utils.R`: Utilities for R scripts
- `uygun/`: Re-analysis of Uygun *et al*.
- `vaughn/`: Re-analysis of Vaughn *et al*.

Goyal, Jacob, Kump, and Nusbaum have subfolders `diversity-data/`, which are
the workspaces for the 16S Snakefiles. The alpha- and beta-diversity results
are placed there.

## Getting started

1. Install [conda](https://docs.conda.io/)
1. Install the working environment (Qiime and Snakemake): `conda create --name donor-effects --file env.txt`
1. Install R packages using `./install-packages.R` into that environment
1. Run `snakemake` in that environment
