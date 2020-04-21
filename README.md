# Re-evaluating the evidence for fecal microbiota transplantation "super-donors" in inflammatory bowel diseases

Author: Scott Olesen <scott@openbiome.org>

## Files

- `README.md`: This file
- `install_packages.R`: R script that installs packages needed for analyses
- `Snakefile`: A [Snakemake](https://snakemake.readthedocs.io/en/stable/) file with instructions for how to repeat the analyses
- `16S_common.snakefile`: A Snakemake file with instructions for 16S rRNA data processing
- `16S_single.snakefile`: 16S processing for single-ended reads (Goyal *et al*.)
- `16S_paired.snakefile`: 16S processing for paired-ended reads (Jacob *et al*. and Kump *et al*.)
- `costello/`: Re-analysis of Costello *et al*.
- `goyal/`: Re-analysis of Goyal *et al*.
- `jacob/`: Re-analysis of Jacob *et al*.
- `kump/`: Re-analysis of Kump *et al*.
- `moayyedi/`: Re-analysis of Moayyedi *et al*.
- `nishida/`: Re-analysis of Nishida *et al*.
- `paramsothy/`: Re-analysis of Paramsothy *et al*.
- `pools/`: Analysis of the 3 "pool stool" studies
- `rossen/`: Re-analysis of Rossen *et al*.
- `utils.R`: Utilities for R scripts

## Getting started

1. Install Qiime 2 (version 2019.10) and active the environment ([instructions](https://docs.qiime2.org/2019.10/install/native/))
2. Install R packages using `./install_packages.R`
3. Run `snakemake`
