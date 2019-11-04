# Re-analysis of Goyal et al.

- `data/` contains the raw sequencing data
    - `Snakefile` downloads, processes, and analyses the 16S data
    - `make-manifest.R` gets the Qiime 2 manifest file
- `analysis.R` performs the analysis on `data.tsv` and `patients.tsv`
- `clean.R` cleans that data from the `data/` folder
- `data.tsv` is the cleaned data made by `clean.R`
- `patients.tsv` records the outcomes (0/1) of each patient
