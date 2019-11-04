# Re-analysis of Kump et al.

- `data/patients.tsv` has the patient outcomes
- `data/donations.tsv` has the donor, patient, timepoint, and diversity of the donor stool sample
- `data/raw/` contains the raw sequencing data
    - `Snakefile` downloads, processes, and analyses the 16S data
    - `make-manifest.R` gets the Qiime 2 manifest file
- `analyze.R` has the basic statistical tests and the analysis of donor diversities
