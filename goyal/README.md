# Re-analysis of Goyal et al.

- `diversity-data/` contains the raw sequencing data
    - `make-manifest.R` gets the Qiime 2 manifest file
- `analysis.R` performs the analysis on `data.tsv` and `patients.tsv`
- `patient-data.tsv` records the outcomes (0/1) of each patient

In two cases (patients #24 and #29), the donor has two associated samples. Here
we used just the first sample (`D1`).
