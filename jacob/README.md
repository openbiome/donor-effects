# Re-analysis of Jacob et al.

- `data/data.tsv` records the patient outcomes and the material they received
    - column `pool` is an identifier for each unique combination of donors
    - column `donorA` shows (0/1) whether the pool of stool received by that patient included donor A
- `data/raw/` contains the raw sequencing data
    - `Snakefile` downloads, processes, and analyses the 16S data
    - `make-manifest.R` gets the Qiime 2 manifest file
- `analyze.R` has the basic statistical tests
- `diversity-model.R` analyzes the patient outcomes by pool bacterial diversity
