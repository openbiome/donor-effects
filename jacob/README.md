# Re-analysis of Jacob et al.

## Files

- `patient-data.tsv` records the patient outcomes and the material they received
    - column `pool` is an identifier for each unique combination of donors
    - column `donor1` shows (0/1) whether the pool of stool received by that patient included donor 1
- `diversity-data/` contains the raw sequencing data
    - `make-manifest.R` gets the Qiime 2 manifest file
- `analyze.R` has the basic statistical tests

## Notes

- Clinical response at 4 weeks was used as the outcome
