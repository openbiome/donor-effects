# Re-analysis of Paramsothy et al.

- `data.tsv` has the patient data
    - column `treatment` shows treatment
        - originally treated with FMT (`fmt`)
        - received FMT during rescue (`rescue`)
        - never recevied (`no_rescue`)
    - column `outcome` shows whether the patient achieved the primary outcome (0/1)
    - column `A` shows whether stool from donor A was present in the pool that that patient received
- `analyze.R` does the statistical tests
