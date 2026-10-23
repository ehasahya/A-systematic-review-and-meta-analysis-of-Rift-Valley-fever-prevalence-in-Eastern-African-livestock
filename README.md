# R scripts

R code to reproduce the meta-analysis of Rift Valley fever (RVF) prevalence in cattle and small ruminants in Eastern Africa, 2014–2024. Each script runs on its own, reads one or more Excel files from the `data/` folder, and writes its figures and tables to an output folder created next to the input data.

## Scripts

| Script | What it does | Input file(s) | Manuscript output |
|---|---|---|---|
| `RVF_meta_prevalence_pooling_analysis (cattle & small ruminants) 2.1.R` | Pools prevalence for each species and assay (RT-PCR, IgM, IgG); forest plots with 95% CI, 95% prediction interval, I² and τ² | `meta RVF overall.xlsx` | Figures 2–4; Supplementary Table 4 |
| `RVF_meta_prevalence_publication-bias_analyses (SUPPLEMENT, exploratory only) 2.0.R` | Funnel plots and Egger's regression test for each species and assay (exploratory) | `meta RVF overall.xlsx` | Supplementary Figure 1 |
| `country_subgp_analysis_2.0.R` | Subgroup analysis by country, with between-country tests | `country subgroup.xlsx` | Table 1 (country) |
| `symptoms_subgp_analysis_2.0.R` | Subgroup analysis by symptom status (symptomatic vs asymptomatic) | `symptoms.xlsx` | Table 1 (symptoms) |
| `RVF_ecological_primary_SUBGROUP_analyses_2.0.R` | IgG subgroup analyses by five ecological and epidemiological factors derived from Google Earth Engine and FAO EMPRES-i | `floods.xlsx`, `near border.xlsx`, `near national parks.xlsx`, `outbreaks by empres.xlsx`, `urban v9.xlsx` | Table 2 |
| `RVF_ecological_SUBGROUP_analyses - REPORTED-ONLY SENSITIVITY 2.0.R` | The same five IgG subgroup analyses using only information reported in the primary studies | `floods reported.xlsx`, `border reported.xlsx`, `parks reported.xlsx`, `outbreaks reported.xlsx`, `urban or rural reported.xlsx` | Supplementary Table 6 |

## Requirements


- R packages: metafor, clubSandwich, dplyr, readxl, writexl

```r
install.packages(c("metafor", "clubSandwich", "dplyr", "readxl", "writexl"))
```


## How to run

1. Download the repository, including the data folder.
2. Open a script and set the path at the top to your local copy of its input:
   - `in_path` in the pooling and publication-bias scripts
   - `IN_PATH` in the country and symptom scripts
   - `DATA_DIR` (folder holding the five factor files) and `MAIN_PATH` in the two ecological scripts
3. Run the whole script. The scripts are independent and can be run in any order.

## Required columns in the input files

| Script | Columns |
|---|---|
| Pooling, publication bias | `Author`, `Species`, `Test`, `Cases`, `SampleSize` |
| Country | `Author`, `Species`, `Country`, `Test`, `Cases`, `SampleSize` |
| Symptoms | `Author`, `Species`, `Symptoms` (`yes`/`no`), `Test`, `Cases`, `SampleSize` |
| Ecological, primary | `Author`, `Species`, `Cases`, `SampleSize`, `Test` (optional), and the factor column: `Floods`, `Border`, `National_park`, `Outbreak` or `Urban_or_rural` |
| Ecological, sensitivity | as above, with factor column `floods_reported`, `Border_reported`, `Parks_reported`, `outbreak_reported` or `urban_or_rural_reported` |

`Species` must be `cattle` or `small ruminants`, and `Test` must be `RT-PCR`, `IgM` or `IgG`. `Author` identifies the study (e.g., `Bett et al (2018)`) and must be identical for all rows from the same study, because it is used to group estimates within studies.

## Outputs

| Script | Output folder | Files |
|---|---|---|
| Pooling | `R3_pooling_results_v11/` | `forest_<species>_<assay>.png`; `R3_pooled_results.xlsx` (pooled estimates); `R3_estimates_per_study.xlsx` |
| Publication bias | `pubbias_supplement_v4/` | `funnel_<species>_<assay>.png`; `S_publication_bias.xlsx` |
| Country | `country subgroup output/` | `country_<species>_<assay>.png`; `country_level_estimates.xlsx`; `country_between_tests.xlsx` |
| Symptoms | `symptom_R3_output/` | `symptom_<species>_<assay>.png`; `symptom_level_estimates.xlsx`; `symptom_between_tests.xlsx` |
| Ecological, primary | `subgroup_multilevel_v3/` | `forest_<factor>_<species>.png`; `subgroup_level_estimates.xlsx`; `subgroup_between_tests.xlsx`; `reconciliation.xlsx` |
| Ecological, sensitivity | `out/` | `<factor>_<species>.png`; `level_estimates.xlsx`; `between_tests.xlsx`; `reconciliation.xlsx` |

## Methods implemented

- Prevalence is pooled on the logit scale with a 0.5 continuity correction for zero events, using REML.
- A random-effects model is fitted when each study contributes one estimate. A three-level model (estimates nested within studies) is fitted when a study contributes more than one, as in the subgroup analyses.
- Between-subgroup differences are tested with the subgroup as a moderator, and the test is repeated with cluster-robust variance estimation (CR2), a binomial-normal GLMM and the Freeman–Tukey transformation.
- Egger's test is exploratory, i.e., it uses one estimate per row, cannot account for clustering within studies, and has low power when there are few studies and heterogeneity is high.

