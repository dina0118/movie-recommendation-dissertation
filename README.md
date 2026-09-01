# Ordinal Collaborative Filtering for Movie Recommendation

Code for an MSc dissertation testing whether per-user category thresholds in an
ordinal collaborative filtering model improve the calibration of predicted
rating distributions. Six models are compared on MovieLens-32M, from a Gaussian
matrix factorisation baseline through to jointly trained factors with per-user
thresholds, with a further variant adding film content features.

## Requirements

R 4.4 or later. Packages: `data.table`, `ggplot2`, `recosystem`, `ordinal`,
`patchwork`, `scales`.

Data is not included in the repository; see `data/README.md` for download
instructions.

## Structure

| Folder | Contents |
|---|---|
| `data/` | Study period, training window search, eligibility rules |
| `models/` | Baseline factorisation, gradient check, joint training |
| `content/` | TMDb feature construction and content-feature models |
| `analysis/` | Metrics, calibration, decomposition, figures |

## Running order

Scripts are numbered and should be run in that order. Files in `models/` write
cached fits that later scripts read, so `07`, `09` and `10` must complete before
anything in `analysis/`.

1. `data/00_data_selection_eda.R` — selects the training window and applies the
   eligibility rules
2. `models/07_M1_M2_first_run_v2.Rmd` — baseline and frozen-factor ordinal models
3. `models/08_M2d_gradient_check.R` — verifies the analytic gradient before any
   joint model is fitted
4. `models/09_M2d_training.Rmd`, `models/10_M2e_training.Rmd` — joint training
5. `content/` — feature construction, then the content-feature models
6. `analysis/` — evaluation and figures, starting with
   `11_M2e_analysis_and_figures.Rmd`

Several scripts in `analysis/` expect objects created by `11`, and expect the
full six-model probability list to be in memory. Running a content-feature
script overwrites that list, so run the content scripts last.

## Notes

The joint models are implemented directly from the equations in Koren and Sill
(2011), since no open-source implementation exists and no R package combines
ordinal likelihoods with matrix factorisation. Correctness is checked three
ways: against central differences, against `ordinal::clm` with the factors held
fixed, and by recovering known parameters from simulated data.
