# Detecting Distribution Shifts in Federated Learning via Topological Data Analysis

MSc Statistics dissertation, UCL.

Detecting non-IID distribution shift across hospitals ahead of federated learning, using persistence diagrams as a privacy-preserving alternative to raw patient data. This project benchmarks four TDA-based hypothesis tests against a raw-covariate baseline, on both synthetic hospital data and real TCGA cancer data.

## Motivation

In federated learning across hospitals, checking whether client data distributions are non-IID typically requires access to raw patient covariates, which may not be shareable for privacy reasons. This project investigates whether **persistence diagrams**, computed locally at each hospital, can be shared instead and still reveal meaningful distribution shift, without exposing individual patient records.

## Methods compared

- **Inco-variance test** (Krebs & Rademacher): a two-sample test based on the variance of pairwise Wasserstein distances between persistence diagrams.
- **Permutation test on the F-loss statistic** (Robinson & Turner): a permutation test comparing within-group vs. pooled dispersion of persistence diagrams.
- **Persistence landscape permutation test**: a permutation test on the integrated squared distance between groups' mean persistence landscapes.
- **Persistence landscape z-test** : a two-sample z-test on total landscape mass (area under the landscape).
- **Raw-covariate Wasserstein permutation test**: a non-topological baseline run directly on patient covariates, for comparison against the TDA-based tests.

## Repository contents

| File | Description |
|---|---|
| `TDA_Mild_H1_experiments_and_assumption_checks.ipynb` | Synthetic two-hospital experiment with a mild H1 (loop/topological) signal. Runs all four TDA tests plus the raw-covariate baseline, and checks core theoretical assumptions (Hoeffding decomposition, U-statistic variance convergence, third- and fourth-moment convergence). |
| `TDA_Strong_H1_experiments_and_assumption_checks.ipynb` | Same experimental setup as above, but with a strong, well-separated ring/H1 topological signal built into the synthetic data. |
| `power_and_size_analysis_permutation_tests.ipynb` | Power and size (Type I error) analysis for the permutation-based tests: sweeps over effect size (alpha), number of permutations, subsample size, and number of subsamples, and compares the R&T, landscape, and raw-Wasserstein tests directly. |
| `Sphere_and_Torus_Outlier_Analysis.ipynb` | Illustrative examples of persistence diagrams and barcodes on canonical shapes (sphere vs. torus), including sensitivity checks under added Gaussian noise and outlier contamination. |
| `TCGA_Exploratory_Data_Analysis.R` | Exploratory data analysis of the TCGA patient covariates dataset: cleans inconsistent hospital ("site") names, pools patients into the largest hospitals, and visualises demographic, stage, and cancer-type composition by hospital. |

## Data

- **Synthetic data**: two simulated hospitals ("Hospital A" and "Hospital B") with configurable demographic distributions (age, severity, gender) and, in the strong-H1 notebook, an embedded ring/loop structure to create a genuine topological signal.
- **Real data**: patient covariates from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga), grouped by contributing hospital/site.

## Requirements

Notebooks are designed to run in Google Colab. Key Python packages:
- `ripser`, `persim`, `gudhi` — persistent homology and persistence landscapes
- `POT` (`ot`) — optimal transport / Wasserstein distances
- `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `scipy`

The R script requires:
- `ggplot2`

## Author

Beatriz — MSc Statistics, UCL, with the Alan Turing Institute.
