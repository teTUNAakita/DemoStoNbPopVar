# DemoStoNbPopVar

Code supporting the manuscript:

**Demographic stochasticity can drive population variability even in large populations**

Author: Tetsuya Akita  
Journal:  
DOI:  

## Overview

This repository contains the R code used to reproduce the main simulation figures and the supplementary empirical re-analysis presented in the manuscript.

The study develops an age-structured framework linking the effective number of breeders (`Nb`) to temporal variability in mature female abundance. The simulation scripts reproduce the population trajectories and the relationships between `Nb`, recruitment variability, and mature-population variability. An additional script re-analyzes publicly available data from Clarke et al. (2024) to examine the empirical co-occurrence of low effective-to-census ratios and small absolute `Nb`.

All scripts are intended to be run from the repository root.

## Repository structure

```text
DemoStoNbPopVar/
├── README.md
├── FigS1_BevertonHolt.R
├── Fig2_PopulationTrajectories.R
├── Fig3_FigS2_Nb_Variability.R
├── clarke2024_nb_nc_reanalysis.R
└── figure/
```

The `figure/` directory is created automatically by the figure-generating scripts if it does not already exist.

## Requirements

The simulation and figure scripts use base R only.

The Clarke et al. re-analysis additionally requires:

```r
library(dplyr)
```

The code was written for R and is expected to run on recent R versions. No platform-specific paths are used.

## Scripts

### `FigS1_BevertonHolt.R`

Reproduces **Figure S1**, showing the Beverton-Holt recruitment function and the replacement line used in the age-structured simulations.

The script:

- defines the survival, maturation, and body-weight schedules;
- calculates spawning biomass per recruit (`SPR`);
- derives the deterministic equilibrium;
- plots the Beverton-Holt recruitment curve, replacement line, and their equilibrium intersection.

Output:

```text
figure/FigS1_BevertonHolt.pdf
```

### `Fig2_PopulationTrajectories.R`

Reproduces the four panels of **Figure 2**, showing representative trajectories of mature female abundance under contrasting reproductive overdispersion and environmental variability.

The four scenarios are:

| Panel | `phi` | `sigma_E` |
|---|---:|---:|
| a | `1e2` | `0.0` |
| b | `1e2` | `0.6` |
| c | `1e-5` | `0.0` |
| d | `1e-5` | `0.6` |

The model includes female age structure, age-specific maturation, binomial survival, age-dependent reproductive output, Beverton-Holt mean recruitment, Poisson-lognormal recruitment variation, and the effective number of female breeders (`Nb`).

Each panel is saved separately.

Outputs:

```text
figure/Fig2a_phi100_sigmaE0.pdf
figure/Fig2b_phi100_sigmaE06.pdf
figure/Fig2c_phi1e-5_sigmaE0.pdf
figure/Fig2d_phi1e-5_sigmaE06.pdf
```

### `Fig3_FigS2_Nb_Variability.R`

Reproduces **Figure 3** and **Figure S2**.

The script runs stochastic age-structured simulations across 13 values of reproductive overdispersion (`phi`), three levels of environmental variability (`sigma_E = 0, 0.3, 0.6`), and 20 replicate simulations per parameter combination. Each simulation runs for 10,000 years, with the first 500 years discarded as burn-in.

The script calculates mean effective number of female breeders (`Nb`), recruitment variability (`CV[R]`), and mature female population variability (`CV[N]`), and compares the simulated results with the analytical approximations

```text
CV[R] ≈ sqrt{ exp(sigma_E^2) [1 + 1/(Nb - 1)] - 1 }
```

and

```text
CV[N] ≈ gamma * CV[R]
```

where `gamma` is determined by the survival and maturation schedules.

Outputs:

```text
figure/Fig3_Nb_vs_CVN.pdf
figure/FigS2_Nb_vs_CVR.pdf
```

### `clarke2024_nb_nc_reanalysis.R`

Re-analyzes the non-replicated dataset associated with Clarke et al. (2024) to examine whether low `Nb/Nc` ratios tend to coincide with small absolute `Nb`.

The required input file is:

```text
MEC-23-0962_NonReplicated.csv
```

It can be downloaded from Dryad:

https://datadryad.org/dataset/doi:10.5061/dryad.p2ngf1vzm

Dryad DOI:

https://doi.org/10.5061/dryad.p2ngf1vzm

Place the CSV file in the repository root before running the script.

The analysis:

- retains records with both `Nb` and census size (`Nc`);
- excludes censored `Nc` values for which an exact ratio cannot be calculated;
- recalculates `Nb/Nc` directly from the reported `Nb` and `Nc`;
- summarizes the frequency of small absolute `Nb` among populations with low `Nb/Nc`;
- examines alternative `Nb/Nc` thresholds;
- summarizes the taxonomic composition of low-ratio populations.

The thresholds `Nb = 50`, `100`, and `200` are used only as descriptive reference values and are not interpreted as universal biological thresholds.

The empirical `Nb` estimates in the Clarke et al. dataset generally combine both sexes, whereas `Nb` in the theoretical model is defined as the effective number of female breeders. Likewise, `Nc` in the comparative dataset is not necessarily identical to the census number of mature females (`N`) used in the theoretical model. The re-analysis should therefore be interpreted as an empirical analogue rather than a direct quantitative validation of the theory.

## Reproducing the analyses

From the repository root, the figure scripts can be run independently:

```r
source("FigS1_BevertonHolt.R")
source("Fig2_PopulationTrajectories.R")
source("Fig3_FigS2_Nb_Variability.R")
```

After downloading `MEC-23-0962_NonReplicated.csv` from Dryad, the supplementary re-analysis can be run with:

```r
source("clarke2024_nb_nc_reanalysis.R")
```

## Data availability

No new empirical data were generated in this study.

The supplementary empirical re-analysis uses the publicly available dataset associated with Clarke et al. (2024), available from Dryad:

https://doi.org/10.5061/dryad.p2ngf1vzm

## Code availability

All code required to reproduce the simulations, analyses, and figures described above is provided in this repository.

A permanent archived version of the repository will be created upon acceptance of the manuscript.

## Citation

If you use this code, please cite the associated manuscript:

Akita, T.  
**Demographic stochasticity can drive population variability even in large populations.**  
Journal:  
DOI:  

For the empirical re-analysis, please also cite the Clarke et al. (2024) dataset and associated publication.
