# SCR-macrophages

This repository contains scripts to run two-state modeling associated with the SCR Macrophage project. 

Requirements: 
1. `phago_suma_combined.Rds`: This is the only dataset required to run the two-state model and must be downloaded from [Zenodo](10.5281/zenodo.19699424). During peer review, this dataset will be restricted for reviewer access only.
2. `run_two-state_model_v3.Rmd`: This is the only script required to run the two-state model. This script is best run in Rstudio or another UI optimized for running R markdown files.

---



## Overview

Here, we describe the two-state model used to predict how SCR motif combinations influence macrophage phenotype. The model treats each phenotype measurement (surface marker MFI or phagocytosis TFI) as the output of a population of cells distributed between two states — an OFF state and an ON state. The fraction of cells in each state is governed by an equilibrium constant $K$, which depends on the identity and copy number of motifs encoded in the SCR.

Phenotype measurements for six readouts (CD163, CD80, CD206, CD40, PDL1, phagocytosis) were each modeled independently. The same mathematical framework and fitting procedure was applied to every phenotype.

---

## Assumptions

1. **Two-state population.** Each cell is either in an OFF (low) or ON (high) state. The observed fluorescence value for a sample reflects the weighted average of the two-state outputs across the cell population.

2. **Additive motif contributions.** Each motif contributes independently and additively to the total energy difference $\Delta E$ between states. There are no pairwise interaction terms between motifs.

3. **Fixed bounds.** The lower bound $f_{low}$ and upper bound $f_{high}$ of each phenotype are fixed as the minimum and maximum mean observed values across all samples in the dataset. These anchors define the full dynamic range of the model.

4. **Copy-number scaling.** A motif present in $n$ copies contributes $n \times \Delta E_s$ to the total energy, i.e., contributions scale linearly with copy number.

---

## Model Equations

### Two-State Equilibrium

Cells are in either an OFF state or an ON state. The equilibrium constant $K$ is defined as:

> ```math
> K = \frac{[ON]}{[OFF]} = \frac{[ON]}{1 - [ON]}
> ```

which can be rearranged to give the fraction of cells in the ON state:

> ```math
> [ON] = \frac{K}{1 + K}
> ```
> <div align="right"><em>(Eq. 1)</em></div>

### Observed Fluorescence

The observed fluorescence of a sample is a weighted average of the fluorescence produced by cells in each state:

> ```math
> f_{obs} = f_{low}(1 - [ON]) + f_{high}[ON]
> ```
> <div align="right"><em>(Eq. 2)</em></div>

where $f_{low}$ and $f_{high}$ are the lowest and highest mean observed fluorescence values for that phenotype across the entire dataset.

### Relating K to Energy

The equilibrium constant is related to the energy difference $\Delta E$ between states:

> ```math
> K = e^{-\Delta E}
> ```
> <div align="right"><em>(Eq. 3)</em></div>

Substituting Eq. 3 into Eq. 1:

> ```math
> [ON] = \frac{1}{1 + e^{\Delta E}}
> ```
> <div align="right"><em>(Eq. 4)</em></div>

### Additive Motif Model for $\Delta E$

The total $\Delta E$ for a given SCR is the sum of an intrinsic baseline and per-copy contributions from each motif:

> ```math
> \Delta E = \Delta E_i + \sum_{j=1}^{9} \Delta E_{s,j} \cdot n_j
> ```
> <div align="right"><em>(Eq. 5)</em></div>

where $\Delta E_i$ is the intrinsic energy when no motifs are present, $\Delta E_{s,j}$ is the per-copy contribution of motif $j$, and $n_j$ is the copy number of motif $j$ in the SCR.

### Full Prediction Equation

Combining Eqs. 2, 4, and 5 gives the single prediction equation used in fitting:

> ```math
> f_{pred} = f_{low} + \frac{f_{high} - f_{low}}{1 + \exp\left(\Delta E_i + \sum_{j=1}^{9} \Delta E_{s,j} \cdot n_j\right)}
> ```
> <div align="right"><em>(Eq. 6)</em></div>

---

## Parameter Estimation

### Fitting Data

Parameters were estimated using all individual technical replicates.

Phenotype-specific bounds $f_{low}$ and $f_{high}$ were computed from the single minimum and maximum mean value across all samples and held fixed throughout fitting.

### Objective Function

For each phenotype, the 10 model parameters ($\Delta E_i$ and $\Delta E_{s,j}$ for $j = 1 \ldots 9$) were estimated by minimizing the sum of squared residuals in fluorescence space over all individual replicate measurements:

> ```math
> \mathcal{L} = \sum_{i} \left( f_{obs,i} - f_{pred,i} \right)^2
> ```
> <div align="right"><em>(Eq. 7)</em></div>

where $f_{pred,i}$ is evaluated using Eq. 6 with the motif copy numbers $n_{j,i}$ for replicate $i$.

### Optimization Algorithm

Minimization of Eq. 7 was performed using the BFGS quasi-Newton algorithm (`optim()` in base R, `method = "BFGS"`, `maxit = 5000`). Parameters were initialized with $\Delta E_i$ set to the mean observed $\Delta E$ of no-motif (SCR0) replicates (computed by applying Eq. 4 to the no-motif data), and all $\Delta E_{s,j}$ initialized at zero.

---

## Model Evaluation

### In-Sample Fit

After fitting, predicted fluorescence values were computed by applying Eq. 6 to replicate-averaged data. Pearson $r$ and $R^2$ between replicate-averaged observed and predicted fluorescence values were used as in-sample goodness-of-fit metrics.

### Cross-Validation

Generalization was assessed by 10-fold cross-validation applied to the individual replicate data. In each fold, the model was re-fit on 90% of replicates using the same BFGS procedure.

Two metrics are reported from cross-validation:

- **CV-R²**: pooled held-out $R^2$ in fluorescence space, as a diagnostic for overfitting.
- **CV-SD**: the standard deviation of each $\Delta E_{s,j}$ estimate across the 10 folds, used as the error bar on $\log_{10}(K_s)$ plots.

---

## Converting $\Delta E$ to Equilibrium Constants

Each fitted $\Delta E_{s,j}$ is converted to a per-motif equilibrium constant contribution:

> ```math
> K_{s,j} = e^{-\Delta E_{s,j}}
> ```
> <div align="right"><em>(Eq. 8)</em></div>

The intrinsic equilibrium constant (when no motifs are present) is:

> ```math
> K_i = e^{-\Delta E_i}
> ```
> <div align="right"><em>(Eq. 9)</em></div>

**Interpretation:**

| $K_s$ | $\Delta E_s$ | Meaning |
|:---:|:---:|:---|
| $> 1$ | $< 0$ | Motif shifts equilibrium toward ON (high) state |
| $= 1$ | $= 0$ | Motif has no effect on equilibrium |
| $< 1$ | $> 0$ | Motif shifts equilibrium toward OFF (low) state |

$K_s$ values are visualized on a $\log_{10}$ scale. Uncertainty bars reflect ±1 CV-fold SD, propagated from $\Delta E$ space to $\log_{10}(K_s)$ space as:

> ```math
> \text{SD}[\log_{10}(K_s)] = \log_{10}(e) \times \text{SD}[\Delta E_s]
> ```
> <div align="right"><em>(Eq. 10)</em></div>

---

## Predicting Phenotypes for Novel Combinations

The total predicted $K$ for any SCR with motif copy numbers $n_1, \ldots, n_9$ is:

> ```math
> K_{pred} = K_i \prod_{j=1}^{9} K_{s,j}^{n_j}
> ```
> <div align="right"><em>(Eq. 11)</em></div>

The predicted fraction of cells in the ON state is:

> ```math
> [ON]_{pred} = \frac{K_{pred}}{1 + K_{pred}}
> ```
> <div align="right"><em>(Eq. 12)</em></div>

And the predicted fluorescence is:

> ```math
> f_{pred} = f_{low} + [ON]_{pred} \cdot (f_{high} - f_{low})
> ```
> <div align="right"><em>(Eq. 13)</em></div>

Eqs. 11–13 were applied to all 715 possible SCR designs with a total of 0–4 motif copies ($n_j \geq 0$, $\sum_j n_j \leq 4$) to generate the predicted design space shown in Figure 6.

Note that Eq. 11 is equivalent to Eq. 6 — the product-of-$K_s$ form is simply an algebraically equivalent rearrangement that makes the multiplicative contribution of each motif explicit.

---

## Software

All analysis was performed in R (v.4.3.2). Nonlinear optimization used `optim()` from base R. Cross-validation splits were generated with the `rsample` package. Data manipulation used `dplyr`, `tidyr`, `purrr`, and `stringr`. Plots were produced with `ggplot2`, `ggrepel`, and `ggh4x`. See `run_two-state_model_v3.Rmd` for the full analysis script.
