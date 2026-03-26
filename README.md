## FDPCox: Cox regression under federated differential privacy

This repository contains the `R` package `FDPCox`, which implements the methods from the preprint [Optimal Cox regression under federated differential privacy: coefficients and cumulative hazards](https://arxiv.org/abs/2508.19640). After installing the `FDPCox` package, running the following scripts from the `inst` directory reproduces these experiments from the paper: 

- `CDP_sims.R`: central DP simulations (varying sample size, privacy budget and dimension) [Figure 1], and label central DP simulations [Figure 9]
- `FDP_sims.R`: federated DP simulations for sequentially-interactive and fully-interactive mechanisms [Figures 2 and 3]
- `real_data.R`: demonstration of our methods on the `survival::rotterdam` dataset [Figure 4]
- `sensitivity_sims.R`: central DP simulations showing sensitivity to the covariate bound and descent step-size [Figure 5]
- `censoring_sims.R`: simulations varying the rate of the censoring distribution [Figure 6]
- `CDP_numsims.R`: log-log plots for numerical demonstration of phase transitions [Figures 7 and 8]

Note: `tidyverse` is needed to produce the plots. In addition, `real_data.R` uses the `survival` and `abind` packages. 
