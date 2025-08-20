# General-Causal-Discovery

## General Causal Discovery: Extending LiNGAM via KCI
To overcome the difficulties associated with explicit residual estimation in nonlinear causal analysis, we propose General Causal Discovery. Our method introduces a new HSIC-type statistic, taking advantage of the asymmetry in causal relationships highlighted by LiNGAM, without relying on residual estimation. This enables its application to both the determination of causal direction and the inference of causal ordering. In GCD.R, GCD_star and GCD represent the deflation and no-deflation versions of the code, respectively.

## Real datasets

### Bivariate case

For the cause-effect scenario, we use the Tübingen Cause-Effect Pairs (TCEP) dataset, which consists of 108 cause-effect pairs. We exclude the multivariate pairs (52, 53, 54, 55, 71, 105). As a result, our experiments include a total of 102 cause-effect pairs.

As mentioned in the paper, since we use up to the first 500 data points for each cause-effect pair in our experiments rather than random sampling, readers can download the data from the following URL **https://webdav.tuebingen.mpg.de/cause-effect/** and reproduce our results. For this real dataset, we provide a preprocessing script to help users obtain the corresponding effect variable across all datasets.

### Multivariate case

We use the Sachs dataset, and we provide both its ground-truth DAG and the corresponding data.

We recommend using CAM pruning due to its high computational efficiency. After obtaining the estimated DAG, please use the Metrics.R file to evaluate the experimental results.

### Required Softwares

R Language version 4.4.1

RStudio 2024.09.0+375
