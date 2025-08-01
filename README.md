# General-Causal-Discovery

## General Causal Discovery: Extending LiNGAM via KCI
To overcome the difficulties associated with explicit residual estimation in nonlinear causal analysis, we propose General Causal Discovery. Our method introduces a new HSIC-type statistic, taking advantage of the asymmetry in causal relationships highlighted by LiNGAM, without relying on residual estimation. This enables its application to both the determination of causal direction and the inference of causal ordering.

## Real datasets

### Bivariate case

For the cause-effect scenario, we use the Tübingen Cause-Effect Pairs (TCEP) dataset, which consists of 108 cause-effect pairs. We exclude the multivariate pairs (52, 53, 54, 55, 71, 105). As a result, our experiments include a total of 102 cause-effect pairs.

As mentioned in the paper, since we use up to the first 500 data points for each cause-effect pair in our experiments rather than random sampling, readers can download the data from the following URL **https://webdav.tuebingen.mpg.de/cause-effect/** and reproduce our results.

### Multivariate case

We use the Sachs dataset, and we provide both its ground-truth DAG and the corresponding data.

We recommend using CAM pruning due to its high computational efficiency. However, for the Sachs dataset, the best results were achieved using KCI pruning with $\sigma_1 = 0.8$ and $\sigma_2 = 0.4$, and this finding was consistent across both the GCD_star and GCD models as GCD_dag_Sachs.csv.

### Required Softwares

R Language version 4.4.1

RStudio 2024.09.0+375
