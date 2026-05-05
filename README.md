---
title: "AFTCoop"
output: rmarkdown::github_document
---

# AFTCoop

**AFTCoop** is an R package for fitting cooperative AFT survival regression models with two omics views.

## Description

The **AFTCoop** package estimates cooperative AFT survival regression models that allow the integration of two matrices of covariates, **U** and **Z**, to improve the estimate.
It implements the direct cooperative algorithm proposed in Angelini, De Canditiis, De Feis, and Iuliano (Submitted 2025).
It relies on adding a cooperative term to a penalized negative log-likelihood function to enforce agreement between two views representing multi-omics datasets associated with the same patients, such as gene expression and methylation.
Cooperative learning is a novel approach to multi-view data integration, initially proposed by Ding et al. (2022).
In this setting, the AFT negative log-likelihood and LASSO penalization are combined to enforce sparsity, and an agreement penalty is used to encourage agreement between predictions from the two data views while improving survival estimates.
The user can choose between three models: *Weibull*, *lognormal*, or *log-logistic*, which handle different noise distributions in the log-linear AFT representation.
The regression coefficient estimates are obtained using the proximal gradient method, and the regularization parameter is chosen using cross-validation.

Current released version is AFTCoop 0.2.4, see NEWS for the changes.

## 🧪 Installation

For now, install the latest version of **AFTCoop** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("angeclau/AFTCoop")
```

In alternative, since  `install_github()` was deprecated in devtools 2.5., you can alse install **AFTCoop** using 

```r
install.packages("pak")
pak::pak("angeclau/AFTCoop")
```

## 📚 Citation
If you use AFTCoop, please cite:

1. C. Angelini, D. De Canditiis, I. De Feis, A. Iuliano. *Cooperative AFT models for multi-omics data integration* submitted (2025)

### 🏛 Funding
This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.
