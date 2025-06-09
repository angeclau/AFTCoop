---
title: "AFTCoop"
output: rmarkdown::github_document
---

# AFTCoop

**AFTCoop** is an R package for fitting cooperative AFT survival regression models with two omics views.

## Description

The **AFTCoop** package estimates cooperative AFT survival regression models that allow the integration of two matrices of covariates, **U** and **Z**, to improve the estimate. 
It implements the direct cooperative algorithm proposed in Angelini, De Canditiis, De Feis and Iuliano (in preparation 2025). 
It relies on adding a cooperative term into a penalized negative log-likelihood function to enforce the agreement between the two views representing two multi-omics datasets associated to the same patients, such as gene expression and methylation. 
The cooperative learning is a novel approach to multi-view data integration iniatially proposed in Ding et al (2022). 
In this setting, the AFT negative loglikeloood and the LASSO penalization are combined to enforce sparsity, and an agreement penalty is used to encourage the predictions from the two data views to agree while improving the survival estimates. 
The user can choose between three models: *Weibull*, *lognormal*, or *log-logistic*, which handle different noise distributions in the log-linear AFT representation. 
The regression coefficient estimates are obtained using the proximal gradient method, and the regularization parameter is chosen using cross-validation

## 🧪 Installation

For now, install the latest version of **AFTCoop** directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("angeclau/AFTCoop")
```

## 📚 Citation
If you use AFTCoop, please cite:
1. C. Angelini, D. De Canditiis, I. De Feis, A. Iuliano. *Cooperative AFT models for multi-omics data integration* in preparation (2025)

### 🏛 Funding
This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.
