# Modeling Scottish Dialects

This repo contains the code for applying spatial models to questionare data of Scottish dialects. 


## Overview

The project incldes:

- a spatial model for each question.
- a joint spatial model that includes multiple quesitons.

The code produces:

- estimation of the effects of age, gender, and education status on probability of using dialect
- estimation of the spatial distribution of dialect usage frequency together with its uncertainty
- posterior mean and std of the spatial fields
- hierarchical clustering of questions based on similarity of their underlying spatial fields

## Repository structure

```text
scottish_dialect/
├── code/
│   ├── functions.R
│   ├── joint_toy_example.R
│   ├── joint_question_modeling.R
│   ├── separate_question_modeling.R
├── data/
│   ├── csv/
│   │   ├── grammar-xxx.csv
│   │   ├── lexical-xxx.csv
│   │   ├── sounds-about-right-xxx.csv
├── output/
│   ├── figures/
│   ├── tables/
└── README.md
```

## Models
- #### Model A: question-wise simple model that accounts for inhomogeneous sampling efforts and spatial autocorrelation.

Spatial point process for the samples:

> Latent field: $\eta_1(s) = \beta_1 + W_1 (s)$

> Link function: $log(\lambda(s))= \eta_1 (s)$

Spatial model for binarized answers:

> Latent field: $\eta_2(s_i) = \beta_2 + \alpha \cdot W_1(s_i) + W_2(s_i)$

> Link function: $logit(p(s_i))= \eta_2 (s_i)$

- #### Model B: add in factors of age, gender, and education in model A.

> Latent field of the spatial model for answers: 
> 
> $\eta_2(s_i) = \beta_2 + \beta_{age} \cdot {Age}_i + \gamma_{\text{gender}(i)} + \gamma_{\text{uni}(i)} + \alpha \cdot W_1(s_i) + W_2(s_i)$

- #### Model C: joint model for multiple (K) questions that (in addition to model B) accounts for random effect on individuals, with shared age, gender, and education effects.

Spatial point process for the samples:

> Shared latent field: $\eta_{pop}(s) = \beta_{pop} + W_0(s)$

> Link function: $log(\lambda(s))= \eta_{pop} (s)$

> Multiple realizations for likelihood: $S_k \sim \text{PoissonProcess}(\lambda(s)) \quad \text{for } k = 1, \dots, K$

Spatial model for binarized answers:

> Shared latent field: $\eta_{k}(s_{i,k}) = \beta_k + \beta_{\text{age}} \cdot \text{Age}_{i,k} + \gamma_{\text{gender}(i,k)} + \gamma_{\text{uni}(i,k)} + \upsilon_{\text{pid}(i)} + \alpha_k \cdot W_0(s_{i,k}) + W_k(s_{i,k})$

> Link function: $\text{logit}(p_{i,k}) = \eta_k(s_{i,k})$

- #### Model D: similar to model C, but with question-independent age, gender, and education effects.

> Shared latent field: of the spatial model for answers: $\eta_{k}(s_{i,k}) = \beta_k + \beta_{\text{age}, k} \cdot \text{Age}_{i,k} + \gamma_{\text{gender}, k(i)} + \gamma_{\text{uni}, k(i)} + \upsilon_{\text{pid}(i)} + \alpha_k \cdot W_0(s_{i,k}) + W_k(s_{i,k})$


## How to reproduce the analysis




## Citation


