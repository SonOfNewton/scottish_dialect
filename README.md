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

## How to reproduce the analysis




## Citation


