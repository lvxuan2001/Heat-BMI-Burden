# Heat-BMI-Burden

# README
This repository contains part of the code used to reproduce the core analyses in the manuscript:
"Extreme Heat Amplifies Unequal Fiscal Divergence of Health Insurance Burden from Economic Growth in China"

## 1. Overview
This study aims to evaluate the impact of heatwaves on hospitalization burdens for climate-sensitive diseases and its association with economic growth across Chinese cities from 2020 to 2023.

The core analytical pipeline consists of two sequential phases:
Phase 1 — Assessment of Heat-Related Health Risks: Quantifies the exposure-response relationship between extreme heat and daily hospital admissions using a two-stage time-series design with Distributed Lag Non-Linear Models (DLNMs). This phase calculates the relative risk (RR), attributable fraction (AF), and baseline attributable burden.
Phase 2 — Health Insurance Burden Assessment: Translates the health risks into economic burden by integrating insurance reimbursement policies, macroeconomic data (GDP), and population metrics. It calculates the absolute and relative insurance burdens and analyzes their decoupling from economic growth.

## 2. Data Availability and Sources
The core health and insurance data, including national and city-level hospitalization records and basic medical insurance claim data, were sourced from the National Health Commission of China. This data is strictly confidential and is not included in this repository. To ensure research transparency in compliance with data regulations, detailed methodological protocols and de-identified statistical summaries are available from the corresponding author upon reasonable request.

Meteorological data used to calculate the Heat Index were obtained from the publicly accessible ECMWF ERA5 dataset. (https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators)
Air pollution data, included in the models as confounders, were sourced from the publicly available ChinaHighAirPollutants (CHAP) dataset. (https://weijing-rs.github.io/product.html)
Population data were obtained from the LandScan Global Population Dataset, a licensed resource used for calculating per-capita metrics and weighting analyses. (https://landscan.ornl.gov)
Disease-specific national hospitalization rates and per-capita inpatient expenditures, used for deriving national baseline estimates, were compiled from the publicly available China Health Statistics Yearbook. (https://www.nhc.gov.cn/mohwsbwstjxxzx/tjtjnj/tjsj_list.shtml)
The city-tier classification used in this study is provided in Supplementary Table 7 of the manuscript and is also included in this repository (city_tier.docx) for reference.

## 3. Repository Structure
├── README.md
├── data/                   # Data directory
│   └── city_tier.docx      # City classification table (Supplementary Table 7)
└── scripts/                # Main analysis R scripts
    ├── 00_attrdl.R         # Analysis function
    ├── 01_Phase 1_Assessment of heat-related health risks.R
    └── 02_Phase 2_Health insurance burden assessment.R
    
## 4. Contact
For inquiries regarding the methodological details or to request access to de-identified aggregate results under a data agreement, please contact the corresponding author:
Name: Zhenyu Zhang.
Email: zzy@pku.edu.cn.
