# Bayesian SEM with Power-Scaling Prior Sensitivity Analysis

This repository contains the replication code and scripts for the study:
> **"Prior Sensitivity Analysis in Bayesian SEM for Stunting Correlates: A Systematic Power-Scaling Framework"**

## Abstract
This study contributes a practical framework for assessing parameter robustness in Bayesian Structural Equation Modeling (BSEM) using power-scaling prior sensitivity analysis ($\tau = 2.0, 1, 0.75, 0.50, 0.25$), validated through a Monte Carlo simulation study. Applied to 119 districts/cities in Java, Indonesia, using 2022 data, four latent constructs are modeled: Family Social Protection, Nutrition Consumption Promotion, Specific Nutrition Intervention, and Child Stunting.

## Data Availability
The data analysed in this study are publicly accessible through the Siperindu database (Population Control Early Warning Information System), managed by the Indonesian Ministry of Population and Family Development (Badan Kependudukan dan Keluarga Berencana Nasional [BKKBN]), and can be accessed at https://siperindu.online/dashboard/#indikator. The dataset covers the year 2025..

## Requirements
- R (version 4.2 or higher)
- R packages: `blavaan`, `rstan`, `ggplot2`

## Repository Structure
- `simulation_study.R`: Code for the Monte Carlo simulation study (800 model estimations).
- `empirical_analysis.R`: Main script for the 119 districts empirical model and power-scaling sensitivity scenarios ($\tau = 2.0, 1, 0.75, 0.50, 0.25$).
- `LICENSE`: MIT License governing the use of these scripts.
