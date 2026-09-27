# Bayesian SEM with Power-Scaling Prior Sensitivity Analysis

This repository contains the replication code and scripts for the study:
> **"Robust Bayesian Structural Equation Modeling with Power-Scaling Prior Sensitivity Analysis for Identifying Determinants of Child Stunting in Java, Indonesia"**

## Abstract
This study contributes a practical framework for assessing parameter robustness in Bayesian Structural Equation Modeling (BSEM) using power-scaling prior sensitivity analysis ($\tau = 2.0, 1, 0.75, 0.50, 0.25$), validated through a Monte Carlo simulation study. Applied to 119 districts/cities in Java, Indonesia, using 2022 data, four latent constructs are modeled: Family Social Protection, Nutrition Consumption Promotion, Specific Nutrition Intervention, and Child Stunting.

## Data Availability
Data are available from the corresponding author upon reasonable request.

## Requirements
- R (version 4.2 or higher)
- R packages: `blavaan`, `rstan`, `ggplot2`

## Repository Structure
- `simulation_study.R`: Code for the Monte Carlo simulation study (800 model estimations).
- `empirical_analysis.R`: Main script for the 119 districts empirical model and power-scaling sensitivity scenarios ($\tau = 2.0, 1, 0.75, 0.50, 0.25$).
- `LICENSE`: MIT License governing the use of these scripts.
