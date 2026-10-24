# Replication Code for Power-Scaling Framework Paper

## Description
This repository contains the complete R code and outputs to replicate the analyses in the paper:
*"A Systematic Power-Scaling Framework for Prior Sensitivity in Small-Sample Bayesian SEM: Evidence from Stunting Correlates in Central Java"*

## Repository Structure
* `analysis rev6.R` - Main R script for the empirical baseline model and power-scaling scenarios.
* `simulation_v6c.R` - R script for the main simulation study (600 estimations).
* `simulation_v7c_supplementer.R` - R script for the supplementary simulation study (300 estimations).
* Figure 1 to Figure 10, and Figure S1 to Figure S3 - from the manuscript.
* `LICENSE` - MIT License.

## How to Run
1. Open the R scripts in RStudio (e.g., `analysis rev6.R`).
2. Install required packages: `blavaan`, `rstan`, `MASS`, `dplyr`, `tidyr`, `ggplot2`, `future`, `future.apply`, `readxl`, etc.
3. Run the scripts line by line. Run `01` for the empirical analysis, and `02`/`03` for the simulations.

## Requirements
* R version 4.5 or higher.
* Stan version 2.32.2 or higher.

## Data and Code Availability
The data that support the findings of this study are available from the official Siperindu database (BKKBN), accessible at (https://siperindu.online/dashboard/#indikator). The replication code and figures are openly available in this repository (https://github.com/TRGteam/power-scaling-bsem) under an MIT license.

## Contact
Open an issue on this repository for questions.
