# Robust Control SC42145 – Practical Assignment  
## Control Design for a Floating Wind Turbine

This repository contains the MATLAB implementation for the **Practical Assignment** of the course **Robust Control (SC42145)** at **Delft University of Technology**.  
The project focuses on designing and evaluating advanced control strategies for a **floating offshore wind turbine**, addressing challenges such as platform motion, environmental disturbances, and sensitivity to modelling uncertainties.

The repository includes:
- State-space modelling of the floating wind turbine system  
- PID, Mixed Sensitivity H∞ and Fixed Structure controller design  
- Open-loop and closed-loop validation  
- Time-domain simulations  
- Controller analysis using robustness metrics (gain/phase margins, singular values)  
- MATLAB Live Scripts (`.mlx`) and `.m` files for reproducibility  

---

## Authors
- Friso Vijverberg
- Freya van Apeldoorn
- Douwe Brogtrop 

---

## Part 2 — Robust Analysis and Controller Design

This repository also contains the code and Live Script for Part 2 of the assignment, implemented in `Robust_Analysis_Part2.m` (and a corresponding `.mlx` Live Script). Key features:
- Constructs a generalized plant with input/output uncertainty and weighting functions for mixed-sensitivity synthesis.
- Synthesizes a nominal H-infinity controller using `hinfsyn` and performs robust controller design using D-K iteration with `musyn`.
- Performs Fixed-Structure synthesis via `hinfstruct` for comparison and extracts controller gains.
- Computes robustness metrics (NS, NP, RS, RP) plots mu-bounds.
- Includes time-domain utilities `stepWind` and `sineWind` for nominal and uncertain-sample
Requirements: MATLAB with Control System Toolbox and Robust Control Toolbox (for `hinfsyn`, `musyn`, `hinfstruct`, uncertainty tools).

