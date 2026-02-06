# CCR-in-the-LGM-Model
Credit Counterparty Risk in the LGM Model (Hull–White vs LGM)
Final project for the Financial Engineering course (Politecnico di Milano, A.Y. 2024–2025). 

This project analyzes counterparty credit risk (CCR) for three interest-rate-swap (IRS) portfolios under two one-factor short-rate models: Hull–White (HW) and Linear Gaussian Markov (LGM). Both models are calibrated to the same set of at-the-money diagonal swaptions; HW jointly calibrates mean reversion and a single volatility level, while LGM uses a piecewise-constant volatility term structure with fixed mean reversion. 

Monte Carlo simulations over a 10-year horizon generate exposure profiles—Expected Exposure (EE), Expected Positive Exposure (EPE), Potential Future Exposure (PFE), and peak PFE—under different configurations (with/without netting and collateral) and multiple time grids (quarterly, weekly, daily). Across the tested setups, HW and LGM produce differences below 5%, with a more visible gap only at the first horizon. 


The repository includes:
MATLAB implementation 
Report with methodology, calibration details, and results


# How to run (MATLAB)
Run the main script:
- `runProject.m`
This entry-point script reproduces the calibration and simulation pipeline and generates the exposure profiles (EE, EPE, PFE, peak PFE) under the configurations discussed in the report.
