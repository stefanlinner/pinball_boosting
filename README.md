# Pinball boosting of hedonic regression quantiles

This repository contains material related to the paper "Pinball boosting of hedonic regression quantiles" by Ida Bauer, Harry Haupt, and Stefan Linner. The contents of this repository consist of three main parts: The digital appendix to our paper, the codefile and data to reproduce our application results, and finally the codefiles to reproduce our Monte Carlo simulation study.

## Digital Appendix

The file `Boost_RQ_Appendix.pdf` provides a detailed step-by-step interpretation and comparison of the AL1BRQ and L2BRQ algorithms. It also discusses further details and results of our Monte Carlo study.

## Reproducing application results

The `ad-data.fil` file contains the Canadian housing market data originally used in Anglin and Gencay and also used by our application. Following the `application_hedonic_regression.R` file, you can reproduce the results discussed in "Pinball boosting of hedonic regression quantiles".

## Reproducing the Monte Carlo study

To reproduce our Monte Carlo study, we need the following files  `simulation_setups.R`, `simulation.R`, and `simulation_results.R`. The `simulation_setups.R` file creates a list of different simulation setups that contain the necessary information for the Monte Carlo study, such as the data generation process. The `simulation.R` file runs the Monte Carlo study and saves the results. Running `simulation.R` can take several days (depending on your computer). Therefore, it is generally recommended to run it on a server rather than on a personal computer. Finally, the `simulation_results.R` file can be used to analyze the saved simulation results and to reproduce the tables presented in the main paper as well as in the digital appendix.
