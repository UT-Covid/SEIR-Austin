# SEIR-Austin

R code for running SEIR model for the city of Austin.

The model is an SEIR-style compartmental model with 9 states by 5 age groups and 2 risk groups for a total of 90 states.

## Input
Two input files are used:
1. Hospital data with admission, recoveries and deaths.
2. Mobility data with 7 variables

## Run
Main run file is run-seir.R

Run occurs in 3 stages:
1. Fit model parameters using POMP mf2 algorithm.
2. Smoothing stage: select possible paths for states of the model with probability proportional to calculated likelihood based on full time series.
3. For each of the paths attach a future path.

## Commannd line arguments
