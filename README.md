# SEIR-Austin

R code for running SEIR model for the city of Austin.

The model is an SEIR-style compartmental model with 9 states by 5 age groups and 2 risk groups for a total of 90 states.

## Input
Two input files are used:
1. Hospital data with admission, recoveries and deaths.
2. Mobility data with 7 variables
3. ICU data file for dashboard (not used in main model)

## Run
Main run file is run-seir.R

Run occurs in 3 stages:
1. Fit model parameters using POMP mf2 algorithm.
2. Smoothing stage: select possible paths for states of the model with probability proportional to calculated likelihood based on full time series.
3. For each of the paths attach a future path.

## Commannd line arguments

First argment can be a modifier, or missing.

### Regular run: 
`Rscript run-seir.R Admit_discharge_datafile.xlsx ICY_datafile.xlsx`

### Continuing previous run:
`Rscript run-seir.R continue Rdata_session_file.Rda`
stage of run to continue is determined by name of file before/after stage.

`Rscript run-seir.R  continue Rdata_session_file.Rda after:stage`
continue run after stage. Possible stages: prepare, mf2, smooth, sim, analysis.

### History runs for performance analysis
`Rscript run-seir.R history Admit_discharge_datafile.xlsx [ndays]`
Load Admit_discharge_datafile.xlsx data file, but cut to first [ndays].

### Additional arguments
Additional aguments come in the form of keyword:argument

`no:dash` don't produce dashboard csv files
`no:pdf` don't produce pdf report
`no:mob` run without mobility data

`name.add:string` add string to output file name

