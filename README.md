# Analysis of non-linear station motions in terrestrial reference frame computations

Masters Thesis: David Wallinger, DGFI. davidw01123@gmail.com
2019/2020

## About this repository

This repository contains the MATLAB code for the master's thesis 
"Approximation of Non-Linear Post-Seismic Station Motions in the Context of Geodetic Reference Frames".

# Instructions

A generic workflow consists of the following two steps:

1. Pre-process data by converting the raw XYZ files to structured mat-files which can be read by main.
In this way, a coordinate conversion "XYZ to ENU" is also carried out ("convertData_XYZ2ENU.m").
1. Compute fitted values and optimized parameters ("main.m"). Choose input files, settings, and parameters accordingly.

For the stations processed in this work, the first step has been carried out already with the resulting mat-files
stored in the corresponding data directory. This means that the main script can be executed right away.

## User Input for main

Input data

* *inputFolder* - Data input folder. Set to the folder where the mat-files from step 1 are saved.
* *jumpCSVLocation* - Location of Heaviside Jump Database. The CSV file needs to be appropriately formatted.
* *itrfChangesTextfile* - Location of ITRF release dates. The TXT file needs to be appropriately formatted (YYYY-MM-DD)

Result storage

* *doSaveResults* - PNGs of the time series/fitted values will be stored automatically. default false
* *doStaticLogFile* - Log files will be named the same. This means that if the script is executed again for the same station, the old log file will be overwritten. default true

Station to be processed

* *stationName* - ID of the station to be processed. For the master's thesis data, this is the DOMES number.

Extended Trajectory Model (ETM) input parameters

* *polynDeg* - Vector for the degree of the Polynomial Trend Model per coordinate. Recommended [1 1 1]
* *osc* - Cell with vectors for the periods in years of the Seasonal Oscillation Model. Recommended {[0.5 1], [0.5 1], [0.5 1]} 
(annual and semi-annual signal)
* *doITRFjump* - If Heaviside Jumps should be modeled for new ITRF releases. 
Some station might show a small jump, so by setting the corresponding entry in the vector to true, this can be compensated for. Set per coordinate. Recommended [false false false]
* *doEQjump* - If earthquakes in the time series do not cause a noticeable jump, this can be set to true. Set per coordinate. Recommended [false false false]
* *doRemoveTs* - If transients with a lower bound relaxation time later than either the time of the last observation in the trajectory OR than the time of the subsequent earthquake should be removed from the ETM, this can be set to true. Recommended false.
* *doTsOverlay* - True if the one sided limit for transient relaxation times is to be used, False if the two-sided limit is to be used. Recommended false

Combined Fit

* *estimationOpt* - Decides on how relaxation times will be optimized by the NLO. 1 = per coordinate, 2 = for the horizontal coordinates combined, 3 = for all coordinates combined (only makes sense for XYZ data). Recommended 1

Other

* *tarFct* - String of the target function to be optimized. Only the (overall) RMS was found to provide stable results. Recommended 'rms'

Parameters for the Approximation of PSD

* *transientType* - This 3x2 dimensional cell defines the type of approximation per coordinate. The allowed values are either 'log', 'exp' or 'nil'. 
The first value corresponds to the short-term transient displacement, the second value to the long-term one.
Recommended {'log','exp';'log','exp';'log','exp'}
* *tauVec1* - Grid Search parameters for the short-term transient displacement in years. Lower bound:Sampling interval:Upper bound. Recommended years(days(1:10:50))
* *tauVec2* - Grid Search parameters for the long-term transient displacement in years. Lower bound:Sampling interval:Upper bound. Recommended years(days(51:25:365*2))
* *lowLimit* - NLO feasible set parameters. Lower bounds of the short-term and the long-term transient displacement in years, respectively. Recommended [ 0.01/365.25, 51/365.25]
* *uppLimit* - NLO feasible set parameters. Upper bounds of the short-term and the long-term transient displacement in years, respectively. Recommended [50/365.25, 5000/365.25]

The remaining variables should not have to be changed (or require more in depth knowledge of the script - documentation will follow)

Please report suggested improvements or bugs.