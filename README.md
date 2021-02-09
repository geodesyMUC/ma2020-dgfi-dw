# Approximation of Non-Linear Post-Seismic Station Motions in the Context of Geodetic Reference Frames

Masters Thesis: David W., [DGFI](https://www.dgfi.tum.de/en/) & [MUAS](https://www.geo.hm.edu/), 2020
## About this repository

This repository contains the MATLAB code, folder structure, and preprocessed files for the Masters Thesis 
"Approximation of Non-Linear Post-Seismic Station Motions in the Context of Geodetic Reference Frames". Results are created by running the "main" script.

# Instructions

A generic workflow consists of the following two steps:

1. "**convertData_XYZ2ENU.m**" - Pre-process data by converting the raw ASCII files containing the XYZ coordinates to structured mat-files which can be read by the main script.
In the process, a coordinate conversion "XYZ to ENU" is also carried out. Both ENU and XYZ files are saved which can then be processed in the second step.
1. "**main.m**" - Compute fitted values and optimized parameters, create plots . Choose input files, directories, settings, and parameters accordingly.

For the stations processed in this work, the first step has been carried out already with the resulting mat-files
stored in the corresponding data directory. This means that the main script can be executed right away.

## User Input for main

Input data

* *inputFolder* - Data input folder. Set to the folder where the mat-files from step 1 are saved. The preprocessing script "convertData_XYZ2ENU.m" creates the mat-files for ENU and XYZ coordinates, so just this path variable has to be changed to either ENU or XYZ.
* *jumpCSVLocation* - Location of Heaviside Jump Database. The CSV file needs to be appropriately formatted. The file for the jumps and earthquakes processed for this work is saved in the src-folder and is loaded by default when running "main". It can be changed and extended, but the formatting must not be changed.
* *itrfChangesTextfile* - Location of ITRF release dates. The TXT file needs to be appropriately formatted (YYYY-MM-DD datetime, newline delimited). The file for the current relevant ITRF releases is located in the src-folder and loaded by default when running "main"..

Result storage

* *doSaveResults* - PNGs of the time series/fitted values and the residual plots will be stored automatically. default false
* *doStaticLogFile* - Log file name will be kept the same for each station. This means that if the script is executed again for the same station, the old log file will be overwritten. default true

Station to be processed

* *stationName* - ID of the station to be processed, corresponding to the name of the mat file in the denoted input directory. For the master's thesis data, this is the DOMES number.

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

* *tarFct* - String of the target function to be optimized. Only the (overall) RMS was found to provide stable results. Recommended 'rms', but can be changed to 'rms35d', 'rms105d', 'rms1y', or 'rms2y'.

Parameters for the Approximation of PSD

* *transientType* - This 3x2 dimensional cell defines the type of approximation per coordinate. The allowed values are either 'log', 'exp' or 'nil'. 
The first value corresponds to the short-term transient displacement, the second value to the long-term one.
Recommended {'log','exp';'log','exp';'log','exp'} or {'log','log';'log','log';'log','log'} (**2 transients** per eq), or {'log','nil';'log','nil';'log','nil'} (1 transient per eq)
* *tauVec1* - Grid Search parameters for the short-term transient displacement in years. Lower bound:Sampling interval:Upper bound. Recommended years(days(1:10:50))
* *tauVec2* - Grid Search parameters for the long-term transient displacement in years. Lower bound:Sampling interval:Upper bound. Recommended years(days(51:25:365*2))
* *lowLimit* - NLO feasible set parameters. Lower bounds of the short-term and the long-term transient displacement in years, respectively. Recommended [ 0.01/365.25, 51/365.25]
* *uppLimit* - NLO feasible set parameters. Upper bounds of the short-term and the long-term transient displacement in years, respectively. Recommended [50/365.25, 5000/365.25]

The remaining variables such as IRLS related variables should **only** be changed with care. They might result in undesired output and require more in depth knowledge of the script.


## Folders
* *data_psd* : Data input for main
	- XYZ: Data input for main, XYZ coordinates (Files already provided for existing data. If new data is to be processed, the preprocessing script needs to be run with the new data as Input)
	- ENU: Data input for main, ENU coordinates (Files already provided for existing data. If new data is to be processed, the preprocessing script needs to be run with the new data as Input)
	- XYZ_plots: Raw data plots (Plots already provided for existing data)
	- ENU_plots: Raw data plots (Plots already provided for existing data)
* *data_sirgasStationsENU* : SIRGAS stations data input for main, might need rework
* *data_template* : Test data input for main
* *myfunctions* : Functions used by main
* *raw_data* : Raw data from the DGFI. A conversion to the main input file format (mat files) has already been carried out, the results are saved in the data_psd directories
* *results_logs* : Folder for log files, will be created if script runs for the first time
* *saved_plots* : Folder for plots, will be created if script runs for the first time
* *src* : Important supplementary input for main such as the jump table (jumps and earthquakes (trigger transients!) in the station trajectories)
* *visualize_results_for_thesis* : Functions that were used to create plots for the Masters Thesis. They might need rework in oder to work properly. Use with care.

## Data Template

It is possible to create own, user-defined input for the main script WITHOUT having to use the preprocessing script. A reason might be that raw data is delivered in a different format than the formatted .x, .y., and .z files.
How this can be done is shown in the "*Create_Data_template*" script. Following the same steps, any time series data can be formatted in a way that it can be read and processed by the main script.

## Feedback

Please report suggested improvements or bugs to the author.
