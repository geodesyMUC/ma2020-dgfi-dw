This ReadMe was written on the 2.9.2019 by David Wallinger @DGFI-TUM

Schematic Workflow of the TSA Trend Computation -----------------------------------------

Input is a file containing a yyyy-mm-dd Timestamp, Station Name, Phi, Lambda and h

1) Run "TSApart1_TransfStationsBLh.m"
	This script converts the Input Coordinates (Latitude, Longitude and h) to
	E(ast), N(orth) and U(p) Coordinates, relative to the FIRST observation of a
	station time series.
	
	Output comprises 1 ".mat" file for every station, 1 time series plot for every
	station.

2A) Run "TSApart2A_TrendComputation.mat"
	This script uses the output from 1) and 
	calculates the trend, depending on the given parameters, for one single
	station. The station needs to be specified in the script code.
	
2B) Run "TSApart2B_allStationsAnalysis.mat"
	This script calculates trends for all ".mat" station coordinate files created
	from 1). Choose Least Square parameters carefully, as they are used globally for
	ALL stations

Notes ------------------------------------------------------------------------------------

Jump Table:
	The jump table ".ods" can be edited in LibreOffice Calc.
	Use the "Save Copy As..." option to save a ".csv" copy, which is used by scripts
	2A and 2B. Make sure to select the correct delimiter ";".
	
	All scripts can be run without a jump table input. Just set the Jump Table Variable 
	to '' (empty string). No Jumps will be considered this way.

Plate Boundaries:
	Plate Boundary txt Files are used for the map visualization in the "analyzeComputedTrends_Position"
	script. They can be downloaded at:
	http://geoscience.wisc.edu/~chuck/MORVEL/PltBoundaries.html

Least Squares:
	The parameters to be estimated in the trend, which can be adjusted in the main scripts ("TSApart2A" or 
	"TSApart2B"), comprise:
		- Polynomial degree of the long-term trend
		- Oscillations with periods in years (365.25d)
		- Parameter T for the logarithmic transient. Default value = 1 (1 year) gives acceptable
		  results (see Bevis Brown 2014)
		- Logical Variable to compensate for ITRF Changes or not. For some stations, this can have a
		  negative effect when set to "true" (when there are other unconsidered factors impacting the TS).

	The parameters which can be set in the Least Squares Algorithm Script ("computeTrendIRLS.m") comprise:
		- Number of iterations to be carried out in the Iteratively Reweighted Least Squares.
		  A value between 1 and 5 (usually sufficient for convergence) is recommended.
		  If set to 0, only one Least Squares Solution will be computed.
		- Factor for the standard deviation used for outlier detection/removal.
	  	  A value between 3 (sensitive, many outliers) and 5 (insensitive, little outliers) is recommended. 
		- The L_p norm parameter p determines the influence of the IRLS reweighting. A value between 1.5
		  and 3 is recommended. 1.5 sets for a robust estimation and decreases the RMS for most station.
		  2 means IRLS has no effect. Larger than 2 values ensue an unstable solution.
		  https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares

Functions used (stored in the "myfunction" directory)
	- computeTrendIRLS: Trend Computation using (Iteratively Reweighted) Least Squares
	- fTSA_TrendComputation: Similar to "TSApart2A_TrendComputation", but written as a
	  callable function used by "TSApart2B_allStationsAnalysis"
	- importfileJumpCS: Imports Jump Table to Matlab Workspace
	- importTectonic: Imports Tectonic Plate Lat,Lon file for visualization purposes
	- importStationTrend: Imports computed trends & parameters for visualization purposes
	- myupdatefcn: Matlab Figure (plot) Data Cursor Add-on
	- VisualizeTS_ENU2: Visualizes Time Series E-N-U and Jumps
	- VisualizeTS_ENU: Visualizes Time Series E-N-U
	- VisualizeTS_Trend_Outliers_ITRF_ENU: Visualizes Time Series E-N-U, computed trend, jumps, ITRF jumps
	  and detected outliers


