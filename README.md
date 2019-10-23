# mouse-VOR
A collection of scripts used to analyze and parametrize mouse horizontal VOR for sinusoidal and impulse motions

Run script MouseDataProcessing.m to get data from unsegmented eye velocity to cycle averaged and parametrized.
Run script MouseTableMaker.m to collect and organize metrics for each mouse.
Run script MouseStats.m to run statistics on the tables made by MouseTableMaker.m

All other functions are called by these three scripts. Any function that ends with "Summary" (e.g. MouseSineCleanDataSummary.m) can be used to plot the data at an intermediate step in the process (as a segment, after desaccading/cycle selection, or after analysis). 
