
All files are zipped in 36006963.zip. Data files are in the data folder, figures in the figures folder, and all dependencies in the functions folder.

There are three files to run.

sp500_LP.m will generate replication results for Table 1 and Figures 1, 2, 3, 4, 5 and 6 using the SP500 price index for stock returns data. IRF and historical decomposition for real oil price are given by Figures 1 and 2. IRF for real stock returns are given by Figures 3 and 4 (Figure 4 for the restricted bivariate VAR model), Local Projection IRF for real oil price is given by Figure 5, and Local Projection IRF for real stock returns is given by Figure 6. 

sp500tr.m will generate IRF for real stock returns in Figure 7 using the SP500 total return index for the stock returns data instead of the SP500 price index.

sp500_extended.m will generate results for Table 2 and Figures 8, 9 and 10 using the SP500 index for stock returns data but now with the sampling period from December 2006 to March 2025.

All three files will also generate a figure for the forecast error variance decomposition and report the results corresponding to Table 1 and Table 2 in the command window for the given conditions of each file.