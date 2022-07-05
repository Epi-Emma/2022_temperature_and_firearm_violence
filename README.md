# 2022_temperature_and_firearm_violence
Data and code to conduct time-series analysis of daily ambient temperature and firearm violence incidents across 100 US cities.

[ADD PUBLICATION CITATION HERE ONCE AVAILABLE... THEN CHANGE NAME OF REPOSITORY TO MATCH CITATION]

--------------
**The files:** 

*GVA_Top100_Cities_2015_2020.rds* is a data file containing both the GVA firearm incident data as well as the daily maximum ambient temperature data. The data is formatted so that each city contains a row for every day from 2015 through 2020. 

*DLNM_0dayLag_Median_Ref.R* contains the analysis code for the primary analysis investigating the relationship between daily temperature and firearm incidents. 

*DLNM_2dayLag_Median_Ref.R* contains the analysis code for one of the sensitivity analyses using a two day lag. This code includes arguments for specifying a lag in the DLNM model.


NOTE: The DLNM analysis code requires that you download the "Attributable Risk from a DLNM" function code from Antonio Gasparrini (attrdl.R): https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata


----------------
This analysis benefited from the R code examples written by Antonio Gasparrini and made available on his personal website (http://www.ag-myresearch.com/)


