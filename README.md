# 2022_temperature_and_firearm_violence
This repository contains code from a time-series analysis of daily ambient temperature and firearm violence incidents across 100 US cities from 2015 - 2020. 

Lyons VH, Gause EL, Spangler KR, Wellenius GA, Jay J. Analysis of Daily Ambient Temperature and Firearm Violence in 100 US Cities. *JAMA Netw Open*. 2022;5(12).


--------------
**The files:** 

- *DLNM_0dayLag_Median_Ref.R* contains the analysis code for the primary analysis investigating the relationship between daily temperature and firearm incidents with a 0-day lag. 

- *DLNM_2dayLag_Median_Ref.R* contains the analysis code for one of the sensitivity analyses using a 2-day lag. This code includes arguments for specifying a lag in the DLNM model.


NOTE: The DLNM analysis code requires that you download the "Attributable Risk from a DLNM" function code from Antonio Gasparrini (attrdl.R): https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata


----------------
This analysis benefited from the R code examples written by Antonio Gasparrini and made available on his [personal website](http://www.ag-myresearch.com/).
