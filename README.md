# 2022_temperature_and_firearm_violence
This repository contains code to replicate a time-series analysis of daily ambient temperature and firearm violence incidents across 100 US cities from 2015 - 2020. 

[ADD PUBLICATION CITATION HERE ONCE AVAILABLE... THEN CHANGE NAME OF REPOSITORY TO MATCH CITATION]

--------------
**The files:** 

*GVA_Incidents_100cities_2015_2020.rds* is a data file that contains a count of Gun Violence Archive recorded firearm incidents per city, per day, from 2015-2020 in the 100 analysis cities.  

*DLNM_0dayLag_Median_Ref.R* contains the analysis code for the primary analysis investigating the relationship between daily temperature and firearm incidents with a 0-day lag. 

*DLNM_2dayLag_Median_Ref.R* contains the analysis code for one of the sensitivity analyses using a 2-day lag. This code includes arguments for specifying a lag in the DLNM model.


NOTE: The DLNM analysis code requires that you download the "Attributable Risk from a DLNM" function code from Antonio Gasparrini (attrdl.R): https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata


----------------
This analysis benefited from the R code examples written by Antonio Gasparrini and made available on his [personal website](http://www.ag-myresearch.com/).

The firearm incidents data was downaloded from [Gun Violence Archive](https://www.gunviolencearchive.org/) (GVA). These are the terms of use from the GVA website:
>Gun Violence Archive (GVA) is a not for profit corporation formed in 2013 to provide free online public access to accurate information about gun-related violence in the United States. GVA will collect and check for accuracy, comprehensive information about gun-related violence in the U.S. and then post and disseminate it online, primarily if not exclusively on this website and summary ledgers at www.facebook.com/gunviolencearchive and on Twitter @gundeaths. It is hoped that this information will inform and assist those engaged in discussions and activities concerning gun violence, including analysis of proposed regulations or legislation relating to gun safety usage. All we ask is to please provide proper credit for use of Gun Violence Archive data and advise us of its use.


