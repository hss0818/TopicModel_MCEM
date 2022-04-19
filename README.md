# TopicModel_MCEM

Source code for the MCMC-EM algorithm for topic model estimation, which is proposed in the paper **Learning Topic Models: Identifiability and Finite-Sample Analysis**.

## Data

https://www1.nyc.gov/site/tlc/about/tlc-trip-record-data.page
`taxi_19.csv` contains processed New York City taxi pickup and dropoff data downloaded from https://s3.amazonaws.com/nyc-tlc/trip+data/yellow_tripdata_2019-01.csv


`MCEM_helper.cpp` contains the helper functions for the implementation of the main function. 

`MCEM.R` contains the main function for the MCMC-EM algorithm proposed in the paper. 

`example.R` contains an example showing how to use the main function to make estimation.

