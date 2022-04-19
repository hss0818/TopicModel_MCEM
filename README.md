# TopicModel_MCEM

Source code for the MCMC-EM algorithm for topic model estimation, which is proposed in the paper **Learning Topic Models: Identifiability and Finite-Sample Analysis**.

## Data

`taxi_19.csv` is the discretized New York City taxi-trip dataset. This dataset contains 7,667,792 New York City yellow cab trips in January 2019. The location information is discretized into 263 taxi zones with 69 in Manhattan, 69 in Queens, 61 in Brooklyn, 43 in Bronx, 20 in Staten Island and 1 in EWR. The original file can be accessed from the webpage https://www1.nyc.gov/site/tlc/about/tlc-trip-record-data.page and downloaded via https://s3.amazonaws.com/nyc-tlc/trip+data/yellow_tripdata_2019-01.csv


## Code

`MCEM_helper.cpp` contains the helper functions for the implementation of the main function. 

`MCEM.R` contains the main function for the MCMC-EM algorithm proposed in the paper. The log-likelihood is printed after each iteration.

`MCEM_k.R` contains the main function for the MCMC-EM algorithm. The log-likelihood is not calculated in the process.

`example.R` contains an toy example showing how to use the main function to make estimation.

`run_kos_nips.R` implements the application to the NIPS dataset and the Daily Kos dataset.

`run_taxi.R` implements the application to the New Yorks City taxi-trip dataset.
