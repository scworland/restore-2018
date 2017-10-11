
## Predictions in ungaged basins workplan FY2018
**Scott Worland, William Asquith, and Rodney Knight**,
**U.S. Geological Survey**,
**October 2017**

### Overview

The primary goal of the RESTORE project is to produce quantative and qualitative measurements of watershed alteration for 1320 watersheds in the RESTORE footprint. Accomplishing this goal is dependent on a suite of basin characterisics for each watershed. These basin characteristics include streamflow statistics calculated from a hydrograph. However, XX% of the watersheds in the RESTORE footprint do not have a streamgage and the streamflow statistics need to be estimated for these basins. There are two primary ways to estimate streamflow statistics in ungaged basins: (1) predict a full hydrograph for a basin and calculate the statistics using the estimated hydrograph, or (2) directly predict the streamflow statistics. We are using the second approach.

### Basic idea behind direct prediction

In the simplest terms, direct prediction involves predicting streamflow statistics using basin characteristics derived primarily from gridded products. Examples of streamflow statistics that can be estimated are mean-annual flow, low streamflow values (i.e., 7Q10, minimum October flow, etc), high streamflow values (i.e., 100 year flood discharge), or any number of other values along a flow duration curve. These predictions are made using some form of multivariate regression. This could range from simple linear regression models such as ordinary least squares or censored regression to non-linear machine learning methods such as random forest models or neural networks. 

### Data needed for direct prediction

All of the data needed for direct prediction can be in the form of one large matrix where the required columns are:

* Unique ID for watershed
* Coordinates (centroid of watershed)
* Response variables: streamflow statistics to calculate (for the basins without a gage the values for this column will be NA)
* Explanatory variables: basin characteristics for every basin (including for basins without the streamflow statistics--we will need them for the final predictions)

The matrix will have the following form:

| ID  | 7Q10 | 7Q2  | MAF | ... | X1  | X2  | X3  | ... |
|-----|------|------|-----|-----|-----|-----|-----|-----|
| 001 | 10.2 | 14.8 | xxx | ... | xxx | xxx | xxx | ... |
| 002 | 4.1  | 11.5 | xxx | ... | xxx | xxx | xxx | ... |
| 003 | 21.7 | 40.0 | xxx | ... | xxx | xxx | xxx | ... |
| 004 | NA   | NA   | NA  | ... | xxx | xxx | xxx | ... |
| ... | ...  | ...  | ... | ... | ... | ... | ... | ... |

### Workplan for direct prediction

Scott Worland and William Asquith will start with the same dataset and independently build prediction models using gaged basins for each response variable. Predictions in ungaged basins will be simulated using leave-one-out cross validation (LOO-CV). The LOO-CV predictions for each model and response variable will be stored in a matrix and various perfomance metrics will be calculated for each model. For example, the 7Q10 LOO-CV prediction matrix will have the following form:

| staid | observed 7Q10 | model1 CV-preds | model2 CV-preds | ... |
|-------|---------------|-----------------|-----------------|-----|
| 001   | 10.2          | 8.8             | 12.6            | ... |
| 002   | 4.1           | 0.7             | 3.3             | ... |
| 003   | 21.7          | 18.0            | 24.5            | ... |
| ...   | ...           | ...             | ...             | ... |

There will be separate matrices for each of the response variables. Several of the performance metrics we might consider are:

* Root mean squared error (RMSE)
* Decomposed RMSE (bias, variance, and covariance)
* Unit-area RMSE
* Nash-Sutcliffe coefficient 
* Coefficient of determination
* Mean percent error
* Custom combination of the above
* etc

### Timeline

Below is a basic outline for FY2018. The specifics of the project will by managed using githubs "project" interface.

| Task             | Oct-Nov | Dec-Jan | Feb-March | April-May | June-July | Aug-Sept |
|------------------|---------|---------|-----------|-----------|-----------|----------|
| Gather data      | xxxxxxx |         |           |           |           |          |
| Build models     |         | xxxxxxx | xxxxxxx   | xxxxxxx   |           |          |
| Analyze models   |         |         |           |           | xxxxxxx   |          |
| Make predictions |         |         |           |           | xxxxxxx   | xxxxxxx  |





