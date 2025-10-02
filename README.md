# FDC2Qt <img src="man/figures/logo.png" align="right" height="120" />

Tool for generating **continuous streamflow series (Qt)** in ungauged basins using **period-of-record flow duration curves** (POR-**FDC**).

## Description
**FDC2Qt** generates long and hydrologically consistent synthetic (daily and hourly) streamflow series in ungauged or scarcely gauged catchments.<br>It takes as **input** morpho-climatic features of gauged and ungauged catchments and hydrometric information of gauged catchments (daily/hourly streamflows, mean annual flows, hourly stream stages, rating curves).<br>Its **methodology** consist in three key steps. In particular, for a given ungauged site it estimates:
  1.  a regional POR-FDC: product between mean annual flow and dimensionless FDC (*Index-Flow Approach*, Castellarin et al., 2004)<br>
       1.1.  mean annual flow: multi-regression model based on basin morpho-climatic features (*Stepwise Regression Analysis*, Efroymson, 1960) <br>
       1.2.  dimensionless FDC: weighted average of gauged-catchment FDCs, weights inversely       proportional to the hydrological distance with respect to the ungauged site (*Region-of-Influence Approach*, Burn, 1990)<br>
1. Elemento principale
    1. Sottoelemento 1
    2. Sottoelemento 2
        1. Sotto-sottoelemento 1
  2. 

## Installation
```r
# Install from GitHub
devtools::install_github("alanspado98/FDC2Qt")
```

## Usage
```r
# Install from GitHub
devtools::install_github("alanspado98/FDC2Qt")
```

## Aknowledgments
The author would like to thank prof. A. Castellarin (University of Bologna) for conceptualization and methodology, the Emilia-Romagna Regional Environmental Agency (ARPAE) for the dataset of the testing case study, and the Po River Basin Authority (AdBPo) and Emilia-Romagna Region (RER) for the financial support.  

## References
Castellarin, A., Vogel, R. M. and Brath, A., 2004. A stochastic index flow model of flow duration curves. Water Resources Research, 40(3). <br>
Efroymson, M. A., 1960. Multiple Regression Analysis. In Mathematical Methods for Digital Computers, pp.191–203.<br>
Burn, D. H., 1990. Evaluation of Regional Flood Frequency Analysis With a Region of Influence Approach. Water Resources Research, 26(10), pp.2257–2265.<br>

## License:

## Contact:
