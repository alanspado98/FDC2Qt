# FDC2Qt <img src="man/figures/logo.png" align="right" height="120" />

Tool for generating **continuous streamflow series (Qt)** in ungauged basins using **period-of-record flow duration curves** (POR-**FDC**).

## Description
**FDC2Qt** generates long and hydrologically consistent synthetic (daily and hourly) streamflow series in ungauged or scarcely gauged catchments.<br>It takes as **input** morpho-climatic features of gauged and ungauged catchments and hydrometric information of gauged catchments (daily/hourly streamflows, mean annual flows, hourly stream stages, rating curves).
Its **methodology** consist in three key steps. In particular, for a given ungauged site it estimates:
  1.  a regional POR-FDC: product between mean annual flow and dimensionless FDC (*Index-Flow Approach*, see Castellarin et al., Water Resour. Res., 2004)<br>
     1.1.  mean annual flow: multi-regression model based on basin morpho-climatic features (*Stepwise Regression Analysis*, see Efroymson, Math. Methods Digit. Comput., 1960)<br>
      1.2.  dimensionless FDC: weighted average of gauged-catchment FDCs, weights inversely proportional to the hydrological distance with respect to the ungauged site (*Region-of-Influence Approach*, Burn, Water Resour. Res., 1990)<br>
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
  
## License:

## Contact:
