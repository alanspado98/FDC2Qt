# FDC2Qt <img src="man/figures/logo.png" align="right" height="130" />

Tool for generating **continuous streamflow series (Qt)** in ungauged basins using **period-of-record flow duration curves** (POR-**FDC**).

## Description
**FDC2Qt** generates hydrologically consistent synthetic streamflow series (daily/hourly) in ungauged or scarcely gauged catchments.<br>It takes as **input** morpho-climatic features of gauged and ungauged catchments and hydrometric information of gauged catchments (daily/hourly streamflows, mean annual flows, hourly stream stages, rating curves).<br>Its **methodology** consist in three key steps. In particular, for a given ungauged site it estimates:
  1.  a regional POR-FDC: product between mean annual flow and dimensionless FDC (*Index-Flow Approach*; Castellarin et al., 2004)
       1. mean annual flow: multi-regression model based on basin morpho-climatic features (*Stepwise Regression Analysis*; Efroymson, 1960)
       2. dimensionless FDC: weighted average of gauged-catchment FDCs, weights inversely proportional to the hydrological distance with respect to the ungauged site (*Region-of-Influence Approach*; Burn, 1990)
  2.  a daily streamflow series: nonlinear interpolation method based on POR-FDCs, after selecting a donor gauged catchment (Smakhtin et al., 1997)
  3.  a hourly streamflow series:
        1. linear interpolation for low flows
        2. second-order polynomials as approximation of rising and falling limbs of flood events, honoring four constraints:
             1. hydrograph continuity
             2. daily flood volumes
             3. hourly flood peak
             4. hydrograph shape (Maione et al., 2003)

<p align="center">
  <img src="man/figures/workflow.png" height="400" />
  <br>
  <em>Figure 1. Workflow of the package.</em>
</p>

Auxiliary functions allow users to perform data quality assessment, by looking at precipitation elasticity of streamflows (Sankarasubramanian et al., 2001) and at the shape of empirical POR-FDCs. <br>     
FDC2Qt was applied to 14 gauged catchments across the Emilia-Romagna region in Northern Italy. LOOCV of hydrological singatures across these sites demonstrates a good accurancy of synthetic daily streamflow series, while a perliminary test against a lumped rainfall-runoff model calibrated using the same regional POR-FDC indicate higher reliability of FDC2Qt hourly streamflow series, especially in flood events.


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
The author want to thank prof. A. Castellarin (University of Bologna) for conceptualization and methodology, the Emilia-Romagna Regional Environmental Agency (ARPAE) for the dataset of the testing case study, and the Po River Basin Authority (AdBPo) and Emilia-Romagna Region (RER) for the financial support.  

## Contact
* Alan Spadoni: alan.spadoni2@unibo.it
* Attilio Castellarin: attilio.castellarin@unibo.it

## References
* Burn, D. H., 1990. *Evaluation of Regional Flood Frequency Analysis With a Region of Influence Approach*. Water Resources Research, 26(10), pp.2257–2265.
* Castellarin, A., Vogel, R. M. and Brath, A., 2004. *A stochastic index flow model of flow duration curves*. Water Resources Research, 40(3).
* Efroymson, M. A., 1960. *Multiple Regression Analysis*. In Mathematical Methods for Digital Computers, pp.191–203.
* Maione, U., Mignosa, P., Tomirotti, M., 2003. *Regional estimation of synthetic design hydrographs*. International Journal of River Basin Management, 1(2), pp.151–163.
* Sankarasubramanian, A., Richard M. V., James F. L., 2001. *Climate Elasticity of Streamflow in the United States*. Water Resources Research, 37 (6), pp.1771–81.
* Smakhtin, V. Y., Hughes, D. A., Creuse-Naudin, E., 1997. *Regionalization of Daily Flow Characteristics in Part of the Eastern Cape*, South Africa. Hydrological Sciences Journal 42(6), pp.919–936.
