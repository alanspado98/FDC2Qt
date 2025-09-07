## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)


#"data" folder with created datasets

usethis::use_data(hydrometric_levels)
#compress="bzip2",“gzip”,“xz” 
usethis::use_data(rating_curves)
usethis::use_data(streamflows)
usethis::use_data(basin_descriptors)