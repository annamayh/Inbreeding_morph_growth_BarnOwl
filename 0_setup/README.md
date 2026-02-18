This folder conatains scripts used to create the dataframes directly used in the analysis in folders 1-4, the output can be found on my FigShare (see previous readme). 

#Allrecords_mkdfs.R# combines all phenotypic and ecological records from the long term study collected by Alexandre Roulin and Bettina Almasi with genomic inbreeding coefficients estimated in this study into one dataframe
#SSGompertz_params.R# estimates population-wide paramters values for non-linear growth of traits to account for age within the LMM (1) and h2 (2) models 
#Calc_FHBD.R and Calc_FuniW# detail how the two genomic inbreeding coefficients were calculated
Finally, #FuniW_split_into_chunks.R# and #extract_HBD_segs_cluster.R# calculate local inbreeding coefficients for certain regions for the genetic architecture analysis (4)
