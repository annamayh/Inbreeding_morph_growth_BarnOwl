.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(corpcor)
library(readr)

##################################################################################
########################### ~~ Mass ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

mass_df=read.table("./input_dfs/mass_all_pheno_df.txt",sep=",", header=T)
mass_df = mass_df[mass_df$gr_stage == "Adult", ]

mass_df$clutch_merge=as.factor(mass_df$clutch_merge)
mass_df$sex=as.factor(mass_df$sex)
mass_df$RingId=as.factor(mass_df$RingId)
mass_df$year=as.factor(mass_df$year)
mass_df$Observer=as.factor(mass_df$Observer)
mass_df$nestboxID=as.factor(mass_df$nestboxID)

mass_df$rank=as.numeric(mass_df$rank)


# mean centering for age so estimates are for an individual of an ave age (not 0 days old)
mass_df$mc_age_acc <- mass_df$age_acc - mean(mass_df$age_acc)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%mass_df[["RingId"]],colnames(grm)%in%mass_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")


mod_mass_GRM.split_stage_h2_adult <- brm(Mass ~  1 +sex+FROH+ # 
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe) + (1|Observer) + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID) + (1|rank), 
                      
                         data = mass_df,
                         control=list(adapt_delta=0.97),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         iter = 35000,
                         warmup = 7500,
                         thin=5
)


saveRDS(mod_mass_GRM.split_stage_h2_adult,file=paste0(scratch,"2.2.mass_split_h2_adult_FROH.RDS")) ##

summary(mod_mass_GRM.split_stage_h2_adult) ###
