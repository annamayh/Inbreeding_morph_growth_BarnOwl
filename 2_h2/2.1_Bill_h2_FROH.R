.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path1-11

library(tidyverse)
library(brms)
library(corpcor)

##################################################################################
########################### ~~ bill ~~ #############################################
#####################################################################################

args = commandArgs(trailingOnly = TRUE)
scratch = as.character(args[1]) 

bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T) ##
  
bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)

bill_df$rank=as.numeric(bill_df$rank)
bill_df$mc_age_acc <- bill_df$age_acc - mean(bill_df$age_acc)

## read in GRM
grm=read_rds("./input_dfs/All3085_AUTOSAUMES_RP502SNPs.RDS")
grm_filt=grm[rownames(grm)%in%bill_df[["RingId"]],colnames(grm)%in%bill_df[["RingId"]]]##filtering grm for ids only in df

grm_filt_pd <- make.positive.definite(grm_filt)
GRM <- as(grm_filt_pd, "dgCMatrix")

## same as priors for simple model + prior for pe ~= 0 because we expect little var 
prior_bill=c(prior(student_t(3, 180,20), class = "Intercept"), ## 
             prior(student_t(3,0,20), class = "sd"),
             prior(normal(0, 5), class = "sd", group="RingId_pe"))

mod_bill_GRM_h2 <- brm(BillLength ~  1 + sex + mc_age_acc + FROH +
                           (1|gr(RingId, cov=Amat)) + (1|RingId_pe)  + (1|clutch_merge) +
                           (1|year) + (1|month) + (1|nestboxID) + (1|rank) +
                           (1|Observer),
                      
                         data = bill_df,
                         control=list(adapt_delta=0.95),
                         data2 = list(Amat = GRM),
                         chains = 4,
                         cores=4,
                         prior=prior_bill, ##
                         iter = 25000,
                         warmup = 10000,
                         thin=5
)


saveRDS(mod_bill_GRM_h2,file=paste0(scratch,"2.1.bill_h2_FROH.RDS")) ##

summary(mod_bill_GRM_h2) ###
