
##################################################################################
######### ~~ Inbreeding depression in growth rates of bill length ~~ ################
#####################################################################################

## This script has been annotated in more detail for reproducibility but annotations apply to other NLMM run 
# Most of these R script were run on the HPC at UNIL hence the file paths #
.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(brms)
library(dplyr)

args = commandArgs(trailingOnly = TRUE) # trailing arguments in slurm script to link to file paths within R script
scratch = as.character(args[1]) # output scratch directory

# Reading in df
bill_df=read.table("./input_dfs/bill_all_pheno_df.txt",sep=",", header=T)

# Sorting out variables
bill_df$clutch_merge=as.factor(bill_df$clutch_merge)
bill_df$sex=as.factor(bill_df$sex)
bill_df$RingId=as.factor(bill_df$RingId)
bill_df$year=as.factor(bill_df$year)
bill_df$Observer=as.factor(bill_df$Observer)
bill_df$nestboxID=as.factor(bill_df$nestboxID)
bill_df$rank=as.numeric(bill_df$rank)

# Setting priors
## NOTE: for a NLMM the priors needed to be relatively strong based on recommendations, otherwise the model is unlikely to converge 
prior_bill_gr<- c(
  prior(normal(184, 20), nlpar = "asym",  class="b", coef="Intercept"),## prior for intercept of asymptote (prior informed by population level info)
  prior(normal(1, 1), nlpar = "b",  class="b", coef="Intercept"), ## prior for intercept of intiation of the growth cure 
  prior(normal(log(0.93), 1), nlpar = "c",  class="b", coef="Intercept"), ## prior for intercept of the growth slope (prior informed by population level info)
  
  prior(normal(0, 10), nlpar = "asym",  class="b"),## prior for general fixed effects of asymptote - this value is the deviation from the intercept  
  prior(normal(0, 1), nlpar = "c",  class="b"), ## prior for general fixed effects in the growth rate parameter
  
  prior(student_t(3, 0, 20), class = "sigma"), # residual 
  prior(student_t(3, 0, 20), class="sd",nlpar = "asym"), # general prior for random effects in asymptote (class=sd)
  prior(student_t(3, 0, 1),  class="sd", nlpar = "c"), # general prior for random effects in growth rate parameter (class =sd)
  
  prior(normal(0, 10), class="sd", group="RingId", nlpar = "asym"), #
  prior(normal(0, 0.5),  class="sd", group="RingId", nlpar = "c")
  
)


growth_bill.mod=brm(
  ## model 
  bf(BillLength ~ asym * exp(-b*(exp(c))^age_days), # define general formula 
     ## 
     asym ~ 1 + FROH + rank + sex + (1|clutch_merge) + (1|Observer) + (1|year) + (1|month) + (1|RingId) + (1|nestboxID), # covariartes to evulate their effect on asymptote
     b ~ 1, ## Initiation of growth curve - fixed 
     c ~ 1 + FROH + rank + sex + (1|RingId), ## covariates to evaulate effect on growth rate parameters (c) 
     nl=TRUE),
  data=bill_df, 
  family = gaussian(),
  chains = 4,
  prior = prior_bill_gr,
  control = list(adapt_delta = 0.99),
  init = 0, 
  cores = 4,
  iter = 35000, 
  warmup = 15000, 
  thin=5
  
)


saveRDS(growth_bill.mod,file=paste0(scratch,"3.1.bill_gr_FROH_subset.RDS")) ## saving model

summary(growth_bill.mod) # checking model fit

