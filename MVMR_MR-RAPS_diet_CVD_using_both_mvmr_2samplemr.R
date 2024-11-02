##This script is to run "single-confounder" multivariable MR on MR-RAPS associations via 2SampleMR and MVMR packages.
##It prepares and formats input data first to run via 2SampleMR package first.
##Then, it uses the input data created via 2SampleMR to format the data and run multivariable MR via MVMR package.
##This approach allows harmonize exposure,confounder and outcome data in MVMR and thereby compare multivariable MR results between 2SampleMR and MVMR. 
##The codes are obtained and sourced from https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/
##The functions and steps are also described and provided by https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#multivariable-mr

#Load packages first 
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)


#Load the single summative MR info table consisting of exposure, outcome and mediator traits (both strong and weak instruments)
tb = read.table("/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_RAPS_diet_CVD/MVMR_MR_RAPS/single_cvd_diet_MR_info_table_both_weak_strong_instruments_5mediators_for_mr_raps.txt",sep='\t',header=T)

#Specify exposure trait (ALCMEAL, CHAMPWH, DRIEDFRU, MUESLI)
exposure = c("DRIEDFRU")

#Specify outcome trait (CAD, PAD, HF)
outcome = c("PAD")

#list confounder(s) of interest
#For MR-RAPS associations in Diet-CVD, only BMI and EA  are tested as potential mediators.
mediator = c("BMI")

#Specify output directory
output_dir = "/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_RAPS_diet_CVD/MVMR_MR_RAPS/"

# Creating a tidy outcome for 2SampleMR results
tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=4) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

#Create tidy outcome for MVMR results
tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

#Create a function "make_mvmr_input" for MVMR 

make_mvmr_input <- function(exposure_dat, outcome.data = outcome_dat){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}



for (e in exposure){
  for (t in outcome){
    for (m in mediator){
      #Load and format the data for exposure, confounder and outcome traits 
      exp_dat = tb %>% filter(Type=='Exposure' & Trait==e) 
      med_dat = tb %>% filter(Type=='Mediator' & Trait==m) 
      exp_data = rbind(exp_dat,med_dat)
      exposure_dat <- format_data(exp_data, type = "exposure",
                                  log_pval = FALSE,
                                  snps = NULL,
                                  header=TRUE,
                                  phenotype_col = "Trait",
                                  chr_col = "CHR",
                                  pos_col = "POS",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  effect_allele_col = "EA",
                                  other_allele_col="NEA",
                                  eaf_col="EAF",
                                  pval_col="P",
                                  samplesize_col = "N")
      out_data = tb %>% filter(Type=='Outcome' & Trait==t) 
      outcome_dat <- format_data(out_data, 
                                 type = "outcome",
                                 log_pval = FALSE,
                                 snps = NULL,
                                 header=TRUE,
                                 phenotype_col = "Trait",
                                 chr_col = "CHR",
                                 pos_col = "POS",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col= "SE",
                                 effect_allele_col = "EA",
                                 other_allele_col= "NEA",
                                 eaf_col = "EAF",
                                 pval_col = "P",
                                 samplesize_col = "N")
      
      #Harmonise so that all are on the same reference allele. Also, this step will remove any SNPs that have incompatible alleles or are duplicates.
      mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
      
      #Create and save data table for SNPs being included and harmonized for MR
      capture.output(df = data.frame(exposure.beta=mvdat$exposure_beta,
                                     exposure.pval=mvdat$exposure_pval,
                                     exposure.se=mvdat$exposure_se,
                                     outcome.beta = mvdat$outcome_beta,
                                     outcome.pval=mvdat$outcome_pval,
                                     outcome.se = mvdat$outcome_se),
                     file = paste0(output_dir,"Multivariable_MR_",e,"_",m,"_",t,"_SNP_data_table.txt")) 
      
      #Finally, perform the multivariable MR analysis
      res <- mv_multiple(mvdat, pval_threshold = 1)
      
      # Creating a tidy outcome
      result_2smr <- res$result %>%
        split_outcome() %>%
        separate(outcome, "outcome", sep="[(]") %>% 
        mutate(outcome=stringr::str_trim(outcome))%>% 
        generate_odds_ratios() %>% 
        select(-id.exposure, -id.outcome) %>% 
        tidy_pvals()
      #Print out 2SampleMR multivariable MR results
      result_2smr
      #Save 2SampleMR results
      write.table(result_2smr, file=paste0(output_dir,"Multivariable_MR_",e,"_",m,"_",t,"_2smr_ivw_results.txt"),sep='\t',row.names=F,quote=F)
      
      #Use the make_mvmr_input function to create input data for MVMR 
      mvmr_input <- make_mvmr_input(exposure_dat, outcome.data=outcome_dat)
      
      #Check the input data 
      glimpse(mvmr_input$XGs)
      glimpse(mvmr_input$Y)
      
      # format data to be in MVMR package-compatible df
      mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  # exposure betas
                              BYG = mvmr_input$YG$beta.outcome,                     # outcome beta
                              seBXGs = mvmr_input$XGs %>% select(contains("se")),  # exposure SEs
                              seBYG = mvmr_input$YG$se.outcome,                     # outcome SEs
                              RSID = mvmr_input$XGs$SNP)    
      
      head(mvmr_out)
      
      #Estimate causal effects using method in MVMR package
      mvmr_res <- ivw_mvmr(r_input=mvmr_out)
      
      #Tidy up the output format
      #Specify the outcome trait!
      result_mvmr <-
        mvmr_res %>% 
        tidy_mvmr_output() %>% 
        mutate(exposure = mvmr_input$exposures,
               outcome = t) %>% 
        select(exposure, outcome, everything()) %>% 
        tidy_pvals()
      
      #Print out multivariable MR results via MVMR package
      result_mvmr
      
      #Save MVMR multivariable MR results
      write.table(result_mvmr, file=paste0(output_dir,'Multivariable_MR_',e,'_',c,'_',t,'_mvmr_ivw_results.txt'), sep='\t',row.names=F, quote=F)
    }
  }
}







