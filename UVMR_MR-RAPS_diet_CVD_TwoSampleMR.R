### Rscript for univariable MR analysis using MR-RAPS.
### Ref panel used for clumping: GRCh37, 1KG Phase3, 2013
### Dietary traits/preferences as exposures.
### T2D and related cardiometabolic traits as outcomes.

##Load packages for MR

library(devtools)
library(tidyverse)
library(TwoSampleMR)


##The genetic instruments (IVs) for exposure/dietary traits were selected via PLINK clumping procedure.

#Load the single summative MR info table consisting of exposure, outcome and potential mediator traits. 
tb = read.table("/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_RAPS_diet_CVD/single_cvd_diet_MR_info_table_both_weak_strong_instruments_for_mr_raps_2.txt",sep="\t",header=T)
#This table contains no duplicate or palindromic variants for exposure traits.
#This table is made for MR-RAPS to test weak instruments of exposure traits of interest to determine whether causal association can be observed in presence of weak instrument bias. 

#Load the table of prevalences for CVDs
prev_table = read.table("/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_diet_CVD/CVD_prevalence.txt",
                        sep='\t',
                        header=T)

#Specify output dir
output_dir = "/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_RAPS_diet_CVD/"

#Specify exposure traits.
exposure = c('STR')

#Specify outcome (T2D) traits
outcome = c('AF)


#Incase using MR-RAPS as a method of choice for analysis causes an error (bug), use the following code for modified version of function mr().
mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj) 
{
  library(TwoSampleMR)
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
method_list_modified <- stringr::str_replace_all(method_list, "mr_raps","mr_raps_modified")
  
mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),function(x1)
  {
    x <- subset(x1, mr_keep)
    
    if (nrow(x) == 0) {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list_modified, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    }
    )
    
    methl <- mr_method_list()
    mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1], 
                         method = methl$name[match(method_list, methl$obj)], 
                         nsnp = sapply(res, function(x) x$nsnp), 
                         b = sapply(res, function(x) x$b), 
                         se = sapply(res, function(x) x$se), 
                         pval = sapply(res, function(x) x$pval))
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    
    return(mr_tab)
  }
  )
  return(mr_tab)
}


#Run this loop for UVMR analysis.
for (e in exposure){
  for(t in outcome){
    #Specify and select exposure and outcome traits of interest.
    #Choose to include both strong and weak instruments or only weak instruments.
    exp_data = tb %>% filter(Type=='Exposure' & Trait==e)
    print(e)
    out_data = tb %>% filter(Type=='Outcome' & Trait==t)
    print(t)
    exp <- format_data(exp_data,
                       type = "exposure",
                       log_pval = FALSE,
                       snps = NULL, 
                       header=TRUE, 
                       phenotype_col = "Phenotype", 
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
    out <- format_data(out_data,
                       type = "outcome",
                       log_pval = FALSE,
                       snps = NULL,
                       header=TRUE, 
                       chr_col = "CHR",
                       pos_col = "POS", 
                       snp_col = "SNP", 
                       beta_col = "BETA", 
                       se_col= "SE", 
                       effect_allele_col = "EA", 
                       other_allele_col= "NEA", 
                       eaf_col = "EAF",
                       pval_col = "P", 
                       samplesize_col = "N", 
                       ncase_col = "CASES", 
                       ncontrol_col = "CONTROLS")
    #Generate harmonized data for exposure-outcome trait pair.
    #action = 2: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative)
    dat = harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
    #Perform MR Steiger test for directionality
    #Calculate variance explained in exposure (using pvals and sample size)
    dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
    
    #Calculate variance explained in binary outcome(using log OR and allele frequencies)
    prev_table2 = prev_table %>% filter(outcome == t) %>% select(prevalence)
    prev = dplyr::pull(prev_table2, prevalence)
    dat$r.outcome <- get_r_from_lor(dat$beta.outcome,
                                    dat$eaf.outcome,
                                    dat$ncase.outcome,
                                    dat$ncontrol.outcome,
                                    prevalence = prev)
    st = mr_steiger2(dat$r.exposure, dat$r.outcome, dat$samplesize.exposure, dat$samplesize.outcome)
    steiger_test = data.frame(exposure = e,
                              outcome = t,
                              snp_r2.exposure = st$r2_exp,
                              snp_r2.outcome = st$r2_out,
                              correct_causal_direction = st$correct_causal_direction,
                              steiger_pval = st$steiger_test)
    write.table(steiger_test, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_steiger_test.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Save the list of SNPs used for the analysis.
    dat2 = dat %>% mutate(Exposure = e, Outcome = t) 
    write.table(dat2,file = paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_SNP_data_table.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Run MR analysis using MR-RAPS as a method. To avoid error from running  mr() with a method of choice as MR-RAPS, use mr_modified().
    res <- mr_modified(dat, method_list = c("mr_raps"))
    res2 = res %>% select(-exposure,-outcome)
    res2$exposure = e
    res2$outcome = t
    write.table(res2, file = paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    #Heterogeneity test.
    het = mr_heterogeneity(dat)
    het2 = het %>% select(-exposure,-outcome)
    het2$exposure = e
    het2$outcome = t
    write.table(het2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_heterogeneity.txt"), sep="\t", quote=F, row.names=F)
    #Pleiotropy test.
    plt = mr_pleiotropy_test(dat)
    plt2 = plt %>% select(-exposure,-outcome)
    plt2$exposure = e
    plt2$outcome = t
    write.table(plt2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_pleiotropy.txt"),sep="\t", quote=F, row.names=F)
    #MR analysis on each SNP individually.
    sin = mr_singlesnp(dat)
    sin2 = sin %>% select(-exposure,-outcome)
    sin2$exposure = e
    sin2$outcome = t
    write.table(sin2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_singleSNP_analysis.txt"),sep="\t", quote=F, row.names=F)
    #Add odds ratios (OR) to MR results. 
    or = generate_odds_ratios(res)
    or2 = or %>% select(-exposure,-outcome)
    or2$exposure = e
    or2$outcome = t
    write.table(or2, file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_OR_with_CI95.txt"),sep="\t", quote=F, row.names=F)
    #Save scatter and forest plots of MR results.
    p1 <- mr_scatter_plot(res2, dat)
    length(p1)
    ggsave(p1[[1]], file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_scatter_plot.pdf"), width=7, height=7)
    #res_single <- mr_singlesnp(dat, all_method=c("mr_raps"))
    #p2 <- mr_forest_plot(res_single)
    #ggsave(p2[[1]], file=paste0(output_dir,"Univariable_MR-RAPS_",e,"_",t,"_forest_plot.pdf"), width=7, height=7)
    options(warn=-1)
  }
}


#Merge MR output files into single files.
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_results.txt   > Univariable_MR-RAPS_dietary_traits_results.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_OR_with_CI95.txt   > Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_heterogeneity.txt   > Univariable_MR-RAPS_dietary_traits_heterogeneity.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_pleiotropy.txt   > Univariable_MR-RAPS_dietary_traits_pleiotropy.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_steiger_test.txt   > Univariable_MR-RAPS_dietary_traits_steiger_test.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_SNP_data_table.txt   > Univariable_MR-RAPS_dietary_traits_SNP_data_table.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR-RAPS_*_*_singleSNP_analysis.txt   > Univariable_MR-RAPS_dietary_traits_singleSNP_analysis.txt"))

#Filter out the pairs that pass Bonferroni-adjusted pval.
#d = read.table(file=paste0(output_dir,"Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"),sep='\t',header=T)

#Filter out exposure-outcome pairs based on the  pass the given Bonferroni-adjusted pval (P < 4.17e-3)
#d2 = d %>% select(-id.exposure,-id.outcome) %>% relocate(exposure,.before=method) %>% relocate(outcome,.after=exposure) %>% filter(pval < 4.17e-3) 

#Save significant associations from MR-RAPS.
#write.table(d2, file=paste0(output_dir,"sig_Univariable_MR-RAPS_dietary_traits_OR_with_CI95.txt"),sep='\t',row.names=F,quote=F)
