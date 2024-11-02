### Re-perform three univariable MR analyses after correction (Ref panel used:GRCh37, 1KG Phase3, 2013)
### T2D and T2D traits as exposures
### Dietary traits as outcomes

##Load packages for MR

library(devtools)
library(tidyverse)
library(TwoSampleMR)


##The genetic instruments (IVs) for exposure/dietary traits are identified via PLINK clumping.

#Load the single summative MR info table consisting of exposure, outcome and confounder traits 
tb = read.table("/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MVMR_diet_CVD/single_table_MR_info_table_phase3_CVD_endpoints_5mediators.txt",sep="\t",header=T)

#Load the table of prevalences for CVDs
prev_table = read.table("/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_diet_CVD/CVD_prevalence.txt",
                  sep='\t',
                  header=T)

#Specify output dir
output_dir = "/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/MR_diet_CVD/"

#Specify exposure traits
exposure = c("ALCMEAL",
             "ALC",
             "BEEF",
             "BUTTER",
             "BUTMARG",
             "CARB",
             "CHAMPWH",
             "CHEESE",
             "COF",
             "COOKEDVEG",
             "CORNFLAK",
             "DRIEDFRU",
             "FAT",
             "FRESHFRU",
             "MUESLI",
             "NONOILYFSH",
             "PORK",
             "POULTRY",
             "PROTEIN",
             "RAWVEG",
             "REDWINE",
             "SKIMMLK",
             "SPREADS",
             "SUGAR",
             "WHITEBRD",
             "WHOLEBRD",
             "WHOLEMLK",
             "ACQ",
             "CAFSWT",
             "COFALC",
             "DESS",
             "FATSALT",
             "LOWCAL",
             "PAL",
             "SAVCAL",
             "SAVOUR",
             "STR",
             "VEG")


#Specify outcome (T2D) traits
outcome = c("VTE",
            "AAA",
            "HF",
            "AF",
            "CAD",
            "PAD",
            "AS")

#Run this loop for UVMR analysis on exposure and outcome traits of interest
for (e in exposure){
  for(t in outcome){
    #Choose to include/exclude palindromic IVs 
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
    
    #Generate harmonized data for exposure-outcome trait pairs
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
    write.table(steiger_test, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_steiger_test.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    
    dat2 = dat %>% mutate(Exposure = e, Outcome = t) 
    write.table(dat2,file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_SNP_data_table.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    res <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    res2 = res %>% select(-exposure,-outcome)
    res2$exposure = e
    res2$outcome = t
    write.table(res2, file = paste0(output_dir,"Univariable_MR_",e,"_",t,"_results.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    het = mr_heterogeneity(dat)
    het2 = het %>% select(-exposure,-outcome)
    het2$exposure = e
    het2$outcome = t
    write.table(het2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_heterogeneity.txt"), sep="\t", quote=F, row.names=F)
    plt = mr_pleiotropy_test(dat)
    plt2 = plt %>% select(-exposure,-outcome)
    plt2$exposure = e
    plt2$outcome = t
    write.table(plt2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_pleiotropy.txt"),sep="\t", quote=F, row.names=F)
    sin = mr_singlesnp(dat)
    sin2 = sin %>% select(-exposure,-outcome)
    sin2$exposure = e
    sin2$outcome = t
    write.table(sin2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_singleSNP_analysis.txt"),sep="\t", quote=F, row.names=F)
    or = generate_odds_ratios(res)
    or2 = or %>% select(-exposure,-outcome)
    or2$exposure = e
    or2$outcome = t
    write.table(or2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_OR_with_CI95.txt"),sep="\t", quote=F, row.names=F)
    leave_ivw = mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
    leave_ivw2 = leave_ivw %>% select(-exposure,-outcome)
    leave_ivw2$exposure = e
    leave_ivw2$outcome = t
    write.table(leave_ivw2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_leaveoneout_ivw.txt"),sep="\t", quote=F, row.names=F)
    leave_wm =mr_leaveoneout(dat, parameters = default_parameters(), method = mr_weighted_median)
    leave_wm2 = leave_wm %>% select(-exposure,-outcome)
    leave_wm2$exposure = e
    leave_wm2$outcome = t
    write.table(leave_wm2, file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_leaveoneout_wm.txt"),sep="\t", quote=F, row.names=F)
    p1 <- mr_scatter_plot(res2, dat)
    length(p1)
    ggsave(p1[[1]], file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_scatter_plot.pdf"), width=7, height=7)
    res_single <- mr_singlesnp(dat, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    p2 <- mr_forest_plot(res_single)
    ggsave(p2[[1]], file=paste0(output_dir,"Univariable_MR_",e,"_",t,"_forest_plot.pdf"), width=7, height=7)
    options(warn=-1)
  }
}


#Combine MR output files into single files
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_results.txt   > Univariable_MR_dietary_traits_cvd_results.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_OR_with_CI95.txt   > Univariable_MR_dietary_traits_cvd_OR_with_CI95.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_heterogeneity.txt   > Univariable_MR_dietary_traits_cvd_heterogeneity.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_pleiotropy.txt   > Univariable_MR_dietary_traits_cvd_pleiotropy.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_ivw.txt   > Univariable_MR_dietary_traits_cvd_leaveoneout_ivw.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_leaveoneout_wm.txt   > Univariable_MR_dietary_traits_cvd_leaveoneout_wm.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_steiger_test.txt   > Univariable_MR_dietary_traits_cvd_steiger_test.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_SNP_data_table.txt   > Univariable_MR_dietary_traits_cvd_SNP_data_table.txt"))
#system(paste0("awk 'NR == 1 || FNR > 1'  Univariable_MR_*_*_singleSNP_analysis.txt   > Univariable_MR_dietary_traits_cvd_singleSNP_analysis.txt"))

#Filter out the pairs that pass Bonferroni-adjusted pval in at least 2 sensitivity analyses.
#d = read.table("Univariable_MR_dietary_traits_cvd_OR_with_CI95.txt",sep='\t',header=T)

#d2 = d %>% select(-id.exposure,-id.outcome) %>% relocate(exposure,.before=method) %>% relocate(outcome,.after=exposure) %>% filter(pval < 1.88e-4) 
#sig = d2 %>% group_by(exposure,outcome) %>% filter(n()>= 2)

#write.table(sig, file="sig_Univariable_MR_dietary_traits_cvd_OR_with_CI95.txt",sep='\t',row.names=F,quote=F)

