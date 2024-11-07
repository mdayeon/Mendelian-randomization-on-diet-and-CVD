## For post-hoc power calculations for all associations tested in MR diet-CVD
## Ref: https://shiny.cnsgenomics.com/mRnd/ ; https://github.com/kn3in/mRnd/blob/master/functions.R

# Provide the following:
# N = sample size (only binary outcomes)
# K = proportion of cases
# Rxz = Proportion of variance explained for the association between the SNP or allele score (Z) and the exposure variable (X)
# OR =  True odds ratio of the outcome variable per standard deviation of the exposure variable

#working directory = '/project/voight_T2D_UM1_FGP/diane531/Diet_CVD/post_hoc_power_calculation'


# Create a table of sample size and proportion of cases for 7 binary outcomes 

sample_size = data.frame(outcome = c("AAA","AF", "AS","CAD","HF","PAD","VTE"),
                         N = c(144970,1030836,446696,1165690, 977323, 174992,190266),
                         K = c(0.035, 0.059, 0.091,0.16, 0.048, 0.14, 0.047))

write.table(sample_size, file='sample_sizes_and_proportion_cases_for_7cvds.txt',sep='\t', row.names=F,quote=F)

# Create a table of 38 dietary traits and Rxz

diet = read.table('38dietary_traits.txt',sep='\t',header=T)
diet2 = diet %>% mutate(Rxz = f*(k / (n-k-1)))


#Merge both binary outcome (CVD) and dietary trait tables together into a single table
exp = diet2
out = sample_size

ass = read.table('266_exposure_outcome_associations.txt',sep='\t',header=T)

ass_exp = left_join(ass,exp2,by='exposure')
ass_exp_out = left_join(ass_exp,out,by='outcome')

write.table(ass_exp_out,file='266_exposure_outcome_associations_for_post_hoc_calculation.txt',sep='\t',row.names=F,quote=F)


# Power calculation, ORbounds for an individual association.
OR_low_high_test = function(N, K, Rxz, OR){
 
  beta_mr <- K * ( OR/ (1 + K * (OR - 1)) -1)
  
  v_mr <- (K * (1-K) - beta_mr^2) / (N * Rxz)
  
  alpha = 0.05
  threschi <- qchisq(1 - alpha, 1)
  NCP = 7.84886
  power <- 1 - pchisq(threschi, 1, NCP)
  
  c1 = N * Rxz
  c2 = K * (1 - K)
  c3 = 7.84886
  c4 = K
  c5 = sqrt((c2*c3)/(c1+c3))
  
  num = 1+(c5/c4)-c5-c4
  den = 1 - c5 - c4
  
  or_high = num/den
  print(or_high)
  
  
  x = -1 * c5 
  y = 1 + (x / K) - x - K
  z = 1 - x - c4
  or_low = y / z
  print(or_low)
  
}



# Power calculation, ORbounds for all 266 associations.
df=ass_exp_out


OR_low_high = function(df){
  exposure = df[1]
  outcome = df[2]
  N = as.numeric(df[14])
  K = as.numeric(df[15])
  Rxz = as.numeric(df[13])
  OR = as.numeric(df[7])
  
  alpha = 0.05
  threschi <- qchisq(1 - alpha, 1)
  NCP = 7.84886
  power <- 1 - pchisq(threschi, 1, NCP)

  beta_mr <- K * (OR/(1 + K * (OR-1))-1)
  
  v_mr <- (K*(1-K)-beta_mr^2)/(N*Rxz)
  
  c1 = N * Rxz
  c2 = K * (1 - K)
  c3 = 7.84886
  c4 = K
  c5 = sqrt((c2*c3)/(c1+c3))
  
  num = 1+(c5/c4)-c5-c4
  den = 1-c5-c4
  
  or_high = num/den
  
  
  x = -1*c5 
  y = 1+(x/K)-x-K
  z = 1-x-c4
  
  or_low = y / z
  
  table = data.frame(exposure = exposure,
                     binary_outcome = outcome,
                     sample_size = N,
                     proportion_cases = K,
                     proportion_variance_exposure = Rxz,
                     Non_centrality_parameter = NCP,
                     alpha = alpha, power = 0.8,
                     Odds_ratio = OR,
                     OR_high = or_high,
                     OR_low = or_low) 
  print(table)
  write.table(table,file=paste0(exposure,"_",outcome,"_power_calculation_ORbound.txt"),sep='\t',row.names=F,quote=F)
}

output = apply(df,1,OR_low_high)