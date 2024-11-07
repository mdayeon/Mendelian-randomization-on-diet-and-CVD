# Mendelian-randomization-on-diet-and-CVD
Mendelian randomization approach to assess causality in diet and 7 common cardiovascular diseases (CVDs). 
38 Dietary preferences as exposures and 7 CVDs as outcomes

The MR approach consists of univariable and multivariable including two-step mediation analyses. 
For univariable MR, TwoSampleMR R package was used. 
For multivariable MR, MVMR along with TwoSampleMR packages were used. 
For two-step MR mediation, RMediation R package was used to estimate indirect (mediated) effect and calculate confidence intervals.

To conduct these analyses the following codes/scripts and tables of genetic instruments were used.

UVMR_diet_CVD_TwoSampleMR.R

UVMR_MR-RAPS_diet_CVD_TwoSampleMR.R

MVMR_diet_CVD_using_both_mvmr_2samplemr.R

MVMR_MR-RAPS_diet_CVD_using_both_mvmr_2samplemr.R

Post_hoc_power_analysis_using_mRnd.R

Table of genetic instruments for MR-RAPS analysis

Table of genetic instruments for UVMR and MVMR analyses
