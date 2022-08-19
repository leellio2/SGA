# SGA
Code to support observational subgroup analysis

Arguments:
dat 		; name of dataframe
treatment	; vector that contains the treatment information, where 0=not treated, 1=treated
cont		; vector of the continuous variables you wish to use: e.g. c("age", "bmi", etc...)
cat		; vector of the dichotomous variables you wish to use: e.g. c("sex","CHF", etc...)
weight		; vector containing the weights
out_type	; 1=just overall, 2=just treated and not treated (DEFAULT), 3=overall and treated/not treated, will default to treated / not_treated
width		; scale of the output png, (recommended "70%")

Example:

source(file=”Simdata.r”)
baseline_table(test_df, A, contnames, catnames, null.weight, 3, "70%")
baseline_table(test_df, A, contnames, catnames, OW, 3, "70%")





