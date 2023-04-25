# Alpha Plots ------------------------------------------------------------
#
## TITLE:   alpha_d_plots.R
#
# PURPOSE: To show the associations of alpha diversity metrices to cognitive outcomes 
#
# DATE:    February 19, 2023
#
# INPUT:   "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv"
#     
# OUTPUT:  Z:IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles/alpha_cogn_plots_Feb2023.pdf

# clear out working space ###
rm(list = ls())

# Load Packages ####
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(dplyr)

# Set working directories ####
d1 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files"
d2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles/alpha_d"
d3 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles\ /alpha_d/"

# Importing data ####
a_meta <- read.csv("/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv") # all participants and timepoints 
# OR read in the entire formating scrifpt using: 
# source("1.make_meta.R")

# Regression loop ####

# Plot theme for graphs to come 
mytheme <- theme(axis.line = element_line(size = 0.5, 
                                          colour = "black"), 
                 panel.background = element_rect(fill = "white"),
                 legend.position = c(.95, .95),
                 legend.justification = c("right", "top"),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6))

## Setting up exposure to loop through ####
ExposureFactors <- c( "gut_richness",
                      "shannon_entropy",
                      "faith_pd")

## setting up outcomes to loop through ####
InfantFactors <- c("bsid_se_cs", "bsid_ab_cs", 
                   "bsid_mot_cs", "bsid_cog_cs",
                   "bsid_lang_cs", "bsid_mot_gm_ss", "bsid_mot_fm_ss")

#output plots from loop to pdf
pdf(file = paste0(d3,"alpha_cogn_plots_April2023.pdf",sep=""))

# loop components 
exposureNamesList <- character(0)
#for model 1
pVal1ExposureList <- numeric(0)
beta1ExposureList <- numeric(0)
beta1SDList <- numeric(0)
beta1ScaledList <- numeric(0)
CI1ExposureLowerList <- numeric(0)
CI1ExposureUpperList <- numeric(0)
CI1UpperScaledList <- numeric(0)
CI1LowerScaledList <- numeric(0)
#for model 2
pVal2ExposureList <- numeric(0)
beta2ExposureList <- numeric(0)
beta2SDList <- numeric(0)
beta2ScaledList <- numeric(0)
CI2ExposureLowerList <- numeric(0)
CI2ExposureUpperList <- numeric(0)
CI2UpperScaledList <- numeric(0)
CI2LowerScaledList <- numeric(0)

Exposure0List <- numeric(0)
Outcome0List <- character(0)
allExposures <- ExposureFactors
allpvals <- numeric(0)
adjpvals <- numeric(0)

#Sets up a loop to run through all cognitive outcomes
for(j in seq(1:length(InfantFactors)))
{
  #For each cognitive outcome, we then loop through each of our exposures & incrementally build a linear model for each exposure/outcome pair
  for(i in seq(1:length(ExposureFactors))) 
    
  {
    
    # Choose the timepoint of interest (remember to change the loop if you change it)
    meta_6m <- a_meta %>% filter(timepoint == 6) #old run without removing 
    
    # Create dateframe without NAs (include all variables to be adjusted for in model)
    thisDataInstance = na.omit(data.frame(thisOutcome=meta_6m[,names(meta_6m) %in% InfantFactors[j]], thisExposureFactor = meta_6m[,names(meta_6m) %in% ExposureFactors[i]],
                                          sex=meta_6m$baby_gender,
                                          birthweight=meta_6m$baby_birthweight_kg,
                                          SES=meta_6m$SES_index_final,
                                          preprebmi=meta_6m$prepreg_bmi_kgm2,
                                          maternalage=meta_6m$mother_age,
                                          gestationalagebl=meta_6m$gestational_age_category,
                                          mode_of_delivery=meta_6m$mode_of_delivery))
    
    #linear model for plot
    modelFormForPlot <- as.formula(paste(InfantFactors[j],"~",ExposureFactors[i]))
    statModForPlot <- lm(modelFormForPlot, data=meta_6m)
    
    #1) Model 1 + infant sex (Full Model)
    modelForm1=as.formula(paste(InfantFactors[j],"~",ExposureFactors[i],"+",
                                "SES_index_final", "+", "breastmilk_per_day", "+", 
                                "as.factor","(","gestational_age_category",")", "+", "prepreg_bmi_kgm2", "+", 
                                "baby_birthweight_kg", "+", "as.factor","(","baby_gender",")", "+",
                                "as.factor", "(", "mode_of_delivery",")", "+", "mother_age"))
    
    #2) Model 2 (Full Model) + Infant Sex Interaction
    modelForm2=as.formula(paste(InfantFactors[j],"~",ExposureFactors[i],"+",
                                "SES_index_final", "+", "breastmilk_per_day", "+",
                                "as.factor","(","gestational_age_category",")", "+",
                                "prepreg_bmi_kgm2", "+", "baby_birthweight_kg", "+",
                                "as.factor","(","baby_gender",")",  "+",
                                "as.factor", "(", "mode_of_delivery",")", "+", "mother_age", "+", "as.factor","(","baby_gender",")",
                                "*",ExposureFactors[i]))
    
    #prints model to be viewed in console as code runs
    print(modelForm1); print(modelForm2)
    
    #Save Beta and Pvalue and put into a running list
    statMod1 <- lm(modelForm1,data=meta_6m) 
    beta1Exposure <- summary(statMod1)$coefficients[2,1] 
    pVal1Exposure <- summary(statMod1)$coefficients[2,4] 
    beta1ExposureList <- beta1ExposureList %>% append(beta1Exposure) 
    pVal1ExposureList <- pVal1ExposureList %>% append(pVal1Exposure)
    
    CI1ExposureLower <- confint(statMod1, level=0.95)[2,1] 
    CI1ExposureUpper <- confint(statMod1, level=0.95)[2,2] 
    CI1ExposureLowerList <- CI1ExposureLowerList %>% append(CI1ExposureLower) 
    CI1ExposureUpperList <- CI1ExposureUpperList %>% append(CI1ExposureUpper) 
    
    Exposure0List <- Exposure0List %>% append(InfantFactors[j])
    Outcome0List <- Outcome0List %>%  append(ExposureFactors[i])
    
    # pulling SD from 6m data set 
    beta1SD <- sd(as.numeric((meta_6m[[ExposureFactors[i]]]) %>% na.omit))
    beta1SDList <- beta1SDList %>% append(beta1SD)
    
    beta1Scaled <- (beta1Exposure * beta1SD) 
    beta1ScaledList <- beta1ScaledList %>% append(beta1Scaled %>% round(2)) 
    
    CI1UpperScaled <- (CI1ExposureUpper * beta1SD) 
    CI1UpperScaledList <- CI1UpperScaledList %>% append(CI1UpperScaled %>% round(2))
    
    CI1LowerScaled <- (CI1ExposureLower * beta1SD) 
    CI1LowerScaledList <- CI1LowerScaledList %>% append(CI1LowerScaled %>% round(2))
    
    # repeating for model 2
    statMod2 <- lm(modelForm2,data=meta_6m)
    beta2Exposure <- summary(statMod2)$coefficients[10,1] 
    pVal2Exposure <- summary(statMod2)$coefficients[10,4] 
    beta2ExposureList <- beta2ExposureList %>% append(beta2Exposure %>% round(2))
    pVal2ExposureList <- pVal2ExposureList %>% append(pVal2Exposure)
    
    CI2ExposureLower <- confint(statMod2, level=0.95)[2,1] 
    CI2ExposureUpper <- confint(statMod2, level=0.95)[2,2] 
    CI2ExposureLowerList <- CI2ExposureLowerList %>% append(CI2ExposureLower) 
    CI2ExposureUpperList <- CI2ExposureUpperList %>% append(CI2ExposureUpper) 
    
    # repeating for model 2. pulling SD from 6m data set  
    beta2SD <- sd(as.numeric((meta_6m[[ExposureFactors[i]]]) %>% na.omit))
    beta2SDList <- beta2SDList %>% append(beta2SD)
    
    beta2Scaled <- (beta2Exposure * beta2SD) 
    beta2ScaledList <- beta2ScaledList %>% append(beta2Scaled %>% round(2)) 
    
    CI2UpperScaled <- (CI2ExposureUpper * beta2SD) 
    CI2UpperScaledList <- CI2UpperScaledList %>% append(CI2UpperScaled %>% round(2))
    
    CI2LowerScaled <- (CI2ExposureLower * beta2SD) 
    CI2LowerScaledList <- CI2LowerScaledList %>% append(CI2LowerScaled %>% round(2))
    

    #Plots show beta and p value for statMod1
    title <- paste("Effect Size = ",
                   format(beta1Exposure,digits=3),
                   ", p-value=",
                   format(pVal1Exposure,digits=3), 
                   ", n = ",nrow(thisDataInstance), sep="")
    
    print(qplot(thisExposureFactor,thisOutcome, data = thisDataInstance) +
            ggtitle(title) + 
            labs(x=ExposureFactors[i], y=InfantFactors[j]) + 
            mytheme + 
            geom_smooth(method='lm'))
    
    
  }
  
}
graphics.off()

#Table of model outcomes ####

LM1 <- cbind(Exposure0List, Outcome0List, pVal1ExposureList, beta1ExposureList, CI1ExposureLowerList, CI1ExposureUpperList) %>% 
  as_tibble()

LM2 <- cbind(Exposure0List, Outcome0List, pVal2ExposureList, beta2ExposureList, CI2ExposureLowerList, CI2ExposureUpperList) %>% 
  as_tibble()

LM1scaled <- cbind(Exposure0List, Outcome0List, pVal1ExposureList, beta1ScaledList, CI1LowerScaledList, CI1UpperScaledList) %>% 
  as_tibble()

LM2scaled <- cbind(Exposure0List, Outcome0List, pVal2ExposureList, beta2ExposureList, CI2LowerScaledList, CI2UpperScaledList) %>% 
  as_tibble()

LM2scaled %>% kable() %>% kable_classic(full_width = F)

# Remove outliers ####
# https://universeofdatascience.com/how-to-remove-outliers-from-data-in-r/

#removing cog_cs outliers 
quartiles <- quantile(a_meta$bsid_cog_cs, probs=c(.25, .75), na.rm = TRUE)
IQR <- IQR(a_meta$bsid_cog_cs)
Lower <- quartiles[1] - 1.5*IQR
Upper <- quartiles[2] + 1.5*IQR 
data_cog_no_outlier <- subset(a_meta, a_meta$bsid_cog_cs > Lower & a_meta$bsid_cog_cs < Upper)
dim(data_cog_no_outlier)

# could also use the Rosners test 
library(EnvStats)
# the parameter k defines the number of potential outlier in a dataset.
rosnerTest(a_meta$gut_richness, k = 3)

# repeating the regression loop without outliers ####

#Sets up a loop to run through all cognitive outcomes
for(j in seq(1:length(InfantFactors)))
{
  #For each cognitive outcome, we then loop through each of our exposures & incrementally build a linear model for each exposure/outcome pair
  for(i in seq(1:length(ExposureFactors))) 
    
  {
    
    # Choose the timepoint of interest (remember to change the loop if you change it)
    meta_6m<- data_cog_no_outlier %>% filter(timepoint == 6)
    
    # Create dateframe without NAs (include all variables to be adjusted for in model)
    thisDataInstance = na.omit(data.frame(thisOutcome=meta_6m[,names(meta_6m) %in% InfantFactors[j]], thisExposureFactor = meta_6m[,names(meta_6m) %in% ExposureFactors[i]],
                                          sex=meta_6m$baby_gender,
                                          birthweight=meta_6m$baby_birthweight_kg,
                                          SES=meta_6m$SES_index_final,
                                          preprebmi=meta_6m$prepreg_bmi_kgm2,
                                          maternalage=meta_6m$mother_age,
                                          gestationalagebl=meta_6m$gestational_age_category,
                                          mode_of_delivery=meta_6m$mode_of_delivery))
    
    #linear model for plot
    modelFormForPlot <- as.formula(paste(InfantFactors[j],"~",ExposureFactors[i]))
    statModForPlot <- lm(modelFormForPlot, data=meta_6m)
    
    #1) Model 1 + infant sex (Full Model)
    modelForm1=as.formula(paste(InfantFactors[j],"~",ExposureFactors[i],"+",
                                "SES_index_final", "+", "breastmilk_per_day", "+", 
                                "as.factor","(","gestational_age_category",")", "+", "prepreg_bmi_kgm2", "+", 
                                "baby_birthweight_kg", "+", "as.factor","(","baby_gender",")", "+",
                                "as.factor", "(", "mode_of_delivery",")", "+", "mother_age"))
    
    #2) Model 2 (Full Model) + Infant Sex Interaction
    modelForm2=as.formula(paste(InfantFactors[j],"~",ExposureFactors[i],"+",
                                "SES_index_final", "+", "breastmilk_per_day", "+",
                                "as.factor","(","gestational_age_category",")", "+",
                                "prepreg_bmi_kgm2", "+", "baby_birthweight_kg", "+",
                                "as.factor","(","baby_gender",")",  "+",
                                "as.factor", "(", "mode_of_delivery",")", "+", "mother_age", "+", "as.factor","(","baby_gender",")",
                                "*",ExposureFactors[i]))
    
    #prints model to be viewed in console as code runs
    print(modelForm1); print(modelForm2)
    
    #Save Beta and Pvalue and put into a running list
    statMod1 <- lm(modelForm1,data=meta_6m) 
    beta1Exposure <- summary(statMod1)$coefficients[2,1] 
    pVal1Exposure <- summary(statMod1)$coefficients[2,4] 
    beta1ExposureList <- beta1ExposureList %>% append(beta1Exposure) 
    pVal1ExposureList <- pVal1ExposureList %>% append(pVal1Exposure)
    
    CI1ExposureLower <- confint(statMod1, level=0.95)[2,1] 
    CI1ExposureUpper <- confint(statMod1, level=0.95)[2,2] 
    CI1ExposureLowerList <- CI1ExposureLowerList %>% append(CI1ExposureLower) 
    CI1ExposureUpperList <- CI1ExposureUpperList %>% append(CI1ExposureUpper) 
    
    Exposure0List <- Exposure0List %>% append(InfantFactors[j])
    Outcome0List <- Outcome0List %>%  append(ExposureFactors[i])
    
    # pulling SD from 6m data set 
    beta1SD <- sd(as.numeric((meta_6m[[ExposureFactors[i]]]) %>% na.omit))
    beta1SDList <- beta1SDList %>% append(beta1SD)
    
    beta1Scaled <- (beta1Exposure * beta1SD) 
    beta1ScaledList <- beta1ScaledList %>% append(beta1Scaled %>% round(2)) 
    
    CI1UpperScaled <- (CI1ExposureUpper * beta1SD) 
    CI1UpperScaledList <- CI1UpperScaledList %>% append(CI1UpperScaled %>% round(2))
    
    CI1LowerScaled <- (CI1ExposureLower * beta1SD) 
    CI1LowerScaledList <- CI1LowerScaledList %>% append(CI1LowerScaled %>% round(2))
    
    # repeating for model 2
    statMod2 <- lm(modelForm2,data=meta_6m)
    beta2Exposure <- summary(statMod2)$coefficients[10,1] 
    pVal2Exposure <- summary(statMod2)$coefficients[10,4] 
    beta2ExposureList <- beta2ExposureList %>% append(beta2Exposure %>% round(2))
    pVal2ExposureList <- pVal2ExposureList %>% append(pVal2Exposure)
    
    CI2ExposureLower <- confint(statMod2, level=0.95)[2,1] 
    CI2ExposureUpper <- confint(statMod2, level=0.95)[2,2] 
    CI2ExposureLowerList <- CI2ExposureLowerList %>% append(CI2ExposureLower) 
    CI2ExposureUpperList <- CI2ExposureUpperList %>% append(CI2ExposureUpper) 
    
    # repeating for model 2. pulling SD from 6m data set  
    beta2SD <- sd(as.numeric((meta_6m[[ExposureFactors[i]]]) %>% na.omit))
    beta2SDList <- beta2SDList %>% append(beta2SD)
    
    beta2Scaled <- (beta2Exposure * beta2SD) 
    beta2ScaledList <- beta2ScaledList %>% append(beta2Scaled %>% round(2)) 
    
    CI2UpperScaled <- (CI2ExposureUpper * beta2SD) 
    CI2UpperScaledList <- CI2UpperScaledList %>% append(CI2UpperScaled %>% round(2))
    
    CI2LowerScaled <- (CI2ExposureLower * beta2SD) 
    CI2LowerScaledList <- CI2LowerScaledList %>% append(CI2LowerScaled %>% round(2))
    
    #Plots show beta and p value for statMod1
    title <- paste("Effect Size = ",
                   format(beta1Exposure,digits=3),
                   ", p-value=",
                   format(pVal1Exposure,digits=3), 
                   ", n = ",nrow(thisDataInstance), sep="")
    
    print(qplot(thisExposureFactor,thisOutcome, data = thisDataInstance) +
            ggtitle(title) + 
            labs(x=ExposureFactors[i], y=InfantFactors[j]) + 
            mytheme + 
            geom_smooth(method='lm'))
    
    
  }
  
}
graphics.off()

#Table of model outcomes without outliers ####

LM1 <- cbind(Exposure0List, Outcome0List, pVal1ExposureList, beta1ExposureList, CI1ExposureLowerList, CI1ExposureUpperList) %>% 
  as_tibble()

LM2 <- cbind(Exposure0List, Outcome0List, pVal2ExposureList, beta2ExposureList, CI2ExposureLowerList, CI2ExposureUpperList) %>% 
  as_tibble()

LM1scaled <- cbind(Exposure0List, Outcome0List, pVal1ExposureList, beta1ScaledList, CI1LowerScaledList, CI1UpperScaledList) %>% 
  as_tibble()

LM2scaled <- cbind(Exposure0List, Outcome0List, pVal2ExposureList, beta2ExposureList, CI2LowerScaledList, CI2UpperScaledList) %>% 
  as_tibble()

LM2scaled %>% kable() %>% kable_classic(full_width = F)
