# HEADER ------------------------------------------------------------
#
# TITLE:   5.descriptive_stats.R
#
# PURPOSE: To better understand the available data at each timepoint 
#
# DATE:    Jan 22, 2023 
#
# INPUT:  
#     
# OUTPUT:  
#          
# SET UP -----------------------------------------------------------

# clear out working space ###
rm(list = ls())

# Load Packages ####
require(ggplot2);require(MASS);
require(stringr);require(withr);
require(phyloseq);require(nlme);library(kableExtra)
require(ggsignif);require(lsr);library(knitr)
require(dplyr);require(tidyverse);library(janitor)
library(FactoMineR);library(factoextra);library(psych)
library(ggpubr);library(ggthemes);library(extrafont);
library(ggthemes);library(rgl); library(knitr);library(tableone);
library(tibble);library(tidyverse);library(sJPlot);

# set directories 
meta2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files"
meta2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles"

# Import and edit metadata ####
mm_meta <- read.csv("/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv")

# make a list of continuous and categorical variables we are interested in: 

#(add the new ones 02/19/2023)
vars<- c("mother_age", "SES_index_final", "prepreg_bmi_kgm2", "baby_gender", "gestational_age_cat",
         "baby_birthweight_kg","breastfeedings_continuous", 'father_age', "baby_gender_cat", "gestational_age_cat", "inf_antibiotics", "mode_of_delivery_cat", "moms_ed_level",'pets','breastmilk_per_day', "martial_status", "cog_quartile", "cog_mean_sd", "cog_low2high", "mot_low2high","bsid_cog_cs","bsid_lang_cs", "bsid_mot_cs", "bsid_mot_gm_ss","bsid_mot_fm_ss","bsid_se_cs")

cont_vars <- c("mother_age", "SES_index_final", "prepreg_bmi_kgm2", "baby_gender", "gestational_age_cat",
                "baby_birthweight_kg","breastfeedings_continuous", 'father_age', 
               "bsid_cog_cs","bsid_lang_cs", "bsid_mot_cs", "bsid_mot_gm_ss","bsid_mot_fm_ss","bsid_se_cs")

cat_vars <- c("baby_gender_cat", "gestational_age_cat", "inf_antibiotics", "mode_of_delivery_cat", "moms_ed_level",'pets','breastmilk_per_day', "martial_status", "cog_quartile", "cog_mean_sd", "cog_low2high", "mot_low2high")

# Creating Table Ones ####

#filter for the time point you are interested in and divide into quartiles 
timepoint6 <- mm_meta[which(mm_meta$timepoint == "6"),]

# creating table 1 for continuous variables at tp of interest
table_6m <- CreateTableOne(vars = vars, data = timepoint6, factorVars = cat_vars)

# creating table 1 comparing between cog quartiles
table_6m_quartiles <- CreateTableOne(vars = vars,  strata =c("cog_quartile"), data = timepoint6, factorVars = cat_vars, addOverall = TRUE)
table_6m_cog_low2high <- CreateTableOne(vars = vars,  strata =c("cog_low2high"), data = timepoint6, factorVars = cat_vars, addOverall = TRUE)
# saving the tables of interest
table_6m_cog_low2high <- as.data.frame.Tableone(table_6m_cog_low2high)
write.csv(Table1_tpts, "table1_timepoints.csv")

# not working, use the janitor package 
timepoint6 %>% 
  tabyl(baby_gender_cat, cog_quartile) %>% 
  adorn_totals(c("row", "col")) %>%
  adorn_percentages("row") %>% 
  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  adorn_ns() %>% 
  kable %>% kable_classic(full_width = FALSE)

# looking at breastfeeding continuous in relation to cog_cs
ggplot(mm_meta, aes(y=bsid_mot_cs, x=breastfeedings_continuous)) +
  geom_point() + geom_smooth(method=lm, level=0.99, se=TRUE, col='blue', size=1)

ggplot(mm_meta, aes(y=bsid_lang_cs, x=breastmilk_per_day)) +
  geom_point() + geom_smooth(method=lm, level=0.99, se=TRUE, col='blue', size=1)

boxplot(mm_meta$bsid_cog_cs ~ mm_meta$breastmilk_per_day)

modelForm <- as.formula(paste("bsid_cog_cs", "~breastfeedings_continuous + SES_index_final + mode_of_delivery_cat + baby_gender + gestational_age_cat"))
modelInfo <- lm(modelForm, data = mm_meta)





