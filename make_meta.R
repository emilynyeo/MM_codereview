# ------------------------------------------------------------
#
# TITLE:   1. metadata_formatting 
#
# PURPOSE: To create a metadata file that has variables of interest onl, additional catagories made, and alpha diversity metrics added, that can be used for subsequent downstream analysis  
#
# DATE:    February 17, 2023
#
# INPUT:   /Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/3_mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv
# OUTPUT:  Z:IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles/metadata

# -----------------------------------------------------------

# clear out working space ###
rm(list = ls())

#load packages
library(tidyverse)
library(dplyr)

#load directories ####
d1 <- "/Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/"
d2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/metadata"
d3 <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/alpha diversity/"

# read in and change variables in meta data ####
meta <- read.csv(paste0(d1,"3_mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"))

#select only columns of interest (21) and remove infants without cognitive scores (300) ####
meta <- meta %>% 
  filter(!is.na(bsid_cog_cs)) %>% 
  select('dyad_id','timepoint','bsid_cog_cs','bsid_lang_cs', 'bsid_lang_cs','bsid_mot_cs','bsid_mot_fm_ss','bsid_mot_gm_ss','bsid_se_cs', 'bsid_ab_cs', 'SES_index_final','baby_gender','baby_birthweight_kg','prepreg_bmi_kgm2','mother_age', 'father_age', 'moms_ed_level','martial_status','pets','inf_antibiotics','breastmilk_per_day','gestational_age_category','mode_of_delivery')

# create variables of interest ####
# Breast feeding per day
meta$breastfeedings_continuous <-ifelse(meta$breastmilk_per_day == 0,0, 
                                 ifelse(meta$breastmilk_per_day == 1,1, 
                                 ifelse(meta$breastmilk_per_day >1 & meta$breastmilk_per_day < 9, 
                                 meta$breastmilk_per_day-1, 8)))

# Mode of delivery
meta$mode_of_delivery_cat <- factor(ifelse(meta$mode_of_delivery %in% 1, 
                                              "Vaginal","C-Section"),
                                       levels = c("Vaginal","C-Section"))
# Baby gender cat 
meta$baby_gender_cat <- factor(ifelse(meta$baby_gender %in% 1, 
                                         "Female","Male"))  
  
# Cognitive category, low, medium, high 
hist(meta$bsid_cog_cs)
meta$cog_mean_sd <- meta %>% 
  #categorizing by the means and SD
  as.factor(ifelse(bsid_cog_cs <81.33925, 'below_mean-1SD', 
            ifelse(bsid_cog_cs <106.3654, 'within_SDs','above_mean+SD'))) %>% 

#categorizing by the <80 literature cutoff for MDI 
meta$cog_low2high <- as.factor(ifelse(meta$bsid_cog_cs < 80, "low", 
                                      ifelse(meta$bsid_cog_cs < 100, "medium", "high")))

#categorizing by quartiles 
quantile(meta$bsid_cog_cs, 
                     probs=c(0.25, 0.5, 0.75), 
                     names=TRUE,
                     na.rm = TRUE) # 25% 50% 75% = 90  95 100

meta$cog_quartile <- ifelse(ntile(n = 4, meta$bsid_cog_cs)==1, "Q1",
                               ifelse(ntile(n=4,meta$bsid_cog_cs)==2, "Q2",
                                      ifelse(ntile(n=4,meta$bsid_cog_cs)==3, "Q3",
                                             ifelse(ntile(n=4,meta$bsid_cog_cs)==4, "Q4", "NA"))))

# force these to be factors
mm_meta$baby_gender <- factor(mm_meta$baby_gender, labels = c("Female", "Male"))
mm_meta$gestational_age_cat <- factor(mm_meta$gestational_age_cat)
mm_meta$baby_antibiotics <- factor(mm_meta$baby_antibiotics)
mm_meta$cog_mean_sd <- factor(mm_meta$cog_mean_sd)


# Import Alpha diversity data
shannon_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_shannon.qza"))
richness_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_richness.qza"))
faith_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_faithspd.qza"))



