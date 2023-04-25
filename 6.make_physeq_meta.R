# ------------------------------------------------------------
#
# TITLE:   1. metadata formatting for phyloseq 
#
# PURPOSE: To create a metadata file that has variables of interest only, additional catagories made, and alpha diversity metrics added, that can be used for subsequent phyloseq analysis  
#
# DATE:    March 11, 2023
#
# INPUT:   /Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/3_mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv
# OUTPUT:  Z:IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles/metadata_phyloseq

# -----------------------------------------------------------

# clear out working space ###
rm(list = ls())

#load packages
library(tidyverse)
library(dplyr)
library(qiime2R)

#load directories ####
d1 <- "/Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/"
d2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/metadata"
d3 <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/alpha diversity/"

# read in and change variables in meta data ####
meta <- read.csv(paste0(d1,"3_mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"))

#select only columns of interest (21) and remove infants without cognitive scores (300) ####
#meta <- meta %>% 
 # filter(!is.na(bsid_cog_cs)) %>% 
  #dplyr::select('merge_id_dyad','dyad_id','timepoint','bsid_cog_cs','bsid_lang_cs', 'bsid_lang_cs','bsid_mot_cs','bsid_mot_fm_ss','bsid_mot_gm_ss','bsid_se_cs', 'bsid_ab_cs', 'SES_index_final','baby_gender','baby_birthweight_kg','prepreg_bmi_kgm2','mother_age', 'father_age', 'moms_ed_level','martial_status','pets','inf_antibiotics','breastmilk_per_day','gestational_age_category','mode_of_delivery')

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


#categorizing cognition by the means and SD
meta$cog_mean_sd <- as.factor(ifelse(meta$bsid_cog_cs < 81.33925, 'below_mean-1SD', 
                                     ifelse(meta$bsid_cog_cs < 106.3654, 'within_SDs', 'above_mean+1SD')))

#categorizing cognition by the literature cutoff for MDI 
meta$cog_low2high <- as.factor(ifelse(meta$bsid_cog_cs < 70, "severe", 
                                      ifelse(meta$bsid_cog_cs < 85, "at risk", 
                                             ifelse(meta$bsid_cog_cs < 100,"mild average", "high"))))

# similarly categorizing motor scores (looking at histogram)
meta$mot_low2high <- as.factor(ifelse(meta$bsid_mot_cs < 80, "low", 
                                      ifelse(meta$bsid_mot_cs < 100, "medium", "high")))

#categorizing cognition by quartiles 
quantile(meta$bsid_cog_cs, 
         probs=c(0.25, 0.5, 0.75), 
         names=TRUE,
         na.rm = TRUE) # 25% 50% 75% = 90  95 100

meta$cog_quartile <- ifelse(ntile(n = 4, meta$bsid_cog_cs)==1, "Q1",
                            ifelse(ntile(n=4,meta$bsid_cog_cs)==2, "Q2",
                                   ifelse(ntile(n=4,meta$bsid_cog_cs)==3, "Q3",
                                          ifelse(ntile(n=4,meta$bsid_cog_cs)==4, "Q4", "NA"))))

# force these to be factors
meta$baby_gender <- factor(meta$baby_gender, labels = c("Female", "Male"))
meta$gestational_age_cat <- factor(meta$gestational_age_cat)
meta$inf_antibiotics <- factor(meta$inf_antibiotics)
meta$cog_mean_sd <- factor(meta$cog_mean_sd)
meta$cog_low2high <- factor(meta$cog_low2high)
meta$cog_quartile <- factor(meta$cog_quartile)

# Import Alpha diversity data
shannon_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_shannon.qza"))
richness_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_richness.qza"))
faith_rar10k_rep2 <- read_qza(paste0(d3,"mm_150bp_deblur_sepp_noMit_noChl_rar10k_alpha_faithspd.qza"))

# Editting SHANNON ENTROPY #### 
shannon_rar10k_rep2 <-shannon_rar10k_rep2$data %>%  #converting from qiime object to R
  as.data.frame() %>%
  rownames_to_column("sample_name")

revisedid<-sapply(as.character(shannon_rar10k_rep2$sample_name),function(x){paste0(
  strsplit(x,"[.]")[[1]][2],"-",strsplit(x,"[.]")[[1]][3],"-",strsplit(x,"[.]")[[1]][5]
)})

shannon_rar10k_rep2$revisedid<-sapply(revisedid, function(x){  # creating timepoints
  if(strsplit(x,"-")[[1]][3]=="Bl") return(gsub("Bl","01",x))
  else if (strsplit(x,"-")[[1]][3]=="6m") return(gsub("6m","06",x))
  else if(strsplit(x,"-")[[1]][3]=="12m'") return(gsub("12m","12",x))
  else return(gsub("24m","24",x))}) #there was no 24 months 

shannon_rar10k_rep2 <- shannon_rar10k_rep2 %>%
  separate(revisedid, c(NA,"IDnum",NA), sep = "-")
shannon_rar10k_rep2$ID=as.numeric(str_remove(shannon_rar10k_rep2$IDnum, "^0+"))

# removing moms from data 
shannon_rar10k_rep2 <-shannon_rar10k_rep2[!grepl("Mom", shannon_rar10k_rep2$sample_name),]  

#converting from long to wide data 
shannon_rar10k_rep2 <- shannon_rar10k_rep2 %>% data.frame()
shannon_rar10k_rep2$timepoint[(grepl("Bl", shannon_rar10k_rep2$sample_name, fixed = TRUE))] = "1"
shannon_rar10k_rep2$timepoint[(grepl("6m", shannon_rar10k_rep2$sample_name, fixed = TRUE))] = "6"
shannon_rar10k_rep2$timepoint[(grepl("12m", shannon_rar10k_rep2$sample_name, fixed = TRUE))] = "12" #there is no 24+ months

shannon_rar10k_rep2 <- shannon_rar10k_rep2 %>% rename("dyad_id" = "IDnum")
shannon_rar10k_rep2$timepoint <-  shannon_rar10k_rep2$timepoint %>% as.numeric()

# Editing RICHNESS data ####

richness_rar10k_rep2 <-richness_rar10k_rep2$data %>%
  as.data.frame() %>%
  rownames_to_column("sample_name")

revisedid<-sapply(as.character(richness_rar10k_rep2$sample_name),function(x){paste0(
  strsplit(x,"[.]")[[1]][2],"-",strsplit(x,"[.]")[[1]][3],"-",strsplit(x,"[.]")[[1]][5]
)})

richness_rar10k_rep2$revisedid<-sapply(revisedid, function(x){
  if(strsplit(x,"-")[[1]][3]=="Bl") return(gsub("Bl","01",x))
  else if (strsplit(x,"-")[[1]][3]=="6m") return(gsub("6m","06",x))
  else if(strsplit(x,"-")[[1]][3]=="12m'") return(gsub("12m","12",x))
  else return(gsub("24m","24",x))})

richness_rar10k_rep2 <- richness_rar10k_rep2 %>%
  separate(revisedid, c(NA,"IDnum",NA), sep = "-")
richness_rar10k_rep2$ID=as.numeric(str_remove(richness_rar10k_rep2$IDnum, "^0+"))

# removing moms from richness data 
richness_rar10k_rep2 <-richness_rar10k_rep2[!grepl("Mom", richness_rar10k_rep2$sample_name),]  

#converting from long to wide data 
richness_rar10k_rep2 <- richness_rar10k_rep2 %>% data.frame()
richness_rar10k_rep2$timepoint[(grepl("Bl", richness_rar10k_rep2$sample_name, fixed = TRUE))] = "1"
richness_rar10k_rep2$timepoint[(grepl("6m", richness_rar10k_rep2$sample_name, fixed = TRUE))] = "6"
richness_rar10k_rep2$timepoint[(grepl("12m", richness_rar10k_rep2$sample_name, fixed = TRUE))] = "12"

richness_rar10k_rep2 <- richness_rar10k_rep2 %>% rename("dyad_id" = "IDnum")
richness_rar10k_rep2$timepoint <-  richness_rar10k_rep2$timepoint %>% as.numeric()

# Editing FAITH PD ####

faith_rar10k_rep2 <-faith_rar10k_rep2$data %>%  #converting from qiime object to R
  as.data.frame() %>%
  rownames_to_column("sample_name")
revisedid<-sapply(as.character(faith_rar10k_rep2$sample_name),function(x){paste0(
  strsplit(x,"[.]")[[1]][2],"-",strsplit(x,"[.]")[[1]][3],"-",strsplit(x,"[.]")[[1]][5])})

faith_rar10k_rep2$revisedid<-sapply(revisedid, function(x){
  if(strsplit(x,"-")[[1]][3]=="Bl") return(gsub("Bl","01",x))
  else if (strsplit(x,"-")[[1]][3]=="6m") return(gsub("6m","06",x))
  else return(gsub("12m","12",x))})

faith_rar10k_rep2 <- faith_rar10k_rep2 %>%
  separate(revisedid, c(NA,"IDnum",NA), sep = "-")
faith_rar10k_rep2$ID=as.numeric(str_remove(faith_rar10k_rep2$IDnum, "^0+"))

# removing moms from faith pd data 
faith_rar10k_rep2 <-faith_rar10k_rep2[!grepl("Mom", faith_rar10k_rep2$sample_name),]  

# creating timepoints for faith pd data
faith_rar10k_rep2 <- faith_rar10k_rep2 %>% data.frame()
faith_rar10k_rep2$timepoint[(grepl("Bl", faith_rar10k_rep2$sample_name, fixed = TRUE))] = "1"
faith_rar10k_rep2$timepoint[(grepl("6m", faith_rar10k_rep2$sample_name, fixed = TRUE))] = "6"
faith_rar10k_rep2$timepoint[(grepl("12m", faith_rar10k_rep2$sample_name, fixed = TRUE))] = "12"

faith_rar10k_rep2 <- faith_rar10k_rep2 %>% rename("dyad_id" = "ID")
faith_rar10k_rep2$timepoint <- faith_rar10k_rep2$timepoint %>% as.numeric()

# Subsetting alpha diversities #### 
faith <- faith_rar10k_rep2 %>% dplyr::select(c(dyad_id,timepoint,faith_pd)) 
shannon <- shannon_rar10k_rep2 %>% dplyr::select(c(dyad_id,timepoint,shannon_entropy))
richness <- richness_rar10k_rep2 %>% dplyr::select(c(dyad_id,timepoint,observed_features))

# Modify dyad_id in alpha-diversity dataframes to numeric 
faith$dyad_id <- faith$dyad_id %>% as.numeric()
shannon$dyad_id <- shannon$dyad_id %>% as.numeric()
richness$dyad_id <- richness$dyad_id %>% as.numeric()

faith$timepoint <- faith$timepoint %>% as.numeric()
shannon$timepoint <- shannon$timepoint %>% as.numeric()
richness$timepoint <- richness$timepoint %>% as.numeric() #there has to be a neater way...

# Join metadata and alpha diversity data #### 
a_meta <- meta %>% left_join(shannon, by=c("dyad_id", "timepoint")) %>%
  left_join(richness,by=c("dyad_id", "timepoint")) %>% 
  left_join(faith,by=c("dyad_id", "timepoint"))
#you'll get warning messages here because some time points don't have alpha measures

a_meta_physeq <- a_meta %>% rename("gut_richness" = "observed_features")

which(duplicated(a_meta_physeq$merge_id_dyad))
# must have unique ID's
#a_meta <- a_meta %>% 
#filter(dyad_id %in% unique(a_meta$dyad_id))

# output trimmed data ####

write.csv(a_meta_physeq, "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/\aa_meta_physeq.csv")
