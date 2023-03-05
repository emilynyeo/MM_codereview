# HEADER ------------------------------------------------------------
#
# TITLE:   taxa_cog_associations.R
#
# PURPOSE: Assess the associatons betweeen infant cognition, at various timepoints and individual taxonomies at each taxonomic level. 
#
# DATE:    Feb 15, 2023 
#
# INPUT:   "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/all_taxa_metatable.csv"
#   OR: source("3.merging_taxa_metadata.R)
#     
# OUTPUT:  
#          
# SET UP -----------------------------------------------------------
#
# git add -A
# git commit -m "some message saying what you did"
# git push
#
# Remember to put " on each side of the message, an do NOT use ! or / or other 
# weird characters.


# clear out working space ###
rm(list = ls())
# Load Packages ####
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(dplyr)

# ReSet working directories ####
d2 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/outputfiles/"
d3 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/"
d4 <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/taxonomy/"

# Load the prior (script 3) created taxa and meta data file (needs to be updated)
taxaMeta <- read.csv(paste0(d3,"all_taxa_metatable_01.16.23.csv"))

# Descriptive statistics for taxonomy data ####

# Create susbets of data based on timepoint (visit number)
taxa_meta1 <- taxaMeta %>% filter(grepl('Inf', mom_baby)) # filtering out the moms data

taxa_meta_1m <- taxaMeta[taxaMeta$timepoint %in% "1",]  #304 observations
taxa_meta_6m <- taxaMeta[taxaMeta$timepoint %in% "6",]  #495 observations
taxa_meta_12m <- taxaMeta[taxaMeta$timepoint %in% "12",] #258 observations

# removing duplicate ID's that arose from merging taxonomy and metadat
taxa_meta_1m <- taxa_meta_1m %>% distinct(study_id, .keep_all = TRUE)
taxa_meta_6m <- taxa_meta_6m %>% distinct(study_id, .keep_all = TRUE)
taxa_meta_12m <- taxa_meta_12m %>% distinct(study_id, .keep_all = TRUE)

# Create taxa_meta the contains only IDs that had a visit at 1 and 6 or 12 or 24 months
taxa_meta_1_and_6m <- filter(taxa_meta_1m, taxa_meta_1m$merge_id_dyad %in% taxa_meta_6m$merge_id_dyad)
taxa_meta_6_and_12m <- filter(taxa_meta_6m, taxa_meta_6m$merge_id_dyad %in% taxa_meta_12m$merge_id_dyad)
# n0. of participants with taxonomy data across all time points 
taxa_meta_1_6_12m <- filter(taxa_meta_1_and_6m, taxa_meta_1_and_6m$merge_id_dyad %in% taxa_meta_12m$merge_id_dyad)

# Nested regression loop for microbiota analysis. Loops through each taxonomic level ####

# Creates a list, makes it a characters and set it to empty
exposureNamesList <- character(0)
variableList <- character(0)

# Need to ALWAYS verify this is using the corr`ect names based on our taxa_meta
taxa_start <- which(colnames(taxaMeta) =="k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces" ) 
allExposures <- names(taxaMeta)[taxa_start:(length(taxaMeta))]
variableList <- c("bsid_se_cs", "bsid_ab_cs", 
                  "bsid_mot_cs", "bsid_cog_cs",
                  "bsid_lang_cs", "bsid_mot_gm_ss", "bsid_mot_fm_ss")

# Regression loop #### (set for time point 6m)
# Diversity and each taxanomic level
taxaLevels <- c("level2","level3","level4","level5","level6")

# Start of loop for obtaining microbiota data
for(taxa in taxaLevels)
{
  taxaTable <- read.delim(paste0(d4,"mm_150bp_deblur_sepp_noMit_noChl_rmsingletons_",taxa,".txt",sep=""),header=TRUE,row.names=1,check.names = FALSE)
  
  # Grabbing just the infant gut microbiota at baseline
  taxaTableInf <- taxaTable[,grepl("Inf.", colnames(taxaTable))]
  
  myT <- t(taxaTableInf) 
  totalReads <- colSums(myT)
  
  # Normalzing sequences at each taxanomic level
  myTLogged <- log10((myT/(rowSums(myT) + 1)) * (sum(rowSums(myT))/nrow(myT)) + 1)
  
  # Filtering normalized sequences at each taxanomic level
  totalReadsFiltered <- totalReads[totalReads >= 10]
  myT3 <- myTLogged[ ,colnames(myTLogged) %in% names(totalReadsFiltered)]
  #myT3 <- myTLogged 
  
  # Filtering normalized sequences (that meet some threshold of presence across samples)
  myT4 <- myT3[ ,(colSums(myT3 == 0) / nrow(myT3)) <= 0.9] # Removing bacteria that are observed in <= 10% of samples
  myT5 <- as.data.frame(myT4)
  
  # Extracting ID number using character positions (e.g., 10 to 13)
  myT5$merge_id_dyad <- substr(row.names(myT5),10,13)
  
  # Merging the microbiota data with the metadata 
  taxaMeta <- merge(taxaMeta,myT5,by="merge_id_dyad")
  
  for(variableOfInterest in variableList) 
  {
    pValList <- numeric(0)
    exposureBetaList <- numeric(0)
    lowerList <- numeric(0)
    upperList <- numeric(0)
    
    pdf(file = paste0(d2,"taxa_associations/microbiome_",taxa,"_",variableOfInterest,"_lm_withCovariates.pdf",sep=""))
    
    # Need to ALWAYS verify this is using the correct names based on our taxa_meta
    #for(i in (length(taxa_meta_6m) + 1):(length(taxa_meta1)))
    for(i in taxa_start:(length(taxa_meta1))) # now loops from line 45 through all taxa
    {
      
      # Printing the name of the taxa we are looking at in each iteration of this nested loop
      print(names(taxa_meta1)[i])
      
      # Creating the model you will run (with covariates) - make sure to consider what x and  y should be
      modelForm <- as.formula(paste("thisMicrobe", "~thisVariable + SES + mode + sex + age")) 
      
      # Creates this data instance to use for linear regression analysis
      
      thisDataInstance <- na.omit(data.frame(merge_id_dyad=taxa_meta1$merge_id_dyad,
                                             thisMicrobe=taxa_meta1[,i],
                                             thisVariable=taxa_meta1[,names(taxa_meta1) %in% variableOfInterest],
                                             age = taxa_meta1$mother_age,
                                             BMI=taxa_meta1$mom_BMI,
                                             SES = taxa_meta1$SES_index_final,
                                             mode = taxa_meta1$mode_of_delivery_cat,
                                             sex = taxa_meta1$baby_gender_cat))
      
      # Running the linear regression analysis
      modelInfo <- lm(modelForm, data = thisDataInstance)
      
      # Creating a list of coefficent names and values
      coefs <- coef(modelInfo)
      
      # Generating the 95% CIs for 'thisMicrobe' beta coefficent 
      ci_est <- confint(modelInfo, 'thisVariable', level=0.95)
      
      # Grabbing and storing the upper and lower 95% CI values in a list
      lower_ci_est <- ci_est[1]
      lowerList[[length(lowerList)+1]] <- lower_ci_est
      upper_ci_est <- ci_est[2]
      upperList[[length(upperList)+1]] <- upper_ci_est
      
      # Grabbing and storing the beta value for "thisMicrobe" in a list
      thisBeta <- format(coefs["thisVariable"],digits=3)
      exposureBetaList[[length(exposureBetaList)+1]] <- thisBeta
      
      # Grabbing and storing the p-value for the beta for "thisMicrobe" in a list
      pVal <- format.pval(coef(summary(modelInfo))["thisVariable",4],digits=3)
      pValList[[length(pValList)+1]] <- pVal
      
      # Creating plots
      title <- paste("beta =",thisBeta,", p-val = ",pVal,sep = "")
      p <- ggplot(thisDataInstance,aes(x = thisVariable,y = thisMicrobe))
      print(p + geom_point( size = 5,color = "grey36") +
              geom_abline(slope = coefs["thisVariable"],intercept = coefs["(Intercept)"] + 
                            mean(thisDataInstance$age)*coefs["age"] +
                            mean(thisDataInstance$mode)*coefs["mode"] +
                            mean(thisDataInstance$sex)*coefs["sex"] +
                            mean(thisDataInstance$SES)*coefs["SES"],size=1.5,linetype = "dashed",colour="grey36") +
              ggtitle(title) +
              xlab(variableOfInterest) +
              ylab(paste(names(taxa_meta1)[i],sep = "")) +
              theme_classic(base_size = 30) +
              theme(axis.line = element_line(size=1),
                    axis.ticks.y = element_line(size=1),
                    axis.ticks.x = element_blank(),
                    axis.text = element_text(face="bold",size=10),
                    text = element_text(face="bold",size=20),
                    legend.position = "bottom",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5)
              ) +
              theme(axis.line.x = element_line(color="black", size = 2),
                    axis.line.y = element_line(color="black", size = 2)
              )
      )
      
    }
    graphics.off()
    adjPvals <- p.adjust(pValList,method = "BH")
    allExposures <- cbind(allExposures,exposureBetaList, lowerList, upperList, pValList,adjPvals);
    exposureNamesList <- cbind(exposureNamesList,paste(variableOfInterest,"beta",sep="_"),paste(variableOfInterest,"lower_estimate_beta",sep="_"),paste(variableOfInterest,"upper_estimate_beta",sep="_"),paste(variableOfInterest,"pval",sep="_"),paste(variableOfInterest,"adj_pval",sep="_"))
    
  }
  allExposures <- data.frame(allExposures)
  names(allExposures) <- c("Taxa",exposureNamesList)
  
write.table(allExposures,paste0(d2,"statisticalTables/momBaseline_",taxa,"_Microbiome_lm_withCovariates.txt",sep=""),quote=FALSE, sep="\t",append=FALSE, col.names=TRUE, row.names=FALSE)
  
}




