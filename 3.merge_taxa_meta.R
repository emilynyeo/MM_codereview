# HEADER ------------------------------------------------------------
#
# TITLE:   3.merging_taxa_metadata.R
#
# PURPOSE: To mereg the mothers milk meta data to taxonomic data of each level 
#
# DATE:    Feb 15, 2023 
#
# INPUT:  "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv"

#Z:/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/taxonomy/mm_150bp_deblur_sepp_noMit_noChl_rmsingletons_level6.txt
    
# OUTPUT:  all_taxa_metatable_02.16.23.csv
#          
# SET UP -----------------------------------------------------------


# clear out working space ###
rm(list = ls())

# Load Packages ####
require(ggplot2);require(MASS);
require(stringr);require(withr);
require(phyloseq);require(nlme);
require(ggsignif);require(lsr);
require(dplyr);require(tidyverse);
library(FactoMineR);library(factoextra);
library(ggpubr);library(ggthemes);library(extrafont);
library(ggthemes);library(rgl);library(knitr)

# Set working directories ####
meta4 <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/taxonomy/"
input <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/"
# depending on the computer, the file paths below might work better
#meta4 <- "/Volumes/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/taxonomy/"
#input <- "/Volumes/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/"

#### Reading in Meta and Exposure Data ####

# Importing metadata #
#long
mm_meta <- read.csv("/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv") #all participants & tpts 
# Read in sequence data #
mm_seqs <- read.delim(paste0(input,"deblur seqs/mm_150bp_deblur_sepp_noMit_noChl.txt")) # OR: 
# mm_seqs <- read.delim("/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/deblur seqs/mm_150bp_deblur_sepp_noMit_noChl.txt")

# Read in mapping file, meta data (not all included in map file), and sequences
mm_map <- read.delim(paste0(input,"metadata/11908_20181022-115423_trim.txt")) #"sample_name" as rows. # OR:
#mm_map <- read.delim("/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/metadata/11908_20181022-115423_trim.txt") 

# Altering map file IDs ####
revisedid<-sapply(as.character(mm_map$sample_name),function(x){paste0(
  strsplit(x,"[.]")[[1]][2],"-",strsplit(x,"[.]")[[1]][3],"-",strsplit(x,"[.]")[[1]][5]
)})

mm_map$merged_id_dyad<-sapply(revisedid, function(x){
  if(strsplit(x,"-")[[1]][3]=="Bl") return(gsub("Bl","01",x))
  else if (strsplit(x,"-")[[1]][3]=="6m") return(gsub("6m","06",x))
  else return(gsub("12m","12",x))})

# merging map and meta ####
meta_trim = merge(mm_meta,mm_map, by = "merge_id_dyad")

# Looking at the number of NA values (already filtered out any infants with missing cog_cs values)
meta_trim %>% summarise_all((funs(sum(is.na(.))))) 

# Altering metadata (trimmed) file IDs to match that of map file 
meta_trim$merge_id_dyad<-substr(meta_trim$merge_id_dyad,4,7)

# Looping through taxonomy ####

# Diversity and each taxanomic level
taxaLevels <- c("level2","level3","level4","level5","level6")

# Start of loop for obtaining microbiota data
for(taxa in taxaLevels)
{
  taxaTable <- read.delim(paste0(meta4,"mm_150bp_deblur_sepp_noMit_noChl_rmsingletons_",taxa,".txt",sep=""),header=TRUE,row.names=1,check.names = FALSE)
  
  # Grabbing just the infant gut microbiota 
  taxaTableInf <- taxaTable[,grepl("Inf.", colnames(taxaTable))]
  myT <- t(taxaTableInf) 
  totalReads <- colSums(myT)
  
  # Normalzing sequences at each taxanomic level
  myTLogged <- log10((myT/(rowSums(myT) + 1)) * (sum(rowSums(myT))/nrow(myT)) + 1)
  
  # Filtering normalized sequences at each taxanomic level
  totalReadsFiltered <- totalReads[totalReads >= 10]
  myT3 <- myTLogged[ ,colnames(myTLogged) %in% names(totalReadsFiltered)]
  
  # Filtering normalized sequences (that meet some threshold of presence across samples)
  myT4 <- myT3[ ,(colSums(myT3 == 0) / nrow(myT3)) <= 0.9] # Removing bacteria that are observed in <= 10% of samples
  myT5 <- as.data.frame(myT4)
  
  # Extracting ID number using character positions (e.g., 10 to 13)
  myT5$merge_id_dyad <- substr(row.names(myT5),10,13)
  
  # Merging the microbiota data with the metadata 
  taxaMeta <- merge(meta_trim,myT5,by="merge_id_dyad")
  #taxaMeta <- taxaMeta %>% relocate(timepoint, .after = dyad_id)
}

write.csv(taxaMeta, file = "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/all_taxa_metatable_02.16.23.csv", append = FALSE, sep = " ",  eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)



