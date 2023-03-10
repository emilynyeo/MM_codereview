# ------------------------------------------------------------
#
# TITLE:   6.  
# AUTHOR: Emily Yeo 
# PURPOSE: To create a taxonomic abundance plot that compares relative taxonomic abundance across infant cognitive groups.
# DATE:    March 8, 2023
# INPUT: source 6.MM_physeq_ampvis_
# OUTPUT:  
# -----------------------------------------------------------

#load directories 
#source("")
d3 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/"
d4 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/output_files/"

#load phyloseq and microshades 
library(phyloseq)
#remotes::install_github("KarstensLab/microshades")
library(microshades)
#BiocManager::install("DESeq2")
library(DESeq2)

#load phyloseq object 
ps_object <- readRDS(paste0(d3, "phyloseq_object.rds"))

#Transform the abundance data from the phyloseq object to be suitable for use with the deseq2 package:
your_deseq2_object <- phyloseq_to_deseq2(ps_object, ~ SampleID)


#FIX FROM HERE DOWN
#Normalize the abundance data in your phyloseq object to relative abundance:
your_phyloseq_object_rel <- transform_sample_counts(ps_object, function(x) x/sum(x))

#Create a MicroShades color dataframe:
color_df <- create_color_dfs(your_phyloseq_object_rel)$color_df

# use prep_mdf funct to prep microshades dataframe, combining color_df ith abundance data
mdf <- prep_mdf(your_phyloseq_object_rel, color_df)

#create abuncance plot 
microshades_plot(mdf)

