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

#load phyloseq object ####
ps_object <- readRDS(paste0(d3, "phyloseq_object.rds"))
ps_obj_ab_tested <- readRDS(paste0(d3, "phyloseq_object_abundance_tested.rds"))
source("6.MM_phyloseq_ampvis_02.13.2022.R") #using physeq_object for now

#Begin using microshades functions to evaluate abundance and apply advanced color organization at the Phylum Family level.

# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
my_mdf_prep <- prep_mdf(my_s3, subgroup_level = "Family")

# Create a color object for the specified data
my_color<- create_color_dfs(my_mdf_prep, group_level = "Phylum", subgroup_level = "Family", cvd = TRUE)

# Extract
my_mdf <- my_color$mdf
my_cdf <- my_color$cdf

# plot abundance
plot_1 <- plot_microshades(my_mdf, my_cdf, group_label = "Phylum Family")

plot_1 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) 

# reorder_samples_by will change the order of samples based on an abundance of a specified subgroup taxonomy
new_sample_order_physeqMM <- reorder_samples_by(my_mdf, my_cdf, order = "Bacteroidaceae", group_level = "Phylum", subgroup_level = "Family", sink_abundant_groups = FALSE)

mdf_new_smample_order <-new_sample_order_physeqMM$mdf
cdf_new_smample_order <-new_sample_order_physeqMM$cdf

plot_2 <- plot_microshades(mdf_new_smample_order, cdf_new_smample_order, group_label = "Phylum Family")

plot_2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6))
