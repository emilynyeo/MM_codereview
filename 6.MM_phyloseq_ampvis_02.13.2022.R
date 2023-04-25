# Emily 11/11/23
# assessing taxonomy associated to cognitive outcomes 
# clear out working space ###
rm(list = ls())
# install.packages("BiocManager")
# BiocManager::install('phyloseq')
BiocManager::install('microbiome')
# Load Packages ####
require(ggplot2);require(MASS);
require(stringr);require(withr);
require(phyloseq);require(nlme);
require(ggsignif);require(lsr);
require(dplyr);require(tidyverse);
library(FactoMineR);library(factoextra);library(BiocManager);
library(ggpubr);library(ggthemes);library(extrafont);
library(ggthemes); library(rgl); library(microbiome);
library(tidyr); library(stringr); library(phyloseq)

# Set working directories ####
meta1 <- "/Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/"
# meta2 <- "/Volumes/ADORLab/__Users/emye7956/MM"
meta4 <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/taxonomy/"
input <- "/Volumes/IPHY/ADORLab/HEI Study/16S Data/MM 16S Processed Files for Mom and Infants 05Jan21/"
d3 <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/"
# Load Meta Data ####
meta_full <- read.csv(paste0(meta1,"3_mothersMilk_metadata_timepointsAsRows_updated051022_Temporary24mDiet.csv"), sep=",", header=TRUE)
meta_physeq <- read.csv(paste0(d3,"a_meta_physeq.csv"), sep=",", header=TRUE)
meta <- as.data.frame(meta_full)
#meta <- read.csv("/Volumes/IPHY/ADORLab/__Users/emye7956/MM/code_review/input_files/a_meta.csv") #all participants & tpts 


# Formatting Taxa table ####
taxaLevels <- c("level2","level3","level4","level5","level6")
for(taxa in taxaLevels)
{
  taxaTable <- read.delim(paste0(meta4,"mm_150bp_deblur_sepp_noMit_noChl_rmsingletons_",taxa,".txt",sep=""),header=TRUE,row.names=1,check.names = FALSE)

taxonomy_table <- data.frame(row.names(taxaTable))
taxonomy_table <- data.frame(str_split_fixed(taxonomy_table$row.names.taxaTable., "[;]", 6))
colnames(taxonomy_table) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus")

taxonomy_table$Kingdom <- substring(taxonomy_table$Kingdom, 4)
taxonomy_table$Phylum <- substring(taxonomy_table$Phylum, 4)
taxonomy_table$Class <- substring(taxonomy_table$Class, 4)
taxonomy_table$Order <- substring(taxonomy_table$Order, 4)
taxonomy_table$Family <- substring(taxonomy_table$Family, 4)
taxonomy_table$Genus <- substring(taxonomy_table$Genus, 4)
row.names(taxonomy_table) <- row.names(taxaTable)

taxonomy_table$Genus <- gsub("[.]", "", taxonomy_table$Genus)
}

# adding an otu name column

#do any rownames match? (no)
setdiff(rownames(taxonomy_table), rownames(taxaTable)) # character(0)
# are all the taxa present? 
all(rownames(taxaTable) == rownames(taxonomy_table)) # true

# adding explicit OTU ID column 
taxonomy_table$OTU_ID <- 1:226
# moving OTU id to front to make row names match
taxonomy_table2 <- taxonomy_table %>%
  dplyr::select(OTU_ID, everything())

# Add row names to be OTU id
row.names(taxonomy_table2) <- taxonomy_table2$OTU_ID

#create/specify Phyloseq taxonomy table 
taxa_physeq <- tax_table(as.matrix(taxonomy_table2))

# check for duplicates - none
duplicated(taxa_names(taxa_physeq))

# Formatting OTU table ####

# Grabbing just the infant gut microbiota at baseline
taxaTableInf <- taxaTable[,grepl("Inf.", colnames(taxaTable))] #306 infants

# Make sure to continue using only infants or else it will have duplicated sample names with moms
otu_normalized <- t(t(taxaTableInf) / colSums(taxaTableInf))

# Altering otu file IDs to match 
revisedid<-sapply(as.character(colnames(otu_normalized)),function(x){paste0(
  strsplit(x,"[.]")[[1]][2],"-",strsplit(x,"[.]")[[1]][3],"-",strsplit(x,"[.]")[[1]][5]
)})

colnames(otu_normalized) <-sapply(revisedid, function(x){
  if(strsplit(x,"-")[[1]][3]=="Bl") return(gsub("Bl","01",x))
  else if (strsplit(x,"-")[[1]][3]=="6m") return(gsub("6m","06",x))
  else return(gsub("12m","12",x))})

# Make dataframe
otu_norm2 <- as.data.frame(otu_normalized)

OTU_physeq <- otu_table(t(otu_norm2), taxa_are_rows = T)

length(which(duplicated(taxa_names(OTU_physeq)) == TRUE)) # 1
length(taxa_names(OTU_physeq)) # 306
# Duplicated taxa names (merge_dyad_id) in OTU_physeq is causing problems. Removing it below:
dupNames <- which(duplicated(taxa_names(OTU_physeq)) == TRUE)
otu_norm3 <- otu_norm2[,(-dupNames)]

# adding OTU id to match taxonomy table 
otu_norm3$OTU_ID <- 1:226
row.names(otu_norm3) <- otu_norm3$OTU_ID

# moving OTU id to front to make row names match (may not be necessary)
otu_norm4 <- otu_norm3 %>%
  dplyr::select(OTU_ID, everything())
row.names(otu_norm4) <- otu_norm4$OTU_ID

# Make sure to remove blanks 
blanks <- otu_norm3[,grepl("BLANK", colnames(otu_norm3))]
otu_norm_without_matches <- otu_norm3 %>% dplyr::select(-one_of(names(blanks)))

# Make phyloseq otu object 
OTU_physeq <- otu_table(as.matrix(otu_norm_without_matches), taxa_are_rows = T)

# Format metadata for phyloseq object ####

#add some variable of interest(can't read in modified a_meta for some reason)

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

# Make sure to metadata sample names match OTU table sample names
meta_t <- as.data.frame(t(meta))
colnames(meta_t) <- meta$merge_id_dyad

meta_3 <- meta_t %>% dplyr::select(one_of(names(otu_norm_without_matches)))
meta_4 <- as.data.frame(t(meta_3))
meta_physeq <- sample_data(meta_4)

length(sample_names(OTU_physeq)) #306 - includes OTU_ID column
length(sample_names(meta_physeq)) #305

# Checking to make sure it will run
sample_names(OTU_physeq) # gives all the dyad IDs
sample_names(taxa_physeq) # NULL - I think this is okay
sample_names(meta_physeq) # these match OTU_physeq

# Should match based on OTU_ID
taxa_names(OTU_physeq) 
taxa_names(taxa_physeq)

# Merge Phyloseq object ####
physeq_object <- phyloseq(OTU_physeq, taxa_physeq, meta_physeq)

#testing normalization (no. of "reads"/sample) to scale relative abundances to 100%
physeq100  = transform_sample_counts(physeq_object, function(x) 100 * x/sum(x))

# Testing some functions in phyloseq
sample_names(physeq_object)
rank_names(physeq_object)
sample_variables(physeq_object)

# Save Phyloseq object as a file for later use ####
file_whole <- paste0(d3, "phyloseq_object.rds")
saveRDS(physeq_object, file = file_whole)

file_100 <- paste0(d3, "phyloseq_object_abundance_tested.rds")
saveRDS(physeq100, file = file_100)

# Subsetting data - for instance, only keep 6 month timepoint
physeq_6m <- subset_samples(physeq2, timepoint ==" 6")
physeq100_6m <- subset_samples(physeq100, timepoint ==" 6")

# Creating stacked bar plots
# Creating SES index to use as an example category
physeq_6m@sam_data$SES_index_final_cat[3 <= physeq_6m@sam_data$SES_index_final & physeq_6m@sam_data$SES_index_final <= 27] = "Very Low"
physeq_6m@sam_data$SES_index_final_cat[27 < physeq_6m@sam_data$SES_index_final & physeq_6m@sam_data$SES_index_final <= 40] = "Medium"
physeq_6m@sam_data$SES_index_final_cat[40 < physeq_6m@sam_data$SES_index_final & physeq_6m@sam_data$SES_index_final <= 66] = "Very High"

# Plot bacterial phyla by abundance in SES groups
physeq_6m_SES <- merge_samples(physeq_6m, "SES_index_final_cat")
plot_bar(physeq_6m_SES, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Plot bacterial phyla by abundance in cognitive groups 

plot_bar(physeq100_6m, x = "bsid_cog_cs", fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  labs(x = "Bayleys composite score")

# testing out a heat map
plot_heatmap(physeq_6m, method = "NMDS", distance = "bray" , fill = "Phylum", sample.label = NULL, taxa.label = "Phylum", low = "#000033", high = "#66CCFF", na.value = "black")

# Ampvis              
library(ampvis2)
amp_otu <- tibble::rownames_to_column(otu_norm2, var = "OTU_ID")
rownames(amp_otu) <- amp_otu[,1]
# create ampvis object
datraw <- amp_load(otutable = otu_norm3,
                   metadata = meta_3,
                   taxonomy = taxonomy_table2)

plot.heat <- datraw %>%
  amp_heatmap(group_by = "dyad_id",
              color_vector = c("white", "royalblue4", "salmon"),
              facet_by = "mode_of_delivery",
              tax_aggregate = "Phylum",
              tax_add = NULL,
              tax_show = 5,
              plot_values_size = 2,
              plot_values = T,
              plot_na = FALSE, 
              plot_colorscale = "sqrt",
              normalise = F,
              round = 2)

#plot bargraph 
plot.bar <- datraw %>%
  group_by("dyad_id") %>%
  pivot_longer(cols = 2:53, names_to = "Sample", values_to = "RA") %>%
  ggplot(aes(x = Sample, y = RA, fill = Species)) +
  geom_bar(stat="identity") +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.5, hjust = 1))

# Beta diversity permanovas 
permanova <- vegan::adonis2(t(datraw$abund) ~ bsid_cog_cs+SES_index_final+prepreg_bmi_kgm2, 
                            data = datraw$metadata, 
                            na.action = na.exclude, 
                            permutations = 999, method = "bray")

permanova %>%
  rownames_to_column("Variable") %>%
  mutate(across(3:6, ~round(., digits = 3))) %>%
  gridExtra::tableGrob(rows = NULL) %>%
  gridExtra::grid.arrange(.)

plot.ord <- datraw %>%
  amp_ordinate(#sample_color_by = "mode_of_delivery",
               sample_shape_by = "SES_index_final",
               x_axis = 1,
               y_axis = 2,
               type = "PCA",
               distmeasure = "bray",
               transform = "log",
               filter_species = 0,
               species_plot = TRUE,
               species_label_taxonomy = "Species",
               species_nlabels = 5,
               species_label_size = 3,
               sample_colorframe = F,
               envfit_factor = c("mode_of_delivery", "SES_index_final"),
               envfit_numeric = NULL,
               envfit_arrowcolor = "Grey",
               envfit_textcolor = "Grey") +
  scale_color_viridis_d() +
  theme_dark()

anosim <- vegan::anosim(t(datraw$abund), datraw$metadata$bsid_cog_cs, permutations = 999, distance = "bray")

# Beta diversity Ellie and Bridget  -----------------------------------------
# calculate Bray-Curtis dissimilarity
dist = vegdist(otu_norm3[,1:305], distance = "BC")
# subset on columns 1-305 because the last column is sample id
adonis2(dist ~ bsid_cog_cs, data = mm_meta, permutations = 1000, na.rm = T)
# get coordinates for your distance matrix
pca <- prcomp(dist)

# plot the ordination
ggbiplot(pca, group = meta$bsid_cog_cs, ellipse = T, var.axes = F) +
  theme_classic() 