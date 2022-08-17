###################### Setting working directory ######################
setwd("/Users/tiagopereira/Documents/doliolid")

###################### Load libraries and additional functions/packages ######################
# Install qiime2R
# download.file("https://github.com/jbisanz/qiime2R/archive/master.zip", "source.zip")
# unzip("source.zip")
# install.packages("qiime2R-master", repos = NULL, type="source")

# To install ggh4x
# install.packages("remotes")
# remotes::install_github("teunbrand/ggh4x")

#Load libraries
library("qiime2R")
library("ape")
library("Biostrings")
library("biomformat")
library("phyloseq")
library("Hmisc")
library("yaml")
library("tidyr")
library("stats")
library("utils")
library("ggplot2")      
library("readxl")       
library("metagMisc")    
library("vegan")
library("remotes")
library("mctoolsr")
library("RColorBrewer")
library("viridis")
library("decontam")     
library("plyr")
library("ALDEx2")
library("tibble")
library("gplots")
#library("dplyr")       # the dplyr may cause conflicts with plyr and related functions
library("ComplexHeatmap")
library("ggpubr")
library("ggh4x")
library("FSA")
library("car")
library("microbiome")
library("devtools")
library("pairwiseAdonis")
library("MicEco")
library("ggVennDiagram")
library("rlist")
library("circlize")
library("RColorBrewer")
library("radiant.data")
library("randomcoloR")

# Load sources
source("scr/miseqR.R")
source("scr/summary.R")

###################### Set plotting theme to blank ######################
theme_set(theme_bw())

###################### Import all input files ######################
# Import ASV table, taxonomy table, and sample table in excel format
# These were uploaded as txt in github. Instead, use read.table to import files
asv_tab <- read_excel("data/FEFE_data.xlsx", sheet = "feature-table") # 3706 observations (AVSs) and 116 varibles (115 samples + ASV IDs)
tax_tab <- read_excel("data/FEFE_data.xlsx", sheet = "taxonomy") # 3706 observations (taxa) and 8 variables (taxonomic ranks)
sample_tab <- read_excel("data/FEFE_data.xlsx", sheet = "sample-table") # 115 observations (samples) and 17 variables (sample factors)

# Convert ASV and taxonomy tables to matrix format; sample table as data frame
asv_tab <- as.matrix(asv_tab)
class(asv_tab) # check object class/type
tax_tab <- as.matrix(tax_tab)
sample_tab <- as.data.frame(sample_tab)

# Rename row names using ASV column for ASV table and taxonomy table; for sample table rename using Sample column
# Takes all rows "," and uses first column as row names
rownames(asv_tab) <- asv_tab[,1]
rownames(tax_tab) <- tax_tab[,1]
rownames(sample_tab) <- sample_tab[,1]

# Exclude first column from all three tables; set ASV table as numeric
asv_tab <- asv_tab[,-1]
class(asv_tab) <- "numeric"
tax_tab <- tax_tab[,-1]
sample_tab <- sample_tab[,-1]
class(sample_tab)

# Transform each matrix to a phyloseq component
ASV = otu_table(asv_tab, taxa_are_rows = TRUE) # note taxa (AVSs) are as rows
TAX = tax_table(tax_tab)
samples = sample_data(sample_tab)

# Create a phyloseq object by merging all three phyloseq components.
# Other components such as a phylogentic tree and representative sequences can also be added.
phy <- phyloseq(ASV, TAX, samples)
phy # 3706 taxa and 115 samples

# Check first few rows of each matrix associated with the phyloseq object
head(sample_data(phy)) # it prints the first six rows and all columns
head(otu_table(phy))
head(tax_table(phy))

#Visualize variables associated with each phyloseq object
sample_names(phy) # name of each sample
rank_names(phy) # all taxonomic ranks from tax_table
sample_variables(phy) # variables from sample data

# Remove low abundance samples if needed. Only two low abundance samples, both negative controls. These will be removed later.
# Check what is the minimum total abundance across all samples in the dataset
smin <- min(sample_sums(phy)) # 72 reads
phy_prune <- prune_samples(sample_sums(phy)>0, phy) # since there are no samples with a total count of 0, nothing changes
smin <- min(sample_sums(phy_prune))
phy_prune # 3706 taxa and 115 samples

# Check for ASVs that have no reads
any(taxa_sums(phy_prune) == 0) # if TRUE, then there are ASVs with 0 reads, otherwise FALSE
sum(taxa_sums(phy_prune) == 0) # gives the number of cases, 0 if none
phy_prune <-prune_taxa(taxa_sums(phy_prune) > 0, phy_prune)
phy_prune # 3706 taxa and 115 samples

# Inspect Library sizes
df_phy_prune <- as.data.frame(sample_data(phy_prune)) # transform sample data into a dataframe
df_phy_prune$LibrarySize <- sample_sums(phy_prune) # creates a new variable and add the sum of sample
df_phy_prune <- df_phy_prune[order(df_phy_prune$LibrarySize),] # reorder the dataframe from low to high abundace
df_phy_prune$Index <- seq(nrow(df_phy_prune)) # create a new variable and fill it with sequential numbers

# Plot library size, i.e., the total abundance found in each sample
# Differentiate true sample vs. control samples
ggplot(data = df_phy_prune, aes(x= Index, y= LibrarySize, color = SampleControl)) +
        geom_point() +
ylab("Library Size") +
labs(color = "Sample Type")
ggsave("results/16S_FEPE_Library-size-control-sample.pdf", width = 6, height = 4, dpi = 150)

# Differentiate sample types (i.e., EG, FG, FP2Hrs, FP24Hrs, SW, and controls)
ggplot(data = df_phy_prune, aes(x= Index, y= LibrarySize, color = SampleType)) +
  geom_point() +
  ylab("Library Size") +
  labs(color = "Sample Type")
ggsave("results/16S_FEPE_Library-size-sample-type.pdf", width = 6, height = 4, dpi = 150)

# Until this point we haven't done any sort of filtering.
# We will proceed with decontam (Davis et al. 2018. Microbiome), which can be applied directly on our phyloseq object.
# Identify Contaminants using package decontam - Frequency method, i.e., based on sample DNA concentrations.
contamdf_freq <- isContaminant(phy_prune, method="frequency", conc="Qubit") # it creates a dataframe using the phyloseq object and the variable Qubit fromt the sample data
class(contamdf_freq) # check class
head(contamdf_freq) # print first six rows
tail(contamdf_freq) # print last six rows

# Count how many FALSEs and TRUEs were recovered with the frequency method
table(contamdf_freq$contaminant) # False: 3596, True: 110
head(which(contamdf_freq$contaminant)) # list first six true contaminants by their row number.
contamdf_freq_true <- contamdf_freq[grep("TRUE", contamdf_freq$contaminant), ] # list which AVSs are considered contaminants
write.csv(contamdf_freq_true, "exported_tables/contamdf_freq_true.csv") # write dataframe to file

# Explore AVSs that were recovered as potential contaminants
# In this case, rows 13 and 52 are provided
plot_frequency(phy_prune, taxa_names(phy_prune)[c(13, 52)], conc="Qubit") + 
        xlab("DNA Concentration (ng/ul)")

# Randomly choosing potential contaminants, a total of four 
set.seed(100)
plot_frequency(phy_prune, taxa_names(phy_prune)[sample(which(contamdf_freq$contaminant),4)], conc="Qubit") + 
        xlab("DNA Concentration (ng/ul)")

# Identify Contaminants - Prevalence method
sample_data(phy_prune)$is.neg <- sample_data(phy_prune)$SampleControl == "Control" # creates a new column "is.neg" and assign TRUE and FALSE for controls and samples, respectively.
tail(sample_data(phy_prune)) # check the last six rows
contamdf_prev <- isContaminant(phy_prune, method="prevalence", neg="is.neg") # apply the isContaminant function using "is.neg" variable
table(contamdf_prev$contaminant) # False: 3705, True: 1
head(which(contamdf_prev$contaminant))
contamdf_prev_true <- contamdf_prev[grep("TRUE", contamdf_prev$contaminant), ] # list which AVSs are considered contaminants
write.csv(contamdf_prev_true, "exported_tables/contamdf_prev_true.csv") # write dataframe to file

# Identify Contaminants - Prevalence -  larger threshold
contamdf_prev05 <- isContaminant(phy_prune, method="prevalence", neg="is.neg", threshold=0.5) # same as previous command, but with larger threshold.
table(contamdf_prev05$contaminant) # False: 3702, True: 4
contamdf_prev05_true <- contamdf_prev05[grep("TRUE", contamdf_prev05$contaminant), ] # list which AVSs are considered contaminants
write.csv(contamdf_prev05_true, "exported_tables/contamdf_prev05_true.csv") # write dataframe to file

# Identify Contaminants using package decontam -  Combined method (i.e., Frequency and Prevalence)
contamdf_comb <- isContaminant(phy_prune, conc="Qubit", neg="is.neg", threshold=0.1, detailed = TRUE, normalize = TRUE, method="combined") # combines prevalence and frequency
table(contamdf_comb$contaminant) # False: 3671, True: 35
head(which(contamdf_comb$contaminant))
contamdf_comb_true <- contamdf_comb[grep("TRUE", contamdf_comb$contaminant), ] # list which AVSs are considered contaminants
write.csv(contamdf_comb_true, "exported_tables/contamdf_comb_true.csv") # write dataframe to file

# Identify Contaminants using package decontam -  Combined method (i.e., Frequency and Prevalence) with larger threshold
contamdf_comb05 <- isContaminant(phy_prune, conc="Qubit", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(contamdf_comb05$contaminant) # False: 3276, True: 430
head(which(contamdf_comb05$contaminant))
contamdf_comb05_true <- contamdf_comb05[grep("TRUE", contamdf_comb05$contaminant), ] # list which AVSs are considered contaminants
write.csv(contamdf_comb05_true, "exported_tables/contamdf_comb05_true.csv") # write dataframe to file

# Make phyloseq object of presence-absence in negative controls and true samples
phy_prune_pa <- transform_sample_counts(phy_prune, function(abund) 1*(abund>0)) # here we are just changing the abundances of the phyloseq object to 0s and 1s
phy_prune_pa_neg <- prune_samples(sample_data(phy_prune_pa)$SampleControl == "Control", phy_prune_pa) # keep only control samples
phy_prune_pa_pos <- prune_samples(sample_data(phy_prune_pa)$SampleControl == "Sample", phy_prune_pa) # keep only true samples

# Make data.frame of prevalence in positive and negative samples
# Here we combine the count of each ASV from real samples (phy_prune_pa_pos), from control samples (phy_prune_pa_neg)
# and whether or not this AVS is a potential contaminant according to a specific decontam method (contamdf_prev$contaminant)
df_phy_prune_pa <- data.frame(pa.pos=taxa_sums(phy_prune_pa_pos), pa.neg=taxa_sums(phy_prune_pa_neg),
                              contaminant=contamdf_prev$contaminant)
ggplot(data=df_phy_prune_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
        xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  labs(color = "Contaminant")
ggsave("results/16S_FEPE_Prevalence01.pdf", width = 6, height = 4, dpi = 150) # only 1 potential contaminant

df_phy_prune_pa_05 <- data.frame(pa.pos=taxa_sums(phy_prune_pa_pos), pa.neg=taxa_sums(phy_prune_pa_neg),
                              contaminant=contamdf_prev05$contaminant)
ggplot(data=df_phy_prune_pa_05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
labs(color = "Contaminant")
ggsave("results/16S_FEPE_Prevalence05.pdf", width = 6, height = 4, dpi = 150) # 4 potential contaminants

df_phy_prune_pa_freq <- data.frame(pa.pos=taxa_sums(phy_prune_pa_pos), pa.neg=taxa_sums(phy_prune_pa_neg),
                                   contaminant=contamdf_freq$contaminant)
ggplot(data=df_phy_prune_pa_freq, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  labs(color = "Contaminant")
ggsave("results/16S_FEPE_Freq.pdf", width = 6, height = 4, dpi = 150) #  110 potential contaminants

df_phy_prune_pa_comb_05 <- data.frame(pa.pos=taxa_sums(phy_prune_pa_pos), pa.neg=taxa_sums(phy_prune_pa_neg),
                                   contaminant=contamdf_comb05$contaminant)
ggplot(data=df_phy_prune_pa_comb_05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  labs(color = "Contaminant")
ggsave("results/16S_FEPE_com05.pdf", width = 6, height = 4, dpi = 150) #  430 potential contaminants

# Remove contaminants from phyloseq object - frequency, prevalance, and combined based-methods, then check/remove low abundance samples
# For the doliolid dataset, the prevalence with a threshold = 0.5 filtering method is the most suitable
phy_prune_noncontam_prev05 <- prune_taxa(!contamdf_prev05$contaminant, phy_prune) # creates a new phyloseq object
phy_prune_noncontam_prev05 # 3702 taxa and 115 samples
smin <- min(sample_sums(phy_prune_noncontam_prev05)) # smin 0
phy_prune_noncontam_prev05 <- prune_samples(sample_sums(phy_prune_noncontam_prev05)>0, phy_prune_noncontam_prev05)
smin <- min(sample_sums(phy_prune_noncontam_prev05)) # smin 22
phy_prune_noncontam_prev05 # 3702 taxa and 114 samples

# After filtering potential contaminants, subset by sample types
# Remove positive/negative controls, and at the same time pruning taxa/samples with zero reads
phy_prune_noncontam_prev05_true <- subset_samples(phy_prune_noncontam_prev05, SampleTarget =="Microbiome") %>% # keeps only real samples
  prune_taxa(taxa_sums(.) > 0, .) %>% # removes any AVS with 0 counts
  prune_samples(sample_sums(.) > 0, .) # removes any sample with 0 counts
phy_prune_noncontam_prev05_true # 3695 taxa and 111 samples

# Check if code worked properly,
# that is, for ASVs and samples with 0 count
any(taxa_sums(phy_prune_noncontam_prev05_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_prev05_true) == 0) # gives the number of cases
any(sample_sums(phy_prune_noncontam_prev05_true) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_prune_noncontam_prev05_true) == 0) # gives the number of cases
smin <- min(sample_sums(phy_prune_noncontam_prev05_true)) # smin = 2175

# Subsetting samples according to Experiments
# by using the filtered phyloseq object (i.e., phy_prune_noncontam_prev05_true)
phy_prune_noncontam_prev05_true_exp1 <- subset_samples(phy_prune_noncontam_prev05_true, Experiment =="Experiment1") %>% # only exp1 kept
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove AVSs with 0 count
  prune_samples(sample_sums(.) > 0, .) # remove sample with 0 count
phy_prune_noncontam_prev05_true_exp1 # 1043 taxa and 24 samples
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp1)) # smin = 3022

phy_prune_noncontam_prev05_true_exp2 <- subset_samples(phy_prune_noncontam_prev05_true, Experiment =="Experiment2") %>% # only exp2 kept
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_exp2 # 1856 taxa and 45 samples
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp2)) # smin = 7833

phy_prune_noncontam_prev05_true_exp3 <- subset_samples(phy_prune_noncontam_prev05_true, Experiment =="Experiment3") %>% # only exp3 kept
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_exp3 # 1872 taxa and 42 samples
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp3)) #smin = 2175

######################### Ordinations ######################### 
# Ordination with nMDS without filtering
# All samples included
# All three experiments together
set.seed(1)
phy_prune_nmds <- ordinate(
  physeq = phy_prune,
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1761449 

# Ordination with nMDS all samples included
# Only potential contaminants filtered. Still has mitochondria and chloroplasts
# All three experiments together
set.seed(1)
phy_prune_noncontam_prev05_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05, 
  method = "NMDS", 
  distance = "bray"
)# Stress: 0.1660181

# Ordination with nMDS with filtering of potential contaminants
# Excluding control samples
# All three experiments together
set.seed(1)
phy_prune_noncontam_prev05_true_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true, 
  method = "NMDS", 
  distance = "bray"
)# Stress: 0.1857099

##################### Ordination with nMDS by experiment #####################
# Exp1 - only true samples
set.seed(1)
phy_prune_noncontam_prev05_true_exp1_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp1, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.07592711 

# Exp2 - only true samples
set.seed(1)
phy_prune_noncontam_prev05_true_exp2_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp2, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.0743162 

# Exp3 - only true samples
set.seed(1)
phy_prune_noncontam_prev05_true_exp3_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp3, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1217803

###################################### nMDS plots with all samples including controls ######################################
#nb.cols_all <- 7 # this can be use to set a specific number of colors
#colors_all <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_all) # uses the previous variable and extract that number of colors from the color palette
# Set colors manually for the different sample types including controls
colors_all <- c("#1B9E77", "#D95F02", "#7570B3", "#D59D08", "#66A61E", "#E7298A", "#666666")

# Plot NMDS without any filtering
# All three experiments together
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune,
  ordination = phy_prune_nmds,
  type = "samples",
  color = "SampleType",
  shape = "Experiment") +
  scale_color_manual(values = colors_all,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "NC", "PC", "SW")) +
  scale_shape_manual(values = c(15, 17, 19, 8),
                     name = "Experiment",
                     breaks = c("Experiment1", "Experiment2", "Experiment3", "no.data"),
                     labels = c("Experiment 1", "Experiment 2", "Experiment 3", "Controls"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 2, y = 2.5, label ="2D Stress: 0.18")
ggsave("results/16S_FEPE_ASV_no-filtering_nmds_all_samples.pdf", width = 8, height = 4, dpi = 150)

# Plot NMDS with filtering of potential contaminants
# All three experiments together
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_prev05,
  ordination = phy_prune_noncontam_prev05_nmds,
  type = "samples",
  color = "SampleType",
  shape = "Experiment") +
  scale_color_manual(values = colors_all,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "NC", "PC", "SW")) +
  scale_shape_manual(values = c(15, 17, 19, 8),
                     name = "Experiment",
                     breaks = c("Experiment1", "Experiment2", "Experiment3", "no.data"),
                     labels = c("Experiment 1", "Experiment 2", "Experiment 3", "Controls"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 2, y = 2.5, label ="2D Stress: 0.17") # in this plot one negative is eliminated
  ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_all_samples.pdf", width = 8, height = 4, dpi = 150)

###################################### nMDS plots excluding controls ######################################
#nb.cols_no_controls <- 5
#colors_no_controls <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_no_controls)
colors_no_controls <- c("#1B9E77", "#D95F02", "#7570B3", "#D59D08", "#666666")

# Plot NMDS with filtering of potential contaminants
# All three experiments together and only true samples
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_prev05_true,
  ordination = phy_prune_noncontam_prev05_true_nmds,
  type = "samples",
  color = "SampleType",
  shape = "Experiment") +
  facet_wrap(~ factor(Experiment, levels = c("Experiment1", "Experiment2", "Experiment3"),
                      labels = c("Experiment 1", "Experiment 2", "Experiment 3")),
                        scales = "free") +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(15, 17, 19),
                     name = "Experiment",
                     breaks = c("Experiment1", "Experiment2", "Experiment3"),
                     labels = c("Experiment 1", "Experiment 2", "Experiment 3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # remove the background of titles
  #annotate("text", x = 2, y = 3.5, label ="2D Stress: 0.19")
  ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_no_controls_facet_by_experiment.pdf", width = 8, height = 4, dpi = 150)

###################################### nMDS plots by experiment prevalence decontam ######################################
# Exp1 - only true samples
theme_set(theme_bw())
nMDS_exp1_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp1,
    ordination = phy_prune_noncontam_prev05_true_exp1_nmds,
    type = "samples",
    shape = factor("PcrReplicate"),
    color = "SampleType") +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.08")
nMDS_exp1_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_no_controls_exp1.pdf", width = 6, height = 4, dpi = 150)

# Exp2 - only true samples
theme_set(theme_bw())
nMDS_exp2_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp2,
    ordination = phy_prune_noncontam_prev05_true_exp2_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.07")
nMDS_exp2_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_no_controls_exp2.pdf", width = 6, height = 4, dpi = 150)

# Exp3 - only true samples
theme_set(theme_bw())
nMDS_exp3_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp3,
    ordination = phy_prune_noncontam_prev05_true_exp3_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.12")
nMDS_exp3_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_no_controls_exp3.pdf", width = 6, height = 4, dpi = 150)

# Combine experiments in one single figure
nMDS_experiments_prev05 <- ggarrange(nMDS_exp1_prev05, nMDS_exp2_prev05, nMDS_exp3_prev05,
                              labels = c("A", "B", "C"), nrow = 1, ncol = 3, align = "hv",
                              legend = "bottom", common.legend = TRUE)
nMDS_experiments_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_all_three_experiments.pdf", width = 10, height = 4, dpi = 150) # save graphic 

######################################### Filter out chloroplasts and mitochondria from ASV table #########################################
phy_prune_noncontam_prev05_true_filt <- phy_prune_noncontam_prev05_true %>% # creates a new phyloseq object
  subset_taxa(
    Family  != "Mitochondria" & # AVSs with mitochondria assignments are eliminated 
      Order   != "Chloroplast" # AVSs with chloroplast assignments are eliminated 
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove AVSs with 0 count
  prune_samples(sample_sums(.) > 0, .) # remove samples with 0 count
phy_prune_noncontam_prev05_true_filt # 3330 taxa and 111 samples

# Check unique taxa according to different taxonomic ranks
rank_names(phy_prune_noncontam_prev05_true_filt)
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Domain") # 2
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Phylum") # 37
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Class") # 87
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Order") # 244
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Family") # 415
get_taxa_unique(phy_prune_noncontam_prev05_true_filt, "Genus") # 836 which is the same for species level

######################################### Filter out chloroplasts and mitochondria from ASV table ###########################
############################################## From each experiment separately ############################################## 
phy_prune_noncontam_prev05_true_exp1_filt <- phy_prune_noncontam_prev05_true_exp1 %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove AVSs with 0 count
  prune_samples(sample_sums(.) > 0, .) # remove samples with 0 count

phy_prune_noncontam_prev05_true_exp1_filt # 1014 taxa and 24 samples
rank_names(phy_prune_noncontam_prev05_true_exp1_filt)
get_taxa_unique(phy_prune_noncontam_prev05_true_exp1_filt, "Phylum") # 29

phy_prune_noncontam_prev05_true_exp2_filt <- phy_prune_noncontam_prev05_true_exp2 %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove AVSs with 0 count
  prune_samples(sample_sums(.) > 0, .) # remove samples with 0 count

phy_prune_noncontam_prev05_true_exp2_filt # 1599 taxa and 45 samples
get_taxa_unique(phy_prune_noncontam_prev05_true_exp2_filt, "Phylum") # 31

phy_prune_noncontam_prev05_true_exp3_filt <- phy_prune_noncontam_prev05_true_exp3 %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>% # remove AVSs with 0 count
  prune_samples(sample_sums(.) > 0, .) # remove samples with 0 count

phy_prune_noncontam_prev05_true_exp3_filt # 1676 taxa and 42 samples
get_taxa_unique(phy_prune_noncontam_prev05_true_exp3_filt, "Phylum") # 31

# Exp1 - only true samples filt mito and chlor
set.seed(1)
phy_prune_noncontam_prev05_true_exp1_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp1_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.07518914

# Exp2 - only true samples filt mito and chlor
set.seed(1)
phy_prune_noncontam_prev05_true_exp2_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp2_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.08374383

# Exp3 - only true samples filt mito and chlor
set.seed(1)
phy_prune_noncontam_prev05_true_exp3_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp3_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1244304

# Plot nMDS with true samples filt mito and chlor
# Exp1 - only true samples
theme_set(theme_bw())
nMDS_exp1_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp1_filt,
    ordination = phy_prune_noncontam_prev05_true_exp1_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.08")
nMDS_exp1_prev05

# Exp2 - only true samples
theme_set(theme_bw())
nMDS_exp2_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp2_filt,
    ordination = phy_prune_noncontam_prev05_true_exp2_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.07")
nMDS_exp2_prev05

# Exp3 - only true samples
theme_set(theme_bw())
nMDS_exp3_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp3_filt,
    ordination = phy_prune_noncontam_prev05_true_exp3_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 0.1, y = 2, label ="2D Stress: 0.12")
nMDS_exp3_prev05

# Combine experiments in one single figure
# FINAL NMDS FOR DOLIOLID PAPER!!
nMDS_experiments_prev05_filt <- ggarrange(nMDS_exp1_prev05, nMDS_exp2_prev05, nMDS_exp3_prev05,
                                     labels = c("A", "B", "C"), nrow = 1, ncol = 3, align = "hv",
                                     legend = "bottom", common.legend = TRUE)
nMDS_experiments_prev05_filt
ggsave("results/16S_FEPE_ASV_noncontam_prev05_all_three_experiments_filt.pdf", width = 10, height = 4, dpi = 200) # save graphic 

############################################## Export phyloseq objects for Picrust2 analysis ##############################################
# Export final phyloseq ASV table according to experiment
# No contaminants, no control samples, no mitochondria and chloroplasts
asv_exp1_picrust = as(otu_table(phy_prune_noncontam_prev05_true_exp1_filt), "matrix")
asv_exp2_picrust = as(otu_table(phy_prune_noncontam_prev05_true_exp2_filt), "matrix")
asv_exp3_picrust = as(otu_table(phy_prune_noncontam_prev05_true_exp3_filt), "matrix")

# Transpose matrix, so AVSs will become columns
if(taxa_are_rows(phy_prune_noncontam_prev05_true_exp1_filt)) {asv_exp1_picrust <- t(asv_exp1_picrust)}
if(taxa_are_rows(phy_prune_noncontam_prev05_true_exp2_filt)) {asv_exp2_picrust <- t(asv_exp2_picrust)}
if(taxa_are_rows(phy_prune_noncontam_prev05_true_exp3_filt)) {asv_exp3_picrust <- t(asv_exp3_picrust)}

# Coerce to a data frame
asv_exp1_picrust_df = as.data.frame(asv_exp1_picrust)
asv_exp2_picrust_df = as.data.frame(asv_exp2_picrust)
asv_exp3_picrust_df = as.data.frame(asv_exp3_picrust)

# Write ASV_df matrices to file as csv
write.csv(asv_exp1_picrust_df, "picrust/asv_exp1_noncontam_prev05_true_exp1_filt_picrust_df.csv")
write.csv(asv_exp2_picrust_df, "picrust/asv_exp2_noncontam_prev05_true_exp1_filt_picrust_df.csv")
write.csv(asv_exp3_picrust_df, "picrust/asv_exp3_noncontam_prev05_true_exp1_filt_picrust_df.csv")

# Write TAX_df matrices to file as csv
tax_exp1_picrust = as(tax_table(phy_prune_noncontam_prev05_true_exp1_filt), "matrix")
tax_exp1_picrust_df = as.data.frame(tax_exp1_picrust)

tax_exp2_picrust = as(tax_table(phy_prune_noncontam_prev05_true_exp2_filt), "matrix")
tax_exp2_picrust_df = as.data.frame(tax_exp2_picrust)

tax_exp3_picrust = as(tax_table(phy_prune_noncontam_prev05_true_exp3_filt), "matrix")
tax_exp3_picrust_df = as.data.frame(tax_exp3_picrust)

write.csv(tax_exp1_picrust_df, "picrust/tax_exp1_noncontam_prev05_true_exp1_filt_picrust_df.csv")
write.csv(tax_exp2_picrust_df, "picrust/tax_exp2_noncontam_prev05_true_exp1_filt_picrust_df.csv")
write.csv(tax_exp3_picrust_df, "picrust/tax_exp3_noncontam_prev05_true_exp1_filt_picrust_df.csv")

############################################## Save clean phyloseq object for later analyses ##############################################
# In this sense, objects can be cleaned from the environment space making memory available
# All three experiments together
saveRDS(phy_prune_noncontam_prev05_true_filt, "data/phy_prune_noncontam_prev05_true_filt.RDS")
# Exp1
saveRDS(phy_prune_noncontam_prev05_true_exp1_filt, "data/phy_prune_noncontam_prev05_true_exp1_filt.RDS")  
# Exp2
saveRDS(phy_prune_noncontam_prev05_true_exp2_filt, "data/phy_prune_noncontam_prev05_true_exp2_filt.RDS")  
# Exp3
saveRDS(phy_prune_noncontam_prev05_true_exp3_filt, "data/phy_prune_noncontam_prev05_true_exp3_filt.RDS")  
############################################## Import clean phyloseq objects ############################################## 
# All three experiments together
phy_prune_noncontam_prev05_true_filt <- readRDS("data/phy_prune_noncontam_prev05_true_filt.RDS")
class(phy_prune_noncontam_prev05_true_filt) # check object class
phy_prune_noncontam_prev05_true_filt # 3330 taxa and 111 samples, same as previously

# All three experiments together - only seawater samples
phy_prune_noncontam_prev05_true_filt_seawater <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType =="Water") %>% # keeps only water samples
  prune_taxa(taxa_sums(.) > 0, .) %>% # removes any AVS with 0 counts
  prune_samples(sample_sums(.) > 0, .) # removes any sample with 0 counts
phy_prune_noncontam_prev05_true_filt_seawater # 827 taxa by 7 taxonomic ranks

# Exp1
phy_prune_noncontam_prev05_true_exp1_filt <- readRDS("data/phy_prune_noncontam_prev05_true_exp1_filt.RDS")
class(phy_prune_noncontam_prev05_true_exp1_filt)
phy_prune_noncontam_prev05_true_exp1_filt # 1014 taxa and 24 samples
# Exp2
phy_prune_noncontam_prev05_true_exp2_filt <- readRDS("data/phy_prune_noncontam_prev05_true_exp2_filt.RDS")
class(phy_prune_noncontam_prev05_true_exp2_filt)
phy_prune_noncontam_prev05_true_exp2_filt # 1599 taxa and 45 samples
# Exp3
phy_prune_noncontam_prev05_true_exp3_filt <- readRDS("data/phy_prune_noncontam_prev05_true_exp3_filt.RDS")
class(phy_prune_noncontam_prev05_true_exp3_filt)
phy_prune_noncontam_prev05_true_exp3_filt # 1676 taxa and 42 samples

#################################  Venn diagrams ################################## 
# Create a dataframe with an ASV variable using original tax_tab
tax_tab_asv <- as.data.frame(tax_table(phy_prune_noncontam_prev05_true_filt)) %>%
  rownames_to_column(var ="ASV")

# Generate Venn diagram using filtered datasets
# Exp1
v_exp1 <- ps_venn(phy_prune_noncontam_prev05_true_exp1_filt, group = "SampleType", type = "counts", plot = TRUE) # creates a plot
v_exp1_list <- ps_venn(phy_prune_noncontam_prev05_true_exp1_filt, group = "SampleType", type = "counts", plot = FALSE) # creates a list
class(v_exp1_list)
v_exp1_list_df <- ldply(v_exp1_list, data.frame) # it creates a dataframe with the lists
head(v_exp1_list_df)
names(v_exp1_list_df) <- c("Experiment", "ASV") # rename column names
v_exp1_list_df_tax <- left_join(v_exp1_list_df, tax_tab_asv)
write.csv(v_exp1_list_df_tax, "exported_tables/v_exp1_list_df_tax.csv")

# Exp2
v_exp2 <- ps_venn(phy_prune_noncontam_prev05_true_exp2_filt, group = "SampleType", type = "counts", plot = TRUE)
v_exp2_list <- ps_venn(phy_prune_noncontam_prev05_true_exp2_filt, group = "SampleType", type = "counts", plot = FALSE)
v_exp2_list_df <- ldply(v_exp2_list, data.frame)
head(v_exp2_list_df)
names(v_exp2_list_df) <- c("Experiment", "ASV")
v_exp2_list_df_tax <- left_join(v_exp2_list_df, tax_tab_asv)
write.csv(v_exp2_list_df_tax, "exported_tables/v_exp2_list_df_tax.csv")

# Exp 3
v_exp3 <- ps_venn(phy_prune_noncontam_prev05_true_exp3_filt, group = "SampleType", type = "counts", plot = TRUE)
v_exp3_list <- ps_venn(phy_prune_noncontam_prev05_true_exp3_filt, group = "SampleType", type = "counts", plot = FALSE)
v_exp3_list_df <- ldply(v_exp3_list, data.frame)
head(v_exp3_list_df)
names(v_exp3_list_df) <- c("Experiment", "ASV")
v_exp3_list_df_tax <- left_join(v_exp3_list_df, tax_tab_asv)
write.csv(v_exp3_list_df_tax, "exported_tables/v_exp3_list_df_tax.csv")

# Combine all three venn diagrams into a single figure
venn_sampletype_exp <- ggarrange(v_exp1, v_exp2, v_exp3, labels = c("A", "B", "C"), nrow = 1, ncol = 3)
venn_sampletype_exp

# Save it as a pdf file and make adjustments using Inkscape
ggsave("results/16S_FEPE_ASV_noncontam_prev05_venn_sampletype_exp.pdf", width = 10, height = 5) # save graphic 

################################ Venn diagrams from each sample group with all three experiments ################################
### Getting respective phyloseq objects
# Empty gut all three experiments
phy_prune_noncontam_prev05_true_filt_eg <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType == "EmptyGut") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_filt_eg # 944 taxa and 39 samples
unique(sample_data(phy_prune_noncontam_prev05_true_filt_eg)$SampleType) # Check outcome
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_eg)) # smin = 2075

# Full gut all three experiments
phy_prune_noncontam_prev05_true_filt_fg <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType == "FullGut") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_filt_fg # 1857 taxa and 39 samples
unique(sample_data(phy_prune_noncontam_prev05_true_filt_fg)$SampleType) # Check outcome
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_fg)) # smin = 7064

# Fecal pellet 2hr all three experiments
phy_prune_noncontam_prev05_true_filt_fp2hrs <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType == "FecalPellets2") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_filt_fp2hrs # 820 taxa and 9 samples
unique(sample_data(phy_prune_noncontam_prev05_true_filt_fp2hrs)$SampleType) # Check outcome
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_fp2hrs)) # smin = 36923

# Fecal pellet 24hr all three experiments
phy_prune_noncontam_prev05_true_filt_fp24hrs <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType == "FecalPellets24") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_filt_fp24hrs # 944 taxa and 39 samples
unique(sample_data(phy_prune_noncontam_prev05_true_filt_fp24hrs)$SampleType) # Check outcome
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_fp24hrs)) # smin = 32829

# Seawater all three experiments
phy_prune_noncontam_prev05_true_filt_sw <- subset_samples(phy_prune_noncontam_prev05_true_filt, SampleType == "Water") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .)
phy_prune_noncontam_prev05_true_filt_sw # 827 taxa and 12 samples
unique(sample_data(phy_prune_noncontam_prev05_true_filt_sw)$SampleType) # Check outcome
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_sw)) # smin = 39249

# Plotting Venn diagrams
# Exploring common ASVs across all three experiments per sample type
# Empty gut
v_all_exp_eg <- ps_venn(phy_prune_noncontam_prev05_true_filt_eg, group = "Experiment", type = "counts", plot = TRUE)
v_all_exp_eg_list <- ps_venn(phy_prune_noncontam_prev05_true_filt_eg, group = "Experiment", type = "counts", plot = FALSE)
v_all_exp_eg_list_df <- ldply(v_all_exp_eg_list, data.frame)
names(v_all_exp_eg_list_df) <- c("Experiment", "ASV")
v_all_exp_eg_list_df_tax <- left_join(v_all_exp_eg_list_df, tax_tab_asv)
write.csv(v_all_exp_eg_list_df_tax, "exported_tables/v_all_exp_eg_list_df_tax.csv")

unique(v_all_exp_eg_list_df_tax$Experiment)
v_all_exp_eg_list_df_tax_common <- subset(v_all_exp_eg_list_df_tax, Experiment == "Experiment1__Experiment2__Experiment3")
phy_prune_noncontam_prev05_true_filt_eg_asv_table <- as.data.frame(otu_table(phy_prune_noncontam_prev05_true_filt_eg)) %>%
  rownames_to_column(var ="ASV")
v_all_exp_eg_list_df_tax_common_abund <- left_join(v_all_exp_eg_list_df_tax_common, phy_prune_noncontam_prev05_true_filt_eg_asv_table)
v_all_exp_eg_list_df_tax_common_abund$Total <- rowSums(v_all_exp_eg_list_df_tax_common_abund [, c(10:48)], na.rm = TRUE)
write.csv(v_all_exp_eg_list_df_tax_common_abund, "exported_tables/v_all_exp_eg_list_df_tax_common_abund.csv")

unique(v_all_exp_eg_list_df_tax_common_abund$Order)
ggplot(v_all_exp_eg_list_df_tax_common_abund, aes(x = Experiment, y = Total/1000000, fill = Order)) + # dividing by 10^6 because data is not between 1 and 0
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # position fill adds to 100%
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  #geom_text(aes(label = ifelse(round(Total*100) >= 5, paste(round(Total*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_x_discrete(
    labels = c("Exp1 Exp2 Exp3"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Experiments") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_v_all_exp_eg_list_df_tax_common_abund_order_level.tiff", width = 4, height = 6, dpi = 150) # save graphic

v_all_exp_fg <- ps_venn(phy_prune_noncontam_prev05_true_filt_fg, group = "Experiment", type = "counts", plot = TRUE)
v_all_exp_fg_list <- ps_venn(phy_prune_noncontam_prev05_true_filt_fg, group = "Experiment", type = "counts", plot = FALSE)
v_all_exp_fg_list_df <- ldply(v_all_exp_fg_list, data.frame)
names(v_all_exp_fg_list_df) <- c("Experiment", "ASV")
v_all_exp_fg_list_df_tax <- left_join(v_all_exp_fg_list_df, tax_tab_asv)
write.csv(v_all_exp_fg_list_df_tax, "exported_tables/v_all_exp_fg_list_df_tax.csv")

unique(v_all_exp_fg_list_df_tax$Experiment)
v_all_exp_fg_list_df_tax_common <- subset(v_all_exp_fg_list_df_tax, Experiment == "Experiment1__Experiment2__Experiment3")
phy_prune_noncontam_prev05_true_filt_fg_asv_table <- as.data.frame(otu_table(phy_prune_noncontam_prev05_true_filt_fg)) %>%
  rownames_to_column(var ="ASV")
v_all_exp_fg_list_df_tax_common_abund <- left_join(v_all_exp_fg_list_df_tax_common, phy_prune_noncontam_prev05_true_filt_fg_asv_table)
v_all_exp_fg_list_df_tax_common_abund$Total <- rowSums(v_all_exp_fg_list_df_tax_common_abund [, c(10:48)], na.rm = TRUE)
write.csv(v_all_exp_fg_list_df_tax_common_abund, "exported_tables/v_all_exp_fg_list_df_tax_common_abund.csv")

v_all_exp_fp2hrs <- ps_venn(phy_prune_noncontam_prev05_true_filt_fp2hrs, group = "Experiment", type = "counts", plot = TRUE)
v_all_exp_fp2hrs_list <- ps_venn(phy_prune_noncontam_prev05_true_filt_fp2hrs, group = "Experiment", type = "counts", plot = FALSE)
v_all_exp_fp2hrs_list_df <- ldply(v_all_exp_fp2hrs_list, data.frame)
names(v_all_exp_fp2hrs_list_df) <- c("Experiment", "ASV")
v_all_exp_fp2hrs_list_df_tax <- left_join(v_all_exp_fp2hrs_list_df, tax_tab_asv)
write.csv(v_all_exp_fp2hrs_list_df_tax, "exported_tables/v_all_exp_fp2hrs_list_df_tax.csv")

unique(v_all_exp_fp2hrs_list_df_tax$Experiment)
v_all_exp_fp2hrs_list_df_tax_common <- subset(v_all_exp_fp2hrs_list_df_tax, Experiment == "Experiment1__Experiment2__Experiment3")
phy_prune_noncontam_prev05_true_filt_fp2hrs_asv_table <- as.data.frame(otu_table(phy_prune_noncontam_prev05_true_filt_fp2hrs)) %>%
  rownames_to_column(var ="ASV")
v_all_exp_fp2hrs_list_df_tax_common_abund <- left_join(v_all_exp_fp2hrs_list_df_tax_common, phy_prune_noncontam_prev05_true_filt_fp2hrs_asv_table)
v_all_exp_fp2hrs_list_df_tax_common_abund$Total <- rowSums(v_all_exp_fp2hrs_list_df_tax_common_abund [, c(10:18)], na.rm = TRUE)
write.csv(v_all_exp_fp2hrs_list_df_tax_common_abund, "exported_tables/v_all_exp_fp2hrs_list_df_tax_common_abund.csv")

v_all_exp_fp24hrs <- ps_venn(phy_prune_noncontam_prev05_true_filt_fp24hrs, group = "Experiment", type = "counts", plot = TRUE)
v_all_exp_fp24hrs_list <- ps_venn(phy_prune_noncontam_prev05_true_filt_fp24hrs, group = "Experiment", type = "counts", plot = FALSE)
v_all_exp_fp24hrs_list_df <- ldply(v_all_exp_fp24hrs_list, data.frame)
names(v_all_exp_fp24hrs_list_df) <- c("Experiment", "ASV")
v_all_exp_fp24hrs_list_df_tax <- left_join(v_all_exp_fp24hrs_list_df, tax_tab_asv)
write.csv(v_all_exp_fp24hrs_list_df_tax, "exported_tables/v_all_exp_fp24hrs_list_df_tax.csv")

unique(v_all_exp_fp24hrs_list_df_tax$Experiment)
v_all_exp_fp24hrs_list_df_tax_common <- subset(v_all_exp_fp24hrs_list_df_tax, Experiment == "Experiment1__Experiment2__Experiment3")
phy_prune_noncontam_prev05_true_filt_fp24hrs_asv_table <- as.data.frame(otu_table(phy_prune_noncontam_prev05_true_filt_fp24hrs)) %>%
  rownames_to_column(var ="ASV")
v_all_exp_fp24hrs_list_df_tax_common_abund <- left_join(v_all_exp_fp24hrs_list_df_tax_common, phy_prune_noncontam_prev05_true_filt_fp24hrs_asv_table)
v_all_exp_fp24hrs_list_df_tax_common_abund$Total <- rowSums(v_all_exp_fp24hrs_list_df_tax_common_abund [, c(10:18)], na.rm = TRUE)
write.csv(v_all_exp_fp24hrs_list_df_tax_common_abund, "exported_tables/v_all_exp_fp24hrs_list_df_tax_common_abund.csv")

v_all_exp_sw <- ps_venn(phy_prune_noncontam_prev05_true_filt_sw, group = "Experiment", type = "counts", plot = TRUE)
v_all_exp_sw_list <- ps_venn(phy_prune_noncontam_prev05_true_filt_sw, group = "Experiment", type = "counts", plot = FALSE)
v_all_exp_sw_list_df <- ldply(v_all_exp_sw_list, data.frame)
names(v_all_exp_sw_list_df) <- c("Experiment", "ASV")
v_all_exp_sw_list_df_tax <- left_join(v_all_exp_sw_list_df, tax_tab_asv)
write.csv(v_all_exp_sw_list_df_tax, "exported_tables/v_all_exp_sw_list_df_tax.csv")

unique(v_all_exp_sw_list_df_tax$Experiment)
v_all_exp_sw_list_df_tax_common <- subset(v_all_exp_sw_list_df_tax, Experiment == "Experiment1__Experiment2__Experiment3")
phy_prune_noncontam_prev05_true_filt_sw_asv_table <- as.data.frame(otu_table(phy_prune_noncontam_prev05_true_filt_sw)) %>%
  rownames_to_column(var ="ASV")
v_all_exp_sw_list_df_tax_common_abund <- left_join(v_all_exp_sw_list_df_tax_common, phy_prune_noncontam_prev05_true_filt_sw_asv_table)
v_all_exp_sw_list_df_tax_common_abund$Total <- rowSums(v_all_exp_sw_list_df_tax_common_abund [, c(10:21)], na.rm = TRUE)
write.csv(v_all_exp_sw_list_df_tax_common_abund, "exported_tables/v_all_exp_sw_list_df_tax_common_abund.csv")

venn_exp <- ggarrange(v_all_exp_eg, v_all_exp_fg, v_all_exp_fp2hrs, v_all_exp_fp24hrs, v_all_exp_sw,
                      labels = c("A", "B", "C", "D", "E"), nrow = 2, ncol = 3,
                      align = "hv")
venn_exp

# Save it as a pdf file and make adjustments using Inkscape
ggsave("results/16S_FEPE_ASV_noncontam_prev05_venn_exp.pdf", width = 12, height = 10) # save graphic 

################################# Alpha diversity #################################
######################### Generate a data.frame with alpha diversity measures #########################
# Calculate alpha-diversity measures (For plot purposes only!)
alpha_div <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_filt, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_filt, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_filt, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_prune_noncontam_prev05_true_filt),
  "SampleType" = phyloseq::sample_data(phy_prune_noncontam_prev05_true_filt)$SampleType,
  "Experiment" = phyloseq::sample_data(phy_prune_noncontam_prev05_true_filt)$Experiment)
alpha_div$Evenness <- alpha_div$Shannon/log(alpha_div$Observed)
head(alpha_div)

names(alpha_div) <- c("Observed", "Shannon", "Simpson", "Reads", "SampleType", "Experiment", "Evenness") # define column names
head(alpha_div)

# Plot alpha diversity measures for all three experiments
# Change names
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")

alpha_div$SampleType <- factor(alpha_div$SampleType,
                               level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

alpha_color <- c("#666666", "#D59D08","#D95F02", "#7570B3", "#1B9E77")

my_comparisons <- list(c("Water", "FullGut"), c("Water", "FecalPellets2"), c("Water","FecalPellets24"), c("Water", "EmptyGut"),
                       c("FullGut", "FecalPellets2"), c("FullGut","FecalPellets24"), c("FullGut", "EmptyGut"),
                       c("FecalPellets2","FecalPellets24"), c("FecalPellets2", "EmptyGut"),
                       c("FecalPellets24", "EmptyGut"))
alpha_div %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>% # to gather all diversity metrics
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>% # to accomodate in this order
  ggplot(aes(x = SampleType, y = value)) +
  geom_boxplot(outlier.color = NA) +
  stat_compare_means(comparisons = my_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 3, hide.ns = TRUE)+
  stat_compare_means(size = 3, label.y = 0)+
  geom_jitter(aes(color = SampleType), height = 0, width = .2) +
  facet_grid(metric ~ Experiment , scales = "free", labeller = labeller(Experiment = experiment.labs)) +
  scale_color_manual(values = alpha_color) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.08, 0.08)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  #theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(legend.position="none") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
ggsave("results/16S_FEPE_alpha_diversity_phy_prune_noncontam_prev05_true_filt.pdf", width = 8, height = 6, dpi = 150) # save graphic

#Summarize alpha diversity measures
summary_alpha <- alpha_div %>%
  group_by(Experiment, SampleType) %>%
  summarise(
    count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_Simpson = mean(Simpson),
    sd_Simpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads))
write.csv(summary_alpha, "exported_tables/summary_alpha.csv")

# Calculating diversity mean/SE using summary function
summary_alpha_all <- summarySE(alpha_div, measurevar=c("Observed"), groupvars=c("SampleType", "Experiment"))
sum_basin_sed

# Alpha-diversity measures for each experiment so KW test can be performed.
# Calculate alpha-diversity measures exp1
alpha_div_exp1 <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp1_filt, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp1_filt, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp1_filt, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_prune_noncontam_prev05_true_exp1_filt),
  "SampleType" = phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType)
alpha_div_exp1$Evenness <- alpha_div_exp1$Shannon/log(alpha_div_exp1$Observed)
head(alpha_div_exp1)

#Summarize alpha diversity measures exp1
summary_alpha_exp1 <- alpha_div_exp1 %>%
  group_by(SampleType) %>%
  summarise(
    count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(InvSimpson),
    sd_invsimpson = sd(InvSimpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads))
write.csv(summary_alpha_exp1, "exported_tables/summary_alpha_exp1.csv")

# KW analysis exp1
kw_exp1_obs <- kruskal.test(Observed ~ SampleType, data = alpha_div_exp1)
kw_exp1_h <- kruskal.test(Shannon ~ SampleType, data = alpha_div_exp1)
kw_exp1_d <- kruskal.test(InvSimpson ~ SampleType, data = alpha_div_exp1)
kw_exp1_j <- kruskal.test(Evenness ~ SampleType, data = alpha_div_exp1)
kw_exp1_n <- kruskal.test(Reads ~ SampleType, data = alpha_div_exp1)

pwt_exp1_obs <- pairwise.t.test(alpha_div_exp1$Observed, alpha_div_exp1$SampleType, p.adjust.method = "BH")
pwt_exp1_h <- pairwise.t.test(alpha_div_exp1$Shannon, alpha_div_exp1$SampleType, p.adjust.method = "BH")
pwt_exp1_d <- pairwise.t.test(alpha_div_exp1$InvSimpson, alpha_div_exp1$SampleType, p.adjust.method = "BH")
pwt_exp1_j <- pairwise.t.test(alpha_div_exp1$Evenness, alpha_div_exp1$SampleType, p.adjust.method = "BH")
pwt_exp1_n <- pairwise.t.test(alpha_div_exp1$Reads, alpha_div_exp1$SampleType, p.adjust.method = "BH")

# Calculate alpha-diversity measures exp2
alpha_div_exp2 <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp2_filt, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp2_filt, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp2_filt, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_prune_noncontam_prev05_true_exp2_filt),
  "SampleType" = phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt)$SampleType)
alpha_div_exp2$Evenness <- alpha_div_exp2$Shannon/log(alpha_div_exp2$Observed)
head(alpha_div_exp2)

#Summarize alpha diversity measures exp2
summary_alpha_exp2 <- alpha_div_exp2 %>%
  group_by(SampleType) %>%
  summarise(
    count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(InvSimpson),
    sd_invsimpson = sd(InvSimpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads))
write.csv(summary_alpha_exp2, "exported_tables/summary_alpha_exp2.csv")

# KW analysis exp2
kw_exp2_obs <- kruskal.test(Observed ~ SampleType, data = alpha_div_exp2)
kw_exp2_h <- kruskal.test(Shannon ~ SampleType, data = alpha_div_exp2)
kw_exp2_d <- kruskal.test(InvSimpson ~ SampleType, data = alpha_div_exp2)
kw_exp2_j <- kruskal.test(Evenness ~ SampleType, data = alpha_div_exp2)
kw_exp2_n <- kruskal.test(Reads ~ SampleType, data = alpha_div_exp2)

pwt_exp2_obs <- pairwise.t.test(alpha_div_exp2$Observed, alpha_div_exp2$SampleType, p.adjust.method = "BH")
pwt_exp2_h <- pairwise.t.test(alpha_div_exp2$Shannon, alpha_div_exp2$SampleType, p.adjust.method = "BH")
pwt_exp2_d <- pairwise.t.test(alpha_div_exp2$InvSimpson, alpha_div_exp2$SampleType, p.adjust.method = "BH")
pwt_exp2_j <- pairwise.t.test(alpha_div_exp2$Evenness, alpha_div_exp2$SampleType, p.adjust.method = "BH")
pwt_exp2_n <- pairwise.t.test(alpha_div_exp2$Reads, alpha_div_exp2$SampleType, p.adjust.method = "BH")

# Calculate alpha-diversity measures exp3
alpha_div_exp3 <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp3_filt, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp3_filt, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_prune_noncontam_prev05_true_exp3_filt, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_prune_noncontam_prev05_true_exp3_filt),
  "SampleType" = phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt)$SampleType)
alpha_div_exp3$Evenness <- alpha_div_exp3$Shannon/log(alpha_div_exp3$Observed)
head(alpha_div_exp3)

#Summarize alpha diversity measures exp3
summary_alpha_exp3 <- alpha_div_exp3 %>%
  group_by(SampleType) %>%
  summarise(
    count = n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(InvSimpson),
    sd_invsimpson = sd(InvSimpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads))
write.csv(summary_alpha_exp3, "exported_tables/summary_alpha_exp3.csv")

# KW analysis exp3
kw_exp3_obs <- kruskal.test(Observed ~ SampleType, data = alpha_div_exp3)
kw_exp3_h <- kruskal.test(Shannon ~ SampleType, data = alpha_div_exp3)
kw_exp3_d <- kruskal.test(InvSimpson ~ SampleType, data = alpha_div_exp3)
kw_exp3_j <- kruskal.test(Evenness ~ SampleType, data = alpha_div_exp3)
kw_exp3_n <- kruskal.test(Reads ~ SampleType, data = alpha_div_exp3)

pwt_exp3_obs <- pairwise.t.test(alpha_div_exp3$Observed, alpha_div_exp3$SampleType, p.adjust.method = "BH")
pwt_exp3_h <- pairwise.t.test(alpha_div_exp3$Shannon, alpha_div_exp3$SampleType, p.adjust.method = "BH")
pwt_exp3_d <- pairwise.t.test(alpha_div_exp3$InvSimpson, alpha_div_exp3$SampleType, p.adjust.method = "BH")
pwt_exp3_j <- pairwise.t.test(alpha_div_exp3$Evenness, alpha_div_exp3$SampleType, p.adjust.method = "BH")
pwt_exp3_n <- pairwise.t.test(alpha_div_exp3$Reads, alpha_div_exp3$SampleType, p.adjust.method = "BH")

# Anova assumptions and test -  not used for this dataset due to the lack of homogeneity on the variances
# one_way_observed_exp1 <- aov(Observed ~ SampleType, data = alpha_div_exp1)
# plot(one_way_observed_exp1, 1) # homogeneity test
# shapiro.test(alpha_div_exp1$Observed) # shapiro test of normality
# one_way_observed_exp1_res <- residuals(object = one_way_observed_exp1) # extract residuals
# shapiro.test(x = one_way_observed_exp1_res) # shapiro test of normality
# levene_observed_exp1 <- leveneTest(Observed ~ SampleType, data = alpha_div_exp1) # Test for homogeneity of variance
# tk_exp1 <- TukeyHSD(one_way_observed_exp1)

############################################## Ordination with nMDS on filtered datasets ##############################################
# All three experiments, no controls, no mitochondria/chloroplast
set.seed(1)
phy_prune_noncontam_prev05_true_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.1913417 

# Exp 1
set.seed(1)
phy_prune_noncontam_prev05_true_exp1_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp1_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.07518914

# Exp 2
set.seed(1)
phy_prune_noncontam_prev05_true_exp2_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp2_filt, 
  method = "NMDS", 
  distance = "bray"
) # 0.08374383

# Exp 3
set.seed(1)
phy_prune_noncontam_prev05_true_exp3_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_prev05_true_exp3_filt, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.1244304
############################################## PERMANOVAs ##############################################
# 1 step: data transformation
phy_prune_noncontam_prev05_true_exp1_filt_comp <- microbiome::transform(phy_prune_noncontam_prev05_true_exp1_filt, "compositional")
phy_prune_noncontam_prev05_true_exp1_filt_clr <- microbiome::transform(phy_prune_noncontam_prev05_true_exp1_filt, "clr")
phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_comp)[1:10, 1:10] # prints the first 10 rows and 10 columns, similar to command head

# 2 step: Generate distance matrix
phy_prune_noncontam_prev05_true_exp1_filt_comp_dist_matrix <- phyloseq::distance(phy_prune_noncontam_prev05_true_exp1_filt_comp, method = "bray") 
phy_prune_noncontam_prev05_true_exp1_filt_clr_dist_matrix <- phyloseq::distance(phy_prune_noncontam_prev05_true_exp1_filt_comp, method = "euclidean") 

# 3 step: ADONIS test
vegan::adonis(phy_prune_noncontam_prev05_true_exp1_filt_comp_dist_matrix ~ phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType, permutations = 9999)
vegan::adonis(phy_prune_noncontam_prev05_true_exp1_filt_clr_dist_matrix ~ phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType, permutations = 9999)

# Step 4: Dispersion test and plot
dispr_comp <- vegan::betadisper(phy_prune_noncontam_prev05_true_exp1_filt_comp_dist_matrix, phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType)
dispr_comp
dispr_clr <- vegan::betadisper(phy_prune_noncontam_prev05_true_exp1_filt_clr_dist_matrix, phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType)
dispr_clr

plot(dispr_comp, main = "", sub = "")
boxplot(dispr_comp, main = "", xlab = "")
permutest(dispr_comp, permutations = 9999)

# Step 5: Perform pairwise comparisons
pairwiseAdonis::pairwise.adonis(phy_prune_noncontam_prev05_true_exp1_filt_comp_dist_matrix, phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType, p.adjust.m = "BH")
pairwiseAdonis::pairwise.adonis(phy_prune_noncontam_prev05_true_exp1_filt_clr_dist_matrix, phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt)$SampleType, p.adjust.m = "BH")

# PCA via phyloseq
phy_prune_noncontam_prev05_true_exp1_filt_ord_clr <- phyloseq::ordinate(phy_prune_noncontam_prev05_true_exp1_filt_clr, "RDA")

# Plot scree plot
phyloseq::plot_scree(phy_prune_noncontam_prev05_true_exp1_filt_ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# Examine eigenvalues and % prop. variance explained
head(phy_prune_noncontam_prev05_true_exp1_filt_ord_clr$CA$eig)     

# Scale axes and plot ordination
phy_prune_noncontam_prev05_true_exp1_filt_ord_clr1 <- phy_prune_noncontam_prev05_true_exp1_filt_ord_clr$CA$eig[1] / sum(phy_prune_noncontam_prev05_true_exp1_filt_ord_clr$CA$eig)
phy_prune_noncontam_prev05_true_exp1_filt_ord_clr2 <- phy_prune_noncontam_prev05_true_exp1_filt_ord_clr$CA$eig[2] / sum(phy_prune_noncontam_prev05_true_exp1_filt_ord_clr$CA$eig)
phyloseq::plot_ordination(phy_prune_noncontam_prev05_true_exp1_filt, phy_prune_noncontam_prev05_true_exp1_filt_ord_clr, type="samples", color="SampleType") + 
  geom_point(size = 2) +
  coord_fixed(phy_prune_noncontam_prev05_true_exp1_filt_ord_clr2 / phy_prune_noncontam_prev05_true_exp1_filt_ord_clr1) +
  stat_ellipse(aes(group = SampleType), linetype = 2)

############################################## nMDS ############################################## 
# All three experiments together
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_prev05_true_filt,
  ordination = phy_prune_noncontam_prev05_true_filt_nmds,
  type = "samples",
  color = "SampleType",
  shape = factor("PcrReplicate")) +
  facet_wrap(~ Experiment, scales = "free", labeller = labeller(Experiment = experiment.labs)) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  #theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = 2, y = 3.5, label ="2D Stress: 0.19") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) 
  ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_filt_mitochondria_no_controls_all_experiments.pdf", width = 7, height = 3, dpi = 150)

########################## nMDS by experiments ########################## 
# Exp 1
theme_set(theme_bw())
nMDS_exp1_filt_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp1_filt,
    ordination = phy_prune_noncontam_prev05_true_exp1_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 1, y = 2, label ="2D Stress: 0.08")
nMDS_exp1_filt_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_filt_mitochondria_no_controls_exp1.pdf", width = 6, height = 4, dpi = 150)

# Exp 2
theme_set(theme_bw())
nMDS_exp2_filt_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp2_filt,
    ordination = phy_prune_noncontam_prev05_true_exp2_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 1, y = 2, label ="2D Stress: 0.08")
nMDS_exp2_filt_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_filt_mitochondria_no_controls_exp2.pdf", width = 6, height = 4, dpi = 150)

# Exp 3
theme_set(theme_bw())
nMDS_exp3_filt_prev05 <-
  phyloseq::plot_ordination(
    physeq = phy_prune_noncontam_prev05_true_exp3_filt,
    ordination = phy_prune_noncontam_prev05_true_exp3_filt_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
  scale_color_manual(values = colors_no_controls,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  annotate("text", x = 0.5, y = 2, label ="2D Stress: 0.12")
nMDS_exp3_filt_prev05
ggsave("results/16S_FEPE_ASV_noncontam_freq_nmds_filt_mitochondria_no_controls_exp3.pdf", width = 6, height = 4, dpi = 150)

# Combine all three experiments
nMDS_experiments_filt_prev05 <- ggarrange(nMDS_exp1_filt_prev05, nMDS_exp2_filt_prev05, nMDS_exp3_filt_prev05,
                                   labels = c("A", "B", "C"), nrow = 1, ncol = 3, legend = "bottom", common.legend = TRUE)
nMDS_experiments_filt_prev05
ggsave("results/16S_FEPE_ASV_noncontam_prev05_nmds_all_three_experiments_filt_mitochondria.pdf", width = 10, height = 4, dpi = 150) # save graphic

############################### Plotting top 20 most abundant at different taxonomic ranks ###############################
# This function can be used to split the ASV matrix into the top N and the less abundant taxa
merge_less_than_top_prev05_top20 <- function(phy_prune_noncontam_prev05_true_filt, top=19){ # defines what are the top taxa, in this case 19
  transformed <- transform_sample_counts(phy_prune_noncontam_prev05_true_filt, function(x) x/sum(x)) # transform abundances of phyloseq object, i.e., between 0-1
  otu.table <- as.data.frame(otu_table(transformed)) # converts ASV matrix into a dataframe
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),] # Sort matrix based of mean values, highest to lowest
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"}
  }
  return(merged)
}

# Same function but for top 12 only
merge_less_than_top_prev05_top12 <- function(phy_prune_noncontam_prev05_true_filt, top=11){
  transformed <- transform_sample_counts(phy_prune_noncontam_prev05_true_filt, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"}
  }
  return(merged)
}

# Same function but for top 8 only
merge_less_than_top_prev05_top8 <- function(phy_prune_noncontam_prev05_true_filt, top=8){
  transformed <- transform_sample_counts(phy_prune_noncontam_prev05_true_filt, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"}
  }
  return(merged)
}

# Same function but for top 4 only
merge_less_than_top_prev05_top4 <- function(phy_prune_noncontam_prev05_true_filt, top=4){
  transformed <- transform_sample_counts(phy_prune_noncontam_prev05_true_filt, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"}
  }
  return(merged)
}
# Collapse phyloseq objects at different taxonomic ranks, from phylum to genus
# Phylum to 20
phy_prune_noncontam_prev05_true_filt_top_phylum <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Phylum")
phy_prune_noncontam_prev05_true_filt_top_phylum_top20 <- merge_less_than_top_prev05_top20(phy_prune_noncontam_prev05_true_filt_top_phylum, top=19)
phy_prune_noncontam_prev05_true_filt_top_phylum_top20_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_phylum_top20) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr = aggregate(Abundance~Experiment+SampleType+Phylum, data=phy_prune_noncontam_prev05_true_filt_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr$Phylum)
# Phylum to 12
phy_prune_noncontam_prev05_true_filt_top_phylum_top12 <- merge_less_than_top_prev05_top12(phy_prune_noncontam_prev05_true_filt_top_phylum, top=11)
phy_prune_noncontam_prev05_true_filt_top_phylum_top12_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_phylum_top12) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr = aggregate(Abundance~Experiment+SampleType+Phylum, data=phy_prune_noncontam_prev05_true_filt_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr$Phylum)
# Phylum to 5
phy_prune_noncontam_prev05_true_filt_top_phylum_top5 <- merge_less_than_top_prev05_top4(phy_prune_noncontam_prev05_true_filt_top_phylum, top=4)
phy_prune_noncontam_prev05_true_filt_top_phylum_top5_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_phylum_top5) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr = aggregate(Abundance~Experiment+SampleType+Phylum, data=phy_prune_noncontam_prev05_true_filt_top_phylum_top5_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr$Phylum)
# Class to 20
phy_prune_noncontam_prev05_true_filt_top_class <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Class")
phy_prune_noncontam_prev05_true_filt_top_class_top20 <- merge_less_than_top_prev05_top20(phy_prune_noncontam_prev05_true_filt_top_class, top=19)
phy_prune_noncontam_prev05_true_filt_top_class_top20_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_class_top20) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_class_top20_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class, data=phy_prune_noncontam_prev05_true_filt_top_class_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_class_top20_agr$Class)
# Class to 12
phy_prune_noncontam_prev05_true_filt_top_class <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Class")
phy_prune_noncontam_prev05_true_filt_top_class_top12 <- merge_less_than_top_prev05_top12(phy_prune_noncontam_prev05_true_filt_top_class, top=11)
phy_prune_noncontam_prev05_true_filt_top_class_top12_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_class_top12) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_class_top12_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class, data=phy_prune_noncontam_prev05_true_filt_top_class_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_class_top12_agr$Class)
# Class to 8
phy_prune_noncontam_prev05_true_filt_top_class <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Class")
phy_prune_noncontam_prev05_true_filt_top_class_top8 <- merge_less_than_top_prev05_top8(phy_prune_noncontam_prev05_true_filt_top_class, top=7)
phy_prune_noncontam_prev05_true_filt_top_class_top8_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_class_top8) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_class_top8_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class, data=phy_prune_noncontam_prev05_true_filt_top_class_top8_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$Class)
# Order to 20
phy_prune_noncontam_prev05_true_filt_top_order <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Order")
phy_prune_noncontam_prev05_true_filt_top_order_top20 <- merge_less_than_top_prev05_top20(phy_prune_noncontam_prev05_true_filt_top_order, top=19)
phy_prune_noncontam_prev05_true_filt_top_order_top20_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_order_top20) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_order_top20_agr = aggregate(Abundance~Experiment+SampleType+Order, data=phy_prune_noncontam_prev05_true_filt_top_order_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_order_top20_agr$Order)
# Order to 12
phy_prune_noncontam_prev05_true_filt_top_order <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Order")
phy_prune_noncontam_prev05_true_filt_top_order_top12 <- merge_less_than_top_prev05_top12(phy_prune_noncontam_prev05_true_filt_top_order, top=11)
phy_prune_noncontam_prev05_true_filt_top_order_top12_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_order_top12) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_order_top12_agr = aggregate(Abundance~Experiment+SampleType+Order, data=phy_prune_noncontam_prev05_true_filt_top_order_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_order_top12_agr$Order)
# Family to 20
phy_prune_noncontam_prev05_true_filt_top_family <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Family")
phy_prune_noncontam_prev05_true_filt_top_family_top20 <- merge_less_than_top_prev05_top20(phy_prune_noncontam_prev05_true_filt_top_family, top=19)
phy_prune_noncontam_prev05_true_filt_top_family_top20_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_family_top20) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_family_top20_agr = aggregate(Abundance~Experiment+SampleType+Family, data=phy_prune_noncontam_prev05_true_filt_top_family_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_family_top20_agr$Family)
# Family to 12
phy_prune_noncontam_prev05_true_filt_top_family <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Family")
phy_prune_noncontam_prev05_true_filt_top_family_top12 <- merge_less_than_top_prev05_top12(phy_prune_noncontam_prev05_true_filt_top_family, top=11)
phy_prune_noncontam_prev05_true_filt_top_family_top12_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_family_top12) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_family_top12_agr = aggregate(Abundance~Experiment+SampleType+Family, data=phy_prune_noncontam_prev05_true_filt_top_family_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$Family)
# Genus to 20
phy_prune_noncontam_prev05_true_filt_top_genus <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Genus")
phy_prune_noncontam_prev05_true_filt_top_genus_top20 <- merge_less_than_top_prev05_top20(phy_prune_noncontam_prev05_true_filt_top_genus, top=19)
phy_prune_noncontam_prev05_true_filt_top_genus_top20_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_genus_top20) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr = aggregate(Abundance~Experiment+SampleType+Genus, data=phy_prune_noncontam_prev05_true_filt_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr$Genus)
# Genus to 12
phy_prune_noncontam_prev05_true_filt_top_genus <- tax_glom(phy_prune_noncontam_prev05_true_filt, "Genus")
phy_prune_noncontam_prev05_true_filt_top_genus_top12 <- merge_less_than_top_prev05_top12(phy_prune_noncontam_prev05_true_filt_top_genus, top=11)
phy_prune_noncontam_prev05_true_filt_top_genus_top12_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_genus_top12) # Melt to long format
phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr = aggregate(Abundance~Experiment+SampleType+Genus, data=phy_prune_noncontam_prev05_true_filt_top_genus_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr$Genus)

############################# Just top 20, no "Others" category being used #############################
# In this case, the abundance of the top taxa do not add to the real total (i.e., 100%)
# When plotting, this can be adjusted with the parameter position = "fill" of the barplot/ggplot2
# Currently, I am not using this strategy
#phy_prune_noncontam_freq_true_filt_top_phylum_20 = tax_glom(phy_prune_noncontam_freq_true_filt, taxrank = "Phylum") # agglomerate at phylum level 
#phy_prune_noncontam_freq_true_filt_top_phylum_20_sorted <- names(sort(taxa_sums(phy_prune_noncontam_freq_true_filt_top_phylum_20), decreasing=TRUE)[1:20]) # sort the top 20 Phyla
#phy_prune_noncontam_freq_true_filt_top_phylum_20_trans = transform_sample_counts(phy_prune_noncontam_freq_true_filt_top_phylum_20, function(x) x/sum(x)) # Transform to rel. abundance
#phy_prune_noncontam_freq_true_filt_top_phylum_20_prune = prune_taxa(phy_prune_noncontam_freq_true_filt_top_phylum_20_sorted, phy_prune_noncontam_freq_true_filt_top_phylum_20_trans) # subtract the top20 taxa from the entire phyloseq object
#phy_prune_noncontam_freq_true_filt_top_phylum_20_prune_df = psmelt(phy_prune_noncontam_freq_true_filt_top_phylum_20_prune) # Melt to long format
#phy_prune_noncontam_freq_true_filt_top_phylum_20_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum, data=phy_prune_noncontam_freq_true_filt_top_phylum_20_prune_df, FUN=mean) # add common factors to use for plotting
############################# Just top 12, no "Others" category being used #############################
#phy_prune_noncontam_prev05_true_filt_top_phylum_12 = tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Phylum") # agglomerate at phylum level
#phy_prune_noncontam_prev05_true_filt_top_phylum_12_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_top_phylum_12), decreasing=TRUE)[1:12]) # sort the top 12 Phyla
#phy_prune_noncontam_prev05_true_filt_top_phylum_12_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_top_phylum_12, function(x) x/sum(x)) # Transform to rel. abundance
#phy_prune_noncontam_prev05_true_filt_top_phylum_12_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_top_phylum_12_sorted, phy_prune_noncontam_prev05_true_filt_top_phylum_12_trans) # subtract the top12 taxa from the entire phyloseq object
#phy_prune_noncontam_prev05_true_filt_top_phylum_12_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_phylum_12_prune) # Melt to long format
#phy_prune_noncontam_prev05_true_filt_top_phylum_12_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum, data=phy_prune_noncontam_prev05_true_filt_top_phylum_12_prune_df, FUN=mean) # add common factors to use for plotting

phy_prune_noncontam_prev05_true_filt_top_class_5 = tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Class") #  agglomerate at class level
phy_prune_noncontam_prev05_true_filt_top_class_5_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_top_class_5), decreasing=TRUE)[1:5])  #sort the top 5 Classes
phy_prune_noncontam_prev05_true_filt_top_class_5_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_top_class_5, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_filt_top_class_5_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_top_class_5_sorted, phy_prune_noncontam_prev05_true_filt_top_class_5_trans)  #subtract the top12 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_filt_top_class_5_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_top_class_5_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class, data=phy_prune_noncontam_prev05_true_filt_top_class_5_prune_df, FUN=mean)  #add common factors to use for plotting

phy_prune_noncontam_prev05_true_filt_eg_top_fam_12 = tax_glom(phy_prune_noncontam_prev05_true_filt_eg, taxrank = "Family")
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12), decreasing=TRUE)[1:12])  #sort the top 5 Classes
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_sorted, phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_trans)  #subtract the top12 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_df, FUN=mean)  #add common factors to use for plotting


phy_prune_noncontam_prev05_true_filt_eg_top_fam_12 = tax_glom(phy_prune_noncontam_prev05_true_filt_eg, taxrank = "Family")
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12), decreasing=TRUE)[1:12])  #sort the top 5 Classes
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_sorted, phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_trans)  #subtract the top12 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_df, FUN=mean)  #add common factors to use for plotting

############################### succession phyloseq object and Plot barcharts 16S  ##############################################################
## All experiments
phy_prune_noncontam_prev05_true_filt_succession <- subset_samples(phy_prune_noncontam_prev05_true_filt,
                                                                  SampleType == c("FullGut", "FecalPellets24", "FecalPellets2"))
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_filt_succession)) # 7064 reads
# Check for ASVs that have no reads
any(taxa_sums(phy_prune_noncontam_prev05_true_filt_succession) == 0) # if TRUE, then there are ASVs with 0 reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_prev05_true_filt_succession) == 0) # gives the number of cases, 0 if none
phy_prune_noncontam_prev05_true_filt_succession <-prune_taxa(taxa_sums(phy_prune_noncontam_prev05_true_filt_succession) > 0, phy_prune_noncontam_prev05_true_filt_succession)
phy_prune_noncontam_prev05_true_filt_succession # 1311 taxa and 20 samples

#Exp 1
phy_prune_noncontam_prev05_true_exp1_filt_succession <- subset_samples(phy_prune_noncontam_prev05_true_exp1_filt,
                                                                  SampleType == c("FullGut", "FecalPellets24", "FecalPellets2"))
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp1_filt_succession)) # 7064 reads
# Check for ASVs that have no reads
any(taxa_sums(phy_prune_noncontam_prev05_true_exp1_filt_succession) == 0) # if TRUE, then there are ASVs with 0 reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_prev05_true_exp1_filt_succession) == 0) # gives the number of cases, 0 if none
phy_prune_noncontam_prev05_true_exp1_filt_succession <-prune_taxa(taxa_sums(phy_prune_noncontam_prev05_true_exp1_filt_succession) > 0, phy_prune_noncontam_prev05_true_exp1_filt_succession)
phy_prune_noncontam_prev05_true_exp1_filt_succession # 336 taxa and 4 samples

#Exp 2
phy_prune_noncontam_prev05_true_exp2_filt_succession <- subset_samples(phy_prune_noncontam_prev05_true_exp2_filt,
                                                                       SampleType == c("FullGut", "FecalPellets24", "FecalPellets2"))
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp2_filt_succession)) # 7064 reads
# Check for ASVs that have no reads
any(taxa_sums(phy_prune_noncontam_prev05_true_exp2_filt_succession) == 0) # if TRUE, then there are ASVs with 0 reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_prev05_true_exp2_filt_succession) == 0) # gives the number of cases, 0 if none
phy_prune_noncontam_prev05_true_exp2_filt_succession <-prune_taxa(taxa_sums(phy_prune_noncontam_prev05_true_exp2_filt_succession) > 0, phy_prune_noncontam_prev05_true_exp2_filt_succession)
phy_prune_noncontam_prev05_true_exp2_filt_succession # 663 taxa and 8 samples

#Exp 3
phy_prune_noncontam_prev05_true_exp3_filt_succession <- subset_samples(phy_prune_noncontam_prev05_true_exp3_filt,
                                                                       SampleType == c("FullGut", "FecalPellets24", "FecalPellets2"))
smin <- min(sample_sums(phy_prune_noncontam_prev05_true_exp3_filt_succession)) # 7064 reads
# Check for ASVs that have no reads
any(taxa_sums(phy_prune_noncontam_prev05_true_exp3_filt_succession) == 0) # if TRUE, then there are ASVs with 0 reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_prev05_true_exp3_filt_succession) == 0) # gives the number of cases, 0 if none
phy_prune_noncontam_prev05_true_exp3_filt_succession <-prune_taxa(taxa_sums(phy_prune_noncontam_prev05_true_exp3_filt_succession) > 0, phy_prune_noncontam_prev05_true_exp3_filt_succession)
phy_prune_noncontam_prev05_true_exp3_filt_succession # 675 taxa and 8 samples

# tax glom at family level
# All experiments
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4 = tax_glom(phy_prune_noncontam_prev05_true_filt_succession, taxrank = "Family")
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4), decreasing=TRUE)[1:4])  #sort the top 4 families
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_sorted, phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_trans)  #subtract the top4 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune_df, FUN=mean)  #add common factors to use for plotting

# All experiments, sewater samples only
phy_prune_noncontam_prev05_true_filt_seawater_fam = tax_glom(phy_prune_noncontam_prev05_true_filt_seawater, taxrank = "Family")
phy_prune_noncontam_prev05_true_filt_seawater_fam_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_filt_seawater_fam), decreasing=TRUE)[1:40])  #sort the top 20 families
phy_prune_noncontam_prev05_true_filt_seawater_fam_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_filt_seawater_fam, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_filt_seawater_fam_prune = prune_taxa(phy_prune_noncontam_prev05_true_filt_seawater_fam_sorted, phy_prune_noncontam_prev05_true_filt_seawater_fam)  #subtract the top20 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_filt_seawater_fam_prune_df = psmelt(phy_prune_noncontam_prev05_true_filt_seawater_fam_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_filt_seawater_fam_prune_agr = aggregate(Abundance~Experiment+SampleType+Description+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_filt_seawater_fam_prune_df, FUN=mean)  #add common factors to use for plotting

#Exp1
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4 = tax_glom(phy_prune_noncontam_prev05_true_exp1_filt_succession, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4), decreasing=TRUE)[1:4])  #sort the top 4 families
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune = prune_taxa(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_sorted, phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_trans)  #subtract the top4 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_df = psmelt(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_df, FUN=mean)  #add common factors to use for plotting

#Exp2
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4 = tax_glom(phy_prune_noncontam_prev05_true_exp2_filt_succession, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4), decreasing=TRUE)[1:4])  #sort the top 4 families
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune = prune_taxa(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_sorted, phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_trans)  #subtract the top4 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_df = psmelt(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_df, FUN=mean)  #add common factors to use for plotting

#Exp3
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4 = tax_glom(phy_prune_noncontam_prev05_true_exp3_filt_succession, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_sorted <- names(sort(taxa_sums(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4), decreasing=TRUE)[1:4])  #sort the top 4 families
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_trans = transform_sample_counts(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4, function(x) x/sum(x))  #Transform to rel. abundance
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune = prune_taxa(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_sorted, phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_trans)  #subtract the top4 taxa from the entire phyloseq object
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_df = psmelt(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune)#  Melt to long format
phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_agr = aggregate(Abundance~Experiment+SampleType+Phylum+Class+Order+Family, data=phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_df, FUN=mean)  #add common factors to use for plotting

phy_prune_noncontam_prev05_true_all_filt_succession_top_fam_4_prune_agr <- rbind(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_agr,
                                                                                 phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_agr,
                                                                                 phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_agr)

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3") # create new names for experiment variable
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3") # assign old name to new name

succession_sample_order <- factor(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune_agr$SampleType, # reorder sample types
                       level = c("FullGut", "FecalPellets2", "FecalPellets24"),
                       labels = c("FG", "FP2Hrs", "FP24Hrs"))

succession_sample_order_all <- factor(phy_prune_noncontam_prev05_true_all_filt_succession_top_fam_4_prune_agr$SampleType, # reorder sample types
                                  level = c("FullGut", "FecalPellets2", "FecalPellets24"),
                                  labels = c("FG", "FP2Hrs", "FP24Hrs"))

succession_sample_order_exp1 <- factor(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_agr$SampleType, # reorder sample types
                                  level = c("FullGut", "FecalPellets2", "FecalPellets24"),
                                  labels = c("FG", "FP2Hrs", "FP24Hrs"))

succession_sample_order_exp2 <- factor(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_agr$SampleType, # reorder sample types
                                       level = c("FullGut", "FecalPellets2", "FecalPellets24"),
                                       labels = c("FG", "FP2Hrs", "FP24Hrs"))

succession_sample_order_exp3 <- factor(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_agr$SampleType, # reorder sample types
                                       level = c("FullGut", "FecalPellets2", "FecalPellets24"),
                                       labels = c("FG", "FP2Hrs", "FP24Hrs"))

succession_top4_family_all <- c("#6A3D9A", "#33A02C", "#4363d8", 
                                "#FB9A99", "#f58231", "#ffe119",
                                "#BF40BF", "#e6beff", "#800000") 

succession_top4_family <- c("#33A02C", "#ffe119",
                         "#BF40BF", "#e6beff")

succession_top4_family_exp1 <- c("#6A3D9A", "#4363d8", 
                                 "#FB9A99", "#e6beff")

succession_top4_family_exp2 <- c("#33A02C", "#ffe119",
                                 "#e6beff", "#800000")

succession_top4_family_exp3 <- c("#33A02C", "#f58231",
                                 "#ffe119", "#BF40BF") 

# Generate 40 distinct colors with randomcoloR
n_colors <- 40
set.seed(1983765)
palette_40 <- distinctColorPalette(n_colors)
palette_40

# Draw colors in pie chart
pie(rep(1, n_colors),
    col = palette_40,
    main = "randomcoloR Package")

# Plot seawater family
ggplot(phy_prune_noncontam_prev05_true_filt_seawater_fam_prune_agr, aes(x = Description, y = Abundance/1000, fill = Family)) + #plotting by sample
  facet_grid(. ~ Experiment, scales = "free",
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance/1000) >= 5, paste(round(Abundance/1000, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = palette_40) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  #scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  #theme(legend.position = "none")+
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_phylum_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 12, height = 6, dpi = 150) # save graphic

# Plot succession family
ggplot(phy_prune_noncontam_prev05_true_all_filt_succession_top_fam_4_prune_agr, aes(x = succession_sample_order_all, y = Abundance, group = Family)) + #plotting by sample
  geom_line(aes(linetype=Family, color = Family))+
  geom_point(aes(color = Family)) +
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  scale_color_manual(values = succession_top4_family_all) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0),
                     limits = c(0, 0.7)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_family_succession_level_all_by_experiments.pdf", width = 7, height = 3, dpi = 150) # save graphic

ggplot(phy_prune_noncontam_prev05_true_filt_succession_top_fam_4_prune_agr, aes(x = succession_sample_order, y = Abundance, group = Family)) + #plotting by sample
  geom_line(aes(linetype=Family, color = Family))+
  geom_point(aes(shape=Family, color = Family)) +
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  scale_color_manual(values = succession_top4_family) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0),
                     limits = c(0, 0.7)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_family_succession_level_all_experiments.pdf", width = 7, height = 3, dpi = 150) # save graphic

# Plot
succession_exp1 <- ggplot(phy_prune_noncontam_prev05_true_exp1_filt_succession_top_fam_4_prune_agr, aes(x = succession_sample_order_exp1, y = Abundance, group = Family)) + #plotting by sample
  geom_line(aes(linetype=Family, color = Family))+
  geom_point(aes(shape=Family, color = Family)) +
  scale_color_manual(values = succession_top4_family_exp1) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
succession_exp1
ggsave("results/16S_FEPE_barplot_family_succession_level_exp1.pdf", width = 5, height = 4, dpi = 150) # save graphic

succession_exp2 <- 
  ggplot(phy_prune_noncontam_prev05_true_exp2_filt_succession_top_fam_4_prune_agr, aes(x = succession_sample_order_exp2, y = Abundance, group = Family)) + #plotting by sample
  geom_line(aes(linetype=Family, color = Family))+
  geom_point(aes(shape=Family, color = Family)) +
  scale_color_manual(values = succession_top4_family_exp2) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
succession_exp2
ggsave("results/16S_FEPE_barplot_family_succession_level_exp2.pdf", width = 5, height = 4, dpi = 150) # save graphic

succession_exp3 <- ggplot(phy_prune_noncontam_prev05_true_exp3_filt_succession_top_fam_4_prune_agr, aes(x = succession_sample_order_exp3, y = Abundance, group = Family)) + #plotting by sample
  geom_line(aes(linetype=Family, color = Family))+
  geom_point(aes(shape=Family, color = Family)) +
  scale_color_manual(values = succession_top4_family_exp3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.position = "none")+
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
succession_exp3
ggsave("results/16S_FEPE_barplot_family_succession_level_exp3.pdf", width = 5, height = 4, dpi = 150) # save graphic

succession_figure <- ggarrange(succession_exp1, succession_exp2, succession_exp3,
                               nrow = 1, ncol = 3, common.legend = "none",
                               align = "h")
succession_figure
############################### Plot barcharts 16S samples per sample/experiment of top20 and top12 taxa ##############################################################
# Define list of colors and adjust accordingly
#nb.cols_no_controls_top20 <- 20
#colors_top20 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols_no_controls_top20)

# List of 20 distinct colors
colors_top20_phylum <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                         "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "gray")

# Plot bar chart 16S samples per sample/experiment - Phylum level - top 20 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3") # create new names for experiment variable
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3") # assign old name to new name
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr$SampleType, # reorder sample types
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Phylum list
phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr$Phylum <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr$Phylum,
                                                                  levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Bdellovibrionota", "Chloroflexi", "Cyanobacteria",
                                                                             "Deinococcota", "Dependentiae", "Desulfobacterota", "Fibrobacterota", "Firmicutes", "Marinimicrobia",
                                                                             "NB1-j", "Patescibacteria", "Planctomycetota", "Proteobacteria", "SAR324 MG-B", "Thermoplasmatota",
                                                                             "Verrucomicrobiota", "Others"))
# Plot
ggplot(phy_prune_noncontam_prev05_true_filt_top_phylum_top20_agr, aes(x = sample_order, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_phylum) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_phylum_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 8, height = 5, dpi = 150) # save graphic

# Plot bar chart 16S samples per sample/experiment - Phylum level  - top 12 taxa
colors_top12_phylum <- c("#46f0f0", "#1F78B4", "#aaffc3", "#33A02C", "#FB9A99", "#E31A1C", "#f032e6", "#FF7F00", "#ffe119", "#6A3D9A", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Phylum list
phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr$Phylum <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr$Phylum,
                                                                           levels = c("Actinobacteriota", "Bacteroidota", "Bdellovibrionota", "Cyanobacteria",
                                                                                      "Desulfobacterota", "Fibrobacterota", "Firmicutes", "Patescibacteria",
                                                                                      "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "Others"))
#Summarize by taxa
#phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr_summary <- phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr %>%
#  group_by(Experiment, SampleType, Phylum) %>%
#  summarise(
#    count = n(),
#    mean_phylum = mean(Abundance),
#    se_phylum = se(Abundance),
#    min_phylum = min(Abundance),
#   max_phylum = max(Abundance))

phylum <- ggplot(phy_prune_noncontam_prev05_true_filt_top_phylum_top12_agr, aes(x = sample_order, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_phylum) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title
  theme(plot.margin = unit(c(1,1,0.5,1), "lines")) +
  theme(panel.spacing.x = unit(0.5, "lines"))
phylum
ggsave("results/16S_FEPE_barplot_phylum_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic

# Plot top 5 phyla using facet
colors_top5_phylum <- c("#1F78B4", "#33A02C", "#ffe119", "#6A3D9A", "gray")


phylum.top5.labs <- c("Bacter.", "Cyanob.", "Planct.", "Proteo.", "Others")
names(phylum.top5.labs) <- c("Bacteroidota", "Cyanobacteria", "Planctomycetota", "Proteobacteria", "Others")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr$Phylum <- factor(phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr$Phylum,
                                                                           levels = c("Bacteroidota", "Cyanobacteria", "Planctomycetota",
                                                                                      "Proteobacteria", "Others"))
phylum_top5 <- ggplot(phy_prune_noncontam_prev05_true_filt_top_phylum_top5_agr, aes(x = sample_order, y = Abundance*100, fill = Phylum)) + #plotting by sample
  facet_grid(Phylum ~ Experiment, scales = "free",
             labeller = labeller(Experiment = experiment.labs,
                                 Phylum = phylum.top5.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top5_phylum) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  #scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  #theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  theme(legend.position = "none") +
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title
  theme(plot.margin = unit(c(1,1,0.5,1), "lines")) +
  theme(panel.spacing.x = unit(0.5, "lines"))
phylum_top5
ggsave("results/16S_FEPE_barplot_phylum_level_phy_prune_noncontam_prev05_true_filt_top5.pdf", width = 7, height = 6, dpi = 200) # save graphic
########################################
# Plot bar chart 16S samples per sample/experiment - Class level - top 20 taxa
colors_top20_class <- c("#9a6324", "#46f0f0", "#e6beff", "#f032e6", "#1F78B4",  "#4363d8", "#33A02C", "#FB9A99", "#e6194b","#6A3D9A", 
                        "#bcf60c", "#aaffc3", "#fabebe","#f58231", "#ffe119", "#008080", "#000075", "#fffac8", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top20_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Class list
phy_prune_noncontam_prev05_true_filt_top_class_top20_agr$Class <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top20_agr$Class,
                                                                         levels = c( "Acidimicrobiia", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Chlamydiae", "Cyanobacteriia",
                                                                                     "Desulfobacterota Class", "Fibrobacteria", "Gammaproteobacteria", "Marinimicrobia Class", "Oligoflexia", "Parcubacteria",
                                                                                     "Phycisphaerae", "Planctomycetes", "Planctomycetota OM190", "SAR324 MG-B Class", "Thermoplasmata", "Verrucomicrobiae", "Others"))

ggplot(phy_prune_noncontam_prev05_true_filt_top_class_top20_agr, aes(x = sample_order, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(Phylum ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_class) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "PF24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  #ggtitle(" 16S Phylum Levle") + # add the title on graphic
  #theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_class_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 8, height = 5, dpi = 150) # save graphic

# Plot bar chart 16S samples per sample/experiment - Class level - top 12 taxa
colors_top12_class <- c("#9a6324", "#46f0f0", "#e6beff", "#1F78B4",  "#33A02C", "#e6194b",
                        "#6A3D9A", "#aaffc3", "#f58231", "#ffe119", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top12_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Class list
phy_prune_noncontam_prev05_true_filt_top_class_top12_agr$Class <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top12_agr$Class,
                                                                         levels = c( "Acidimicrobiia", "Actinobacteria", "Alphaproteobacteria", "Bacteroidia", "Cyanobacteriia",
                                                                                     "Fibrobacteria", "Gammaproteobacteria", "Oligoflexia", "Phycisphaerae",
                                                                                     "Planctomycetes","Verrucomicrobiae", "Others"))

ggplot(phy_prune_noncontam_prev05_true_filt_top_class_top12_agr, aes(x = sample_order, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(Phylum ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_class) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "PF24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  #ggtitle(" 16S Phylum Levle") + # add the title on graphic
  #theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_class_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic

# Plot bar chart 16S samples per sample/experiment - Class level - top 8 taxa
colors_top8_class <- c("#e6beff", "#1F78B4",  "#33A02C", "#6A3D9A",
                       "#f58231", "#ffe119", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

phylum.top8.labs <- c("Bacter.", "Cyanob.", "Planct.", "Proteo.", "Verruc.", "Others")
names(phylum.top8.labs) <- c("Bacteroidota", "Cyanobacteria", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "Others")

# Put "Others" to the final of the Phylum list
phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$Phylum <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$Phylum,
                                                                          levels = c("Bacteroidota", "Cyanobacteria", "Planctomycetota",
                                                                                     "Proteobacteria", "Verrucomicrobiota", "Others"))

# Put "Others" to the final of the Class list
phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$Class <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top8_agr$Class,
                                                                         levels = c( "Alphaproteobacteria", "Bacteroidia", "Cyanobacteriia",
                                                                                     "Gammaproteobacteria", "Phycisphaerae",
                                                                                     "Planctomycetes","Verrucomicrobiae", "Others"))

ggplot(phy_prune_noncontam_prev05_true_filt_top_class_top8_agr, aes(x = sample_order, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(Phylum ~ Experiment, scales = "free",
             labeller = labeller(Experiment = experiment.labs, Phylum = phylum.top8.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top8_class) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent_format(accuracy = 1L), expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "PF24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  #ggtitle(" 16S Phylum Levle") + # add the title on graphic
  #theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_class_level_phy_prune_noncontam_prev05_true_filt_top8.pdf", width = 8, height = 4, dpi = 150) # save graphic


phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr

# Plot bar chart 16S samples per sample/experiment - Class level - top 8 taxa
colors_top6_class <- c("#e6beff", "#1F78B4",  "#33A02C", "#6A3D9A",
                       "#f58231", "#ffe119", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

phylum.top6.labs <- c("Bacter.", "Cyanob.", "Planct.", "Proteo.", "Verruc.", "Others")
names(phylum.top6.labs) <- c("Bacteroidota", "Cyanobacteria", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "Others")

# Put "Others" to the final of the Phylum list
phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr$Phylum <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top6_agr$Phylum,
                                                                         levels = c("Bacteroidota", "Cyanobacteria", "Planctomycetota",
                                                                                    "Proteobacteria", "Verrucomicrobiota", "Others"))

# Put "Others" to the final of the Class list
phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr$Class <- factor(phy_prune_noncontam_prev05_true_filt_top_class_top6_agr$Class,
                                                                        levels = c( "Alphaproteobacteria", "Bacteroidia", "Cyanobacteriia",
                                                                                    "Gammaproteobacteria", "Phycisphaerae",
                                                                                    "Planctomycetes","Verrucomicrobiae", "Others"))

ggplot(phy_prune_noncontam_prev05_true_filt_top_class_5_prune_agr, aes(x = sample_order, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(Phylum ~ Experiment,
             labeller = labeller(Experiment = experiment.labs, Phylum = phylum.top6.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top6_class) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent_format(accuracy = 1L), expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "PF24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  #ggtitle(" 16S Phylum Levle") + # add the title on graphic
  #theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_class_level_phy_prune_noncontam_prev05_true_filt_top6.pdf", width = 8, height = 4, dpi = 150) # save graphic


########################################
colors_top20_order <- c("#9a6324", "#46f0f0", "#3CB371", "#1F78B4", "#6A3D9A", "#e6194b", "#4363d8", "#FB9A99", "#f58231", "#ffe119",
                        "#FFFF99", "#BF40BF", "#aaffc3", "#e6beff", "#fabebe", "#000075", "#008080", "#33A02C", "#800000", "gray")

# Plot bar chart 16S samples per sample/experiment - Order Level  - top 20 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_order_top20_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Order list
phy_prune_noncontam_prev05_true_filt_top_order_top20_agr$Order <- factor(phy_prune_noncontam_prev05_true_filt_top_order_top20_agr$Order,
                                                                       levels = c( "Burkholderiales", "Coxiellales", "Cyanobacteriales", "Cytophagales", "Enterobacterales", "Fibrobacterales",
                                                                                   "Flavobacteriales", "Legionellales", "Phycisphaerales", "Pirellulales", "Planctomycetales", "Pseudomonadales",
                                                                                   "Rhizobiales", "Rhodobacterales", "Rhodospirillales", "SAR11", "Sneathiellales", "Synechococcales", "Verrucomicrobiales", "Others"))

ggplot(phy_prune_noncontam_prev05_true_filt_top_order_top20_agr, aes(x = sample_order, y = Abundance, fill = Order)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_order) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "PF24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_order_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 8, height = 5, dpi = 150) # save graphic

# Plot bar chart 16S samples per sample/experiment - Order Level - top 12 taxa
colors_top12_order <- c("#6A3D9A", "#4363d8", "#FB9A99", "#f58231", "#ffe119","#BF40BF",
                        "#aaffc3", "#e6beff", "#000075", "#008080", "#800000", "gray")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_order_top12_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Order list
phy_prune_noncontam_prev05_true_filt_top_order_top12_agr$Order <- factor(phy_prune_noncontam_prev05_true_filt_top_order_top12_agr$Order,
                                                                         levels = c( "Enterobacterales", "Fibrobacterales","Flavobacteriales", "Legionellales", "Phycisphaerales",
                                                                                     "Pirellulales", "Pseudomonadales", "Rhizobiales", "Rhodobacterales",
                                                                                     "SAR11", "Synechococcales", "Verrucomicrobiales", "Others"))

order <- ggplot(phy_prune_noncontam_prev05_true_filt_top_order_top12_agr, aes(x = sample_order, y = Abundance, fill = Order)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_order) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title
theme(strip.background = element_blank(), strip.text.x = element_blank()) +
theme(plot.margin = unit(c(1,1,0.5,1), "lines")) +
  theme(panel.spacing.x = unit(0.5, "lines"))
order
ggsave("results/16S_FEPE_barplot_order_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic

figure3 <- ggarrange(phylum, order, ncol =1, nrow =2,
                     labels = c("A", "B"),
                     legend = "right",
                     align = "v")
figure3
annotate_figure(figure3, left = text_grob("Relative Abundance (%)", face = "bold", size = 12, rot = 90))
ggsave("results/16S_FEPE_barplot_phylum_order_level_combined_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 9, dpi = 150) # save graphic

########################################
colors_top20_family <- c("#6A3D9A", "#aaffc3", "#33A02C", "#4363d8", "#FB9A99", "#1F78B4", "#e6194b", "#9a6324", "#f58231", "#ffe119",
                        "#FFFF99", "#BF40BF", "#e6beff", "#800000", "#000075", "#46f0f0", "#fabebe", "#3CB371", "#008080",  "gray")

# Plot bar chart 16S samples per sample/experiment - Family - top 20 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_family_top20_agr$SampleType,
                     level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))


# Put "Others" to the final of the Family list
phy_prune_noncontam_prev05_true_filt_top_family_top20_agr$Family <- factor(phy_prune_noncontam_prev05_true_filt_top_family_top20_agr$Family,
                                                                       levels = c("Alteromonadaceae", "Beijerinckiaceae", "Cyanobiaceae", "Flavobacteriaceae", "Legionellaceae", "Marinobacteraceae",
                                                                                  "Marinomonadaceae", "Nitrincolaceae", "Phycisphaeraceae", "Pirellulaceae", "Porticoccaceae", "Pseudoalteromonadaceae",
                                                                                  "Rhodobacteraceae", "Rubritaleaceae", "SAR11-Clade1", "SAR11-Clade2", "Spongiibacteraceae", "Stappiaceae", "Vibrionaceae", "Others")) 

family <- ggplot(phy_prune_noncontam_prev05_true_filt_top_family_top20_agr, aes(x = sample_order, y = Abundance, fill = Family)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_family) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  theme(axis.ticks = element_blank()) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title
  theme(plot.margin = unit(c(1,1,0.5,1), "lines")) +
  theme(panel.spacing.x = unit(0.5, "lines"))
family
ggsave("results/16S_FEPE_barplot_family_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 8, height = 5, dpi = 150) # save graphic
########################################
colors_top12_family <- c("#6A3D9A", "#33A02C", "#4363d8", "#e6194b", "#9a6324", "#ffe119",
                         "#FFFF99", "#BF40BF", "#e6beff", "#800000", "#000075", "gray")


# Plot bar chart 16S samples per sample/experiment - Family - top 12 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))


# Put "Others" to the final of the Family list
phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$Family <- factor(phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$Family,
                                                                           levels = c("Alteromonadaceae", "Cyanobiaceae", "Flavobacteriaceae", "Marinomonadaceae",
                                                                                      "Nitrincolaceae", "Pirellulaceae", "Porticoccaceae", "Pseudoalteromonadaceae",
                                                                                      "Rhodobacteraceae", "Rubritaleaceae", "SAR11-Clade1", "Others")) 

ggplot(phy_prune_noncontam_prev05_true_filt_top_family_top12_agr, aes(x = sample_order, y = Abundance, fill = Family)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_family) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_family_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))


# Put "Others" to the final of the Family list
phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$Family <- factor(phy_prune_noncontam_prev05_true_filt_top_family_top12_agr$Family,
                                                                           levels = c("Alteromonadaceae", "Cyanobiaceae", "Flavobacteriaceae", "Marinomonadaceae",
                                                                                      "Nitrincolaceae", "Pirellulaceae", "Porticoccaceae", "Pseudoalteromonadaceae",
                                                                                      "Rhodobacteraceae", "Rubritaleaceae", "SAR11-Clade1", "Others")) 

ggplot(phy_prune_noncontam_prev05_true_filt_eg_top_fam_12_prune_agr, aes(x = Family, y = Abundance, fill = Family)) + #plotting by sample
  facet_grid(Experiment ~ .,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_family) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_family_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic
########################################
colors_top20_genus <- c("#6A3D9A", "#ffe119", "#fabebe", "#33A02C", "#4363d8", "#1F78B4", "#FB9A99", "#e6194b", "#9a6324", "#aaffc3",
                        "#f58231", "#46f0f0", "#BF40BF", "#e6beff", "#FFFF99", "#800000", "#000075", "#008080", "#3CB371", "gray")

# Plot bar chart 16S samples per sample/experiment - Genus - top 20 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr$SampleType,
                     level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Family list
phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr$Genus <- factor(phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr$Genus,
                                                                         levels = c("Alteromonas", "Blastopirellula", "C1-B045",  "Cyanobium PCC-6307", "Flavobacteriaceae_gen", "Labrenzia", "Legionellaceae_gen", "Marinomonas",
                                                                                    "Neptuniibacter", "Pelagibaca", "Phycisphaeraceae CL500-3", "Pirellulaceae_gen", "Pseudoalteromonas", "Rhodobacteraceae_gen", "Rubripirellula",
                                                                                    "Rubritalea", "SAR11-Clade1a", "Shimia", "Synechococcus CC9902", "Others")) 

genus <- ggplot(phy_prune_noncontam_prev05_true_filt_top_genus_top20_agr, aes(x = sample_order, y = Abundance, fill = Genus)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_genus) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(1,1,0.5,1), "lines")) +
  theme(panel.spacing.x = unit(0.5, "lines"))
genus
ggsave("results/16S_FEPE_barplot_genus_level_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 8, height = 5, dpi = 150) # save graphic

########################################
colors_top12_genus <- c("#6A3D9A", "#ffe119", "#33A02C", "#e6194b",  "#46f0f0", "#BF40BF",
                        "#e6beff", "#FFFF99", "#800000", "#008080", "#3CB371", "gray")

# Plot bar chart 16S samples per sample/experiment - Genus - top 12 taxa
experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")
sample_order <- factor(phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr$SampleType,
                       level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"))

# Put "Others" to the final of the Family list
phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr$Genus <- factor(phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr$Genus,
                                                                         levels = c("Alteromonas", "Blastopirellula", "Cyanobium PCC-6307", "Marinomonas",
                                                                                    "Pirellulaceae_gen", "Pseudoalteromonas", "Rhodobacteraceae_gen", "Rubripirellula",
                                                                                    "Rubritalea", "Shimia", "Synechococcus CC9902", "Others")) 

ggplot(phy_prune_noncontam_prev05_true_filt_top_genus_top12_agr, aes(x = sample_order, y = Abundance, fill = Genus)) + #plotting by sample
  facet_grid(. ~ Experiment,
             labeller = labeller(Experiment = experiment.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12_genus) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    labels = c("SW", "FG", "FP2Hrs", "FP24Hrs",  "EG"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/16S_FEPE_barplot_genus_level_phy_prune_noncontam_prev05_true_filt_top12.pdf", width = 8, height = 4, dpi = 150) # save graphic

figure4 <- ggarrange(family, genus, ncol =1, nrow =2,
                             labels = c("A", "B"),
                             legend = "right",
                             align = "v")
figure4
annotate_figure(figure4, left = text_grob("Relative Abundance (%)", face = "bold", size = 12, rot = 90))

# save as pdf for additional edits using Inkspace
ggsave("results/16S_FEPE_barplot_family_genus_level_combined_phy_prune_noncontam_prev05_true_filt_top20.pdf", width = 10, height = 12, dpi = 150) # save graphic

##################################### Aldex2 on phyoseq objects #####################################

# Collapse tax table according to taxonomic rank for all three experiments
phy_prune_noncontam_prev05_true_filt_rank1 <- tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Phylum")
#phy_prune_noncontam_prev05_true_filt_rank2 <- tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Class")
#phy_prune_noncontam_prev05_true_filt_rank3 <- tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Order")
phy_prune_noncontam_prev05_true_filt_rank4 <- tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Family")
phy_prune_noncontam_prev05_true_filt_rank5 <- tax_glom(phy_prune_noncontam_prev05_true_filt, taxrank = "Genus")

# Collapse tax table according to taxonomic rank for experiment 1
phy_prune_noncontam_prev05_true_exp1_filt_rank1 <- tax_glom(phy_prune_noncontam_prev05_true_exp1_filt, taxrank = "Phylum")
get_taxa_unique(phy_prune_noncontam_prev05_true_exp1_filt_rank1, taxonomic.rank = "Phylum")
phy_prune_noncontam_prev05_true_exp1_filt_rank2 <- tax_glom(phy_prune_noncontam_prev05_true_exp1_filt, taxrank = "Class")
phy_prune_noncontam_prev05_true_exp1_filt_rank3 <- tax_glom(phy_prune_noncontam_prev05_true_exp1_filt, taxrank = "Order")
phy_prune_noncontam_prev05_true_exp1_filt_rank4 <- tax_glom(phy_prune_noncontam_prev05_true_exp1_filt, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp1_filt_rank5 <- tax_glom(phy_prune_noncontam_prev05_true_exp1_filt, taxrank = "Genus")

# Collapse tax table according to taxonomic rank for experiment 2
phy_prune_noncontam_prev05_true_exp2_filt_rank1 <- tax_glom(phy_prune_noncontam_prev05_true_exp2_filt, taxrank = "Phylum")
phy_prune_noncontam_prev05_true_exp2_filt_rank2 <- tax_glom(phy_prune_noncontam_prev05_true_exp2_filt, taxrank = "Class")
phy_prune_noncontam_prev05_true_exp2_filt_rank3 <- tax_glom(phy_prune_noncontam_prev05_true_exp2_filt, taxrank = "Order")
phy_prune_noncontam_prev05_true_exp2_filt_rank4 <- tax_glom(phy_prune_noncontam_prev05_true_exp2_filt, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp2_filt_rank5 <- tax_glom(phy_prune_noncontam_prev05_true_exp2_filt, taxrank = "Genus")

# Collapse tax table according to taxonomic rank for experiment 3
phy_prune_noncontam_prev05_true_exp3_filt_rank1 <- tax_glom(phy_prune_noncontam_prev05_true_exp3_filt, taxrank = "Phylum")
#phy_prune_noncontam_prev05_true_exp3_filt_rank2 <- tax_glom(phy_prune_noncontam_prev05_true_exp3_filt, taxrank = "Class")
#phy_prune_noncontam_prev05_true_exp3_filt_rank3 <- tax_glom(phy_prune_noncontam_prev05_true_exp3_filt, taxrank = "Order")
phy_prune_noncontam_prev05_true_exp3_filt_rank4 <- tax_glom(phy_prune_noncontam_prev05_true_exp3_filt, taxrank = "Family")
phy_prune_noncontam_prev05_true_exp3_filt_rank5 <- tax_glom(phy_prune_noncontam_prev05_true_exp3_filt, taxrank = "Genus")

# Apply Aldex2 on phyoseq objects for the different ranks all three experiments
aldex2_phy_group <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank1)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_filt_rank1)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_rank1_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank1))

aldex2_fam_group <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank4)),
                                  phyloseq::sample_data(phy_prune_noncontam_prev05_true_filt_rank1)$SampleType,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_rank4_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank4))

aldex2_genus_group <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank5)),
                                  phyloseq::sample_data(phy_prune_noncontam_prev05_true_filt_rank1)$SampleType,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_rank5_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_filt_rank5))

# Apply Aldex2 on phyoseq objects for the different ranks of experiment 1
aldex2_phy_group_exp1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank1)),
                                phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt_rank1)$SampleType,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp1_rank1_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank1))
group_exp1_rank1_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp1_filt_rank1))
taxa_info_rank1_exp1 <- group_exp1_rank1_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank1_exp1 <- group_exp1_rank1_otu_table %>% rownames_to_column(var = "OTU")

aldex2_class_group_exp1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank2)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt_rank2)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp1_rank2_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank2))
group_exp1_rank2_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp1_filt_rank2))
taxa_info_rank2_exp1 <- group_exp1_rank2_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank2_exp1 <- group_exp1_rank2_otu_table %>% rownames_to_column(var = "OTU")

aldex2_order_group_exp1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank3)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt_rank3)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp1_rank3_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank3))
group_exp1_rank3_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp1_filt_rank3))
taxa_info_rank3_exp1 <- group_exp1_rank3_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank3_exp1 <- group_exp1_rank3_otu_table %>% rownames_to_column(var = "OTU")

aldex2_fam_group_exp1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank4)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt_rank4)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp1_rank4_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank4))
group_exp1_rank4_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp1_filt_rank4))
taxa_info_rank4_exp1 <- group_exp1_rank4_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank4_exp1 <- group_exp1_rank4_otu_table %>% rownames_to_column(var = "OTU")

aldex2_genus_group_exp1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank5)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp1_filt_rank5)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp1_rank5_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp1_filt_rank5))
group_exp1_rank5_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp1_filt_rank5))
taxa_info_rank5_exp1 <- group_exp1_rank5_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank5_exp1 <- group_exp1_rank5_otu_table %>% rownames_to_column(var = "OTU")


# Apply Aldex2 on phyoseq objects for the different ranks of experiment 2
aldex2_phy_group_exp2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank1)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt_rank1)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp2_rank1_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank1))
group_exp2_rank1_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp2_filt_rank1))
taxa_info_rank1_exp2 <- group_exp2_rank1_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank1_exp2 <- group_exp2_rank1_otu_table %>% rownames_to_column(var = "OTU")


aldex2_class_group_exp2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank2)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt_rank2)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp2_rank2_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank2))
group_exp2_rank2_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp2_filt_rank2))
taxa_info_rank2_exp2 <- group_exp2_rank2_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank2_exp2 <- group_exp2_rank2_otu_table %>% rownames_to_column(var = "OTU")


aldex2_order_group_exp2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank3)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt_rank3)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp2_rank3_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank3))
group_exp2_rank3_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp2_filt_rank3))
taxa_info_rank3_exp2 <- group_exp2_rank3_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank3_exp2 <- group_exp2_rank3_otu_table %>% rownames_to_column(var = "OTU")

aldex2_fam_group_exp2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank4)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt_rank4)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp2_rank4_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank4))
group_exp2_rank4_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp2_filt_rank4))
taxa_info_rank4_exp2 <- group_exp2_rank4_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank4_exp2 <- group_exp2_rank4_otu_table %>% rownames_to_column(var = "OTU")

aldex2_genus_group_exp2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank5)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp2_filt_rank5)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp2_rank5_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp2_filt_rank5))
group_exp2_rank5_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp2_filt_rank5))
taxa_info_rank5_exp2 <- group_exp2_rank5_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank5_exp2 <- group_exp2_rank5_otu_table %>% rownames_to_column(var = "OTU")

# Apply Aldex2 on phyoseq objects for the different ranks of experiment 3
aldex2_phy_group_exp3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank1)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt_rank1)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp3_rank1_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank1))
group_exp3_rank1_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp3_filt_rank1))
taxa_info_rank1_exp3 <- group_exp3_rank1_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank1_exp3 <- group_exp3_rank1_otu_table %>% rownames_to_column(var = "OTU")

aldex2_class_group_exp3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank2)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt_rank2)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp3_rank2_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank2))
group_exp3_rank2_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp3_filt_rank2))
taxa_info_rank2_exp3 <- group_exp3_rank2_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank2_exp3 <- group_exp3_rank2_otu_table %>% rownames_to_column(var = "OTU")

aldex2_order_group_exp3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank3)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt_rank3)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp3_rank3_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank3))
group_exp3_rank3_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp3_filt_rank3))
taxa_info_rank3_exp3 <- group_exp3_rank3_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank3_exp3 <- group_exp3_rank3_otu_table %>% rownames_to_column(var = "OTU")

aldex2_fam_group_exp3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank4)),
                                       phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt_rank4)$SampleType,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp3_rank4_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank4))
group_exp3_rank4_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp3_filt_rank4))
taxa_info_rank4_exp3 <- group_exp3_rank4_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank4_exp3 <- group_exp3_rank4_otu_table %>% rownames_to_column(var = "OTU")

aldex2_genus_group_exp3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank5)),
                                         phyloseq::sample_data(phy_prune_noncontam_prev05_true_exp3_filt_rank5)$SampleType,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
group_exp3_rank5_otu_table <- data.frame(phyloseq::otu_table(phy_prune_noncontam_prev05_true_exp3_filt_rank5))
group_exp3_rank5_tax_table <- data.frame(phyloseq::tax_table(phy_prune_noncontam_prev05_true_exp3_filt_rank5))
taxa_info_rank5_exp3 <- group_exp3_rank5_tax_table %>% rownames_to_column(var = "OTU")
otu_info_rank5_exp3 <- group_exp3_rank5_otu_table %>% rownames_to_column(var = "OTU")

# Write Aldex2 results to file as csv, all three experiments
aldex2_phy_group <- rownames_to_column(aldex2_phy_group, var ="OTU")
#aldex2_class_group <- rownames_to_column(aldex2_class_group, var ="OTU")
#aldex2_order_group <- rownames_to_column(aldex2_order_group, var ="OTU")
aldex2_fam_group <- rownames_to_column(aldex2_fam_group, var ="OTU")
aldex2_genus_group <- rownames_to_column(aldex2_genus_group, var ="OTU")

write.csv(aldex2_phy_group, "aldex/aldex2_phy_group.csv")
#write.csv(aldex2_class_group, "aldex/aldex2_class_group.csv")
#write.csv(aldex2_order_group, "aldex/aldex2_order_group.csv")
write.csv(aldex2_fam_group, "aldex/aldex2_fam_group.csv")
write.csv(aldex2_genus_group, "aldex/aldex2_genus_group.csv")

# Write Aldex2 results to file as csv, experiment 1
aldex2_phy_group_exp1 <- rownames_to_column(aldex2_phy_group_exp1, var ="OTU")
aldex2_class_group_exp1 <- rownames_to_column(aldex2_class_group_exp1, var ="OTU")
aldex2_order_group_exp1 <- rownames_to_column(aldex2_order_group_exp1, var ="OTU")
aldex2_fam_group_exp1 <- rownames_to_column(aldex2_fam_group_exp1, var ="OTU")
aldex2_genus_group_exp1 <- rownames_to_column(aldex2_genus_group_exp1, var ="OTU")

write.csv(aldex2_phy_group_exp1, "aldex/aldex2_phy_group_exp1.csv")
write.csv(aldex2_class_group_exp1, "aldex/aldex2_class_group_exp1.csv")
write.csv(aldex2_order_group_exp1, "aldex/aldex2_order_group_exp1.csv")
write.csv(aldex2_fam_group_exp1, "aldex/aldex2_fam_group_exp1.csv")
write.csv(aldex2_genus_group_exp1, "aldex/aldex2_genus_group_exp1.csv")

# Write Aldex2 results to file as csv, experiment 2
aldex2_phy_group_exp2 <- rownames_to_column(aldex2_phy_group_exp2, var ="OTU")
aldex2_class_group_exp2 <- rownames_to_column(aldex2_class_group_exp2, var ="OTU")
aldex2_order_group_exp2 <- rownames_to_column(aldex2_order_group_exp2, var ="OTU")
aldex2_fam_group_exp2 <- rownames_to_column(aldex2_fam_group_exp2, var ="OTU")
aldex2_genus_group_exp2 <- rownames_to_column(aldex2_genus_group_exp2, var ="OTU")

write.csv(aldex2_phy_group_exp2, "aldex/aldex2_phy_group_exp2.csv")
write.csv(aldex2_class_group_exp2, "aldex/aldex2_class_group_exp2.csv")
write.csv(aldex2_order_group_exp2, "aldex/aldex2_order_group_exp2.csv")
write.csv(aldex2_fam_group_exp2, "aldex/aldex2_fam_group_exp2.csv")
write.csv(aldex2_genus_group_exp2, "aldex/aldex2_genus_group_exp2.csv")

# Write Aldex2 results to file as csv, experiment 3
aldex2_phy_group_exp3 <- rownames_to_column(aldex2_phy_group_exp3, var ="OTU")
aldex2_class_group_exp3 <- rownames_to_column(aldex2_class_group_exp3, var ="OTU")
aldex2_order_group_exp3 <- rownames_to_column(aldex2_order_group_exp3, var ="OTU")
aldex2_fam_group_exp3 <- rownames_to_column(aldex2_fam_group_exp3, var ="OTU")
aldex2_genus_group_exp3 <- rownames_to_column(aldex2_genus_group_exp3, var ="OTU")

write.csv(aldex2_phy_group_exp3, "aldex/aldex2_phy_group_exp3.csv")
write.csv(aldex2_class_group_exp3, "aldex/aldex2_class_group_exp3.csv")
write.csv(aldex2_order_group_exp3, "aldex/aldex2_order_group_exp3.csv")
write.csv(aldex2_fam_group_exp3, "aldex/aldex2_fam_group_exp3.csv")
write.csv(aldex2_genus_group_exp3, "aldex/aldex2_genus_group_exp3.csv")

# Import Aldex2 results and rename the X variable by OTU, all three experiments
aldex2_phy_group_result <- read.csv("aldex/aldex2_phy_group.csv")
aldex2_phy_group_result <- aldex2_phy_group_result[, -1]
#aldex2_class_group_result <- read.csv("aldex/aldex2_class_group.csv")
#aldex2_class_group_result <- aldex2_class_group_result[, -1]
#aldex2_order_group_result <- read.csv("aldex/aldex2_order_group.csv")
#aldex2_order_group_result <- aldex2_order_group_result[, -1]
aldex2_fam_group_result <- read.csv("aldex/aldex2_fam_group.csv")
aldex2_fam_group_result <- aldex2_fam_group_result[, -1]
aldex2_genus_group_result <- read.csv("aldex/aldex2_genus_group.csv")
aldex2_genus_group_result <- aldex2_genus_group_result[, -1]

# Import Aldex2 results and rename the X variable by OTU, experiment 1
aldex2_phy_group_exp1_result <- as.data.frame(read.csv("aldex/aldex2_phy_group_exp1.csv", sep=",", header=T))
aldex2_class_group_exp1_result <- as.data.frame(read.csv("aldex/aldex2_class_group_exp1.csv", sep=",", header=T))
aldex2_order_group_exp1_result <- as.data.frame(read.csv("aldex/aldex2_order_group_exp1.csv", sep=",", header=T))
aldex2_fam_group_exp1_result <- as.data.frame(read.csv("aldex/aldex2_fam_group_exp1.csv", sep=",", header=T))
aldex2_genus_group_exp1_result <- as.data.frame(read.csv("aldex/aldex2_genus_group_exp1.csv", sep=",", header=T))

# Import Aldex2 results and rename the X variable by OTU, experiment 2
aldex2_phy_group_exp2_result <- as.data.frame(read.csv("aldex/aldex2_phy_group_exp2.csv", sep=",", header=T))
aldex2_class_group_exp2_result <- as.data.frame(read.csv("aldex/aldex2_class_group_exp2.csv", sep=",", header=T))
aldex2_order_group_exp2_result <- as.data.frame(read.csv("aldex/aldex2_order_group_exp2.csv", sep=",", header=T))
aldex2_fam_group_exp2_result <- as.data.frame(read.csv("aldex/aldex2_fam_group_exp2.csv", sep=",", header=T))
aldex2_genus_group_exp2_result <- as.data.frame(read.csv("aldex/aldex2_genus_group_exp2.csv", sep=",", header=T))

# Import Aldex2 results and rename the X variable by OTU, experiment 3
aldex2_phy_group_exp3_result <- as.data.frame(read.csv("aldex/aldex2_phy_group_exp3.csv", sep=",", header=T))
aldex2_class_group_exp3_result <- as.data.frame(read.csv("aldex/aldex2_class_group_exp3.csv", sep=",", header=T))
aldex2_order_group_exp3_result <- as.data.frame(read.csv("aldex/aldex2_order_group_exp3.csv", sep=",", header=T))
aldex2_fam_group_exp3_result <- as.data.frame(read.csv("aldex/aldex2_fam_group_exp3.csv", sep=",", header=T))
aldex2_genus_group_exp3_result <- as.data.frame(read.csv("aldex/aldex2_genus_group_exp3.csv", sep=",", header=T))

#Clean up presentation, all three experiments
taxa_info <- data.frame(tax_table(phy_prune_noncontam_prev05_true_filt))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

otu_count <- data.frame(otu_table(phy_prune_noncontam_prev05_true_filt))
otu_count <- otu_count %>% rownames_to_column(var = "OTU")

sample_tab_group <- data.frame(sample_data(phy_prune_noncontam_prev05_true_filt))
sample_tab_group_factors <- sample_tab_group[, c(1:3)]

#Clean up presentation, experiment 1
sample_tab_group_exp1 <- data.frame(sample_data(phy_prune_noncontam_prev05_true_exp1_filt))
sample_tab_group_exp1_factors <- sample_tab_group_exp1[, c(1:3)]

#Clean up presentation, experiment 2
sample_tab_group_exp2 <- data.frame(sample_data(phy_prune_noncontam_prev05_true_exp2_filt))
sample_tab_group_exp2_factors <- sample_tab_group_exp2[, c(1:3)]

#Clean up presentation, experiment 3
sample_tab_group_exp3 <- data.frame(sample_data(phy_prune_noncontam_prev05_true_exp3_filt))
sample_tab_group_exp3_factors <- sample_tab_group_exp3[, c(1:3)]

# Filter Aldex2 results by sig kw.eBH and join the taxonomic information, all three experiments
sig_aldex2_phy_group_result <- aldex2_phy_group_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_group_result <- left_join(sig_aldex2_phy_group_result, taxa_info)
write.csv(sig_aldex2_phy_group_result, "aldex/sig_aldex2_phy_group_result.csv")

#sig_aldex2_class_group_result <- aldex2_class_group_result %>%
#  filter(kw.ep < 0.05) %>%
#  arrange(kw.ep, kw.eBH) %>%
#  dplyr::select(OTU, kw.ep, kw.eBH)
#sig_aldex2_class_group_result <- left_join(sig_aldex2_class_group_result, taxa_info)
#write.csv(sig_aldex2_class_group_result, "aldex/sig_aldex2_class_group_result.csv")

#sig_aldex2_order_group_result <- aldex2_order_group_result %>%
#  filter(kw.ep < 0.05) %>%
#  arrange(kw.ep, kw.eBH) %>%
#  dplyr::select(OTU, kw.ep, kw.eBH)
#sig_aldex2_order_group_result <- left_join(sig_aldex2_order_group_result, taxa_info)
#write.csv(sig_aldex2_order_group_result, "aldex/sig_aldex2_order_group_result.csv")

sig_aldex2_fam_group_result <- aldex2_fam_group_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_group_result <- left_join(sig_aldex2_fam_group_result, taxa_info)
write.csv(sig_aldex2_fam_group_result, "aldex/sig_aldex2_fam_group_result.csv")
sig_aldex2_fam_group_result_top30 <- sig_aldex2_fam_group_result %>% top_n(-30, kw.ep)
write.csv(sig_aldex2_fam_group_result_top30, "aldex/sig_aldex2_fam_group_result_top30.csv")

sig_aldex2_genus_group_result <- aldex2_genus_group_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_group_result <- left_join(sig_aldex2_genus_group_result, taxa_info)
write.csv(sig_aldex2_genus_group_result, "aldex/sig_aldex2_genus_group_result.csv")
sig_aldex2_genus_group_result_top30 <- sig_aldex2_genus_group_result %>% top_n(-30, kw.ep)
write.csv(sig_aldex2_genus_group_result_top30, "aldex/sig_aldex2_genus_group_result_top30.csv")

# Filter Aldex2 results by sig kw.eBH and join the taxonomic information, experiment 1
sig_aldex2_phy_group_exp1_result <- aldex2_phy_group_exp1_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_group_exp1_result <- left_join(sig_aldex2_phy_group_exp1_result, taxa_info_rank1_exp1)
sig_aldex2_phy_group_exp1_result_top10 <- sig_aldex2_phy_group_exp1_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_phy_group_exp1_result, "aldex/sig_aldex2_phy_group_exp1_result.csv")

sig_aldex2_class_group_exp1_result <- aldex2_class_group_exp1_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_group_exp1_result <- left_join(sig_aldex2_class_group_exp1_result, taxa_info_rank2_exp1)
write.csv(sig_aldex2_class_group_exp1_result, "aldex/sig_aldex2_class_group_exp1_result.csv")

sig_aldex2_order_group_exp1_result <- aldex2_order_group_exp1_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_group_exp1_result <- left_join(sig_aldex2_order_group_exp1_result, taxa_info_rank3_exp1)
sig_aldex2_order_group_exp1_result_top20 <- sig_aldex2_order_group_exp1_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_order_group_exp1_result, "aldex/sig_aldex2_order_group_exp1_result.csv")

sig_aldex2_fam_group_exp1_result <- aldex2_fam_group_exp1_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_group_exp1_result <- left_join(sig_aldex2_fam_group_exp1_result, taxa_info_rank4_exp1)
sig_aldex2_fam_group_exp1_result_top20 <- sig_aldex2_fam_group_exp1_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_fam_group_exp1_result, "aldex/sig_aldex2_fam_group_exp1_result.csv")

sig_aldex2_genus_group_exp1_result <- aldex2_genus_group_exp1_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_group_exp1_result <- left_join(sig_aldex2_genus_group_exp1_result, taxa_info_rank5_exp1)
sig_aldex2_genus_group_exp1_result_top20 <- sig_aldex2_genus_group_exp1_result %>% top_n(-20, kw.ep)
sig_aldex2_genus_group_exp1_result_top10 <- sig_aldex2_genus_group_exp1_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_genus_group_exp1_result, "aldex/sig_aldex2_genus_group_exp1_result.csv")

# Filter Aldex2 results by sig kw.eBH and join the taxonomic information, experiment 2
sig_aldex2_phy_group_exp2_result <- aldex2_phy_group_exp2_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_group_exp2_result <- left_join(sig_aldex2_phy_group_exp2_result, taxa_info_rank1_exp2)
sig_aldex2_phy_group_exp2_result_top10 <- sig_aldex2_phy_group_exp2_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_phy_group_exp2_result, "aldex/sig_aldex2_phy_group_exp2_result.csv")

sig_aldex2_class_group_exp2_result <- aldex2_class_group_exp2_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_group_exp2_result <- left_join(sig_aldex2_class_group_exp2_result, taxa_info_rank2_exp2)
sig_aldex2_class_group_exp2_result_top20 <- sig_aldex2_class_group_exp2_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_class_group_exp2_result, "aldex/sig_aldex2_class_group_exp2_result.csv")

sig_aldex2_order_group_exp2_result <- aldex2_order_group_exp2_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_group_exp2_result <- left_join(sig_aldex2_order_group_exp2_result, taxa_info_rank3_exp2)
sig_aldex2_order_group_exp2_result_top20 <- sig_aldex2_order_group_exp2_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_order_group_exp2_result, "aldex/sig_aldex2_order_group_exp2_result.csv")

sig_aldex2_fam_group_exp2_result <- aldex2_fam_group_exp2_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_group_exp2_result <- left_join(sig_aldex2_fam_group_exp2_result, taxa_info_rank4_exp2)
sig_aldex2_fam_group_exp2_result_top20 <- sig_aldex2_fam_group_exp2_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_fam_group_exp2_result, "aldex/sig_aldex2_fam_group_exp2_result.csv")

sig_aldex2_genus_group_exp2_result <- aldex2_genus_group_exp2_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_group_exp2_result <- left_join(sig_aldex2_genus_group_exp2_result, taxa_info_rank5_exp2)
sig_aldex2_genus_group_exp2_result_top20 <- sig_aldex2_genus_group_exp2_result %>% top_n(-20, kw.ep)
sig_aldex2_genus_group_exp2_result_top10 <- sig_aldex2_genus_group_exp2_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_genus_group_exp2_result, "aldex/sig_aldex2_genus_group_exp2_result.csv")

# Filter Aldex2 results by sig kw.eBH and join the taxonomic information, experiment 3
sig_aldex2_phy_group_exp3_result <- aldex2_phy_group_exp3_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_group_exp3_result <- left_join(sig_aldex2_phy_group_exp3_result, taxa_info_rank1_exp3)
sig_aldex2_phy_group_exp3_result_top10 <- sig_aldex2_phy_group_exp3_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_phy_group_exp3_result, "aldex/sig_aldex2_phy_group_exp3_result.csv")

sig_aldex2_class_group_exp3_result <- aldex2_class_group_exp3_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_group_exp3_result <- left_join(sig_aldex2_class_group_exp3_result, taxa_info_rank2_exp3)
sig_aldex2_class_group_exp3_result_top20 <- sig_aldex2_class_group_exp3_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_class_group_exp3_result, "aldex/sig_aldex2_class_group_exp3_result.csv")

sig_aldex2_order_group_exp3_result <- aldex2_order_group_exp3_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_group_exp3_result <- left_join(sig_aldex2_order_group_exp3_result, taxa_info_rank3_exp3)
sig_aldex2_order_group_exp3_result_top20 <- sig_aldex2_order_group_exp3_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_order_group_exp3_result, "aldex/sig_aldex2_order_group_exp3_result.csv")

sig_aldex2_fam_group_exp3_result <- aldex2_fam_group_exp3_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_group_exp3_result <- left_join(sig_aldex2_fam_group_exp3_result, taxa_info_rank4_exp3)
sig_aldex2_fam_group_exp3_result_top20 <- sig_aldex2_fam_group_exp3_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_fam_group_exp3_result, "aldex/sig_aldex2_fam_group_exp3_result.csv")

sig_aldex2_genus_group_exp3_result <- aldex2_genus_group_exp3_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_group_exp3_result <- left_join(sig_aldex2_genus_group_exp3_result, taxa_info_rank5_exp3)
sig_aldex2_genus_group_exp3_result_top20 <- sig_aldex2_genus_group_exp3_result %>% top_n(-20, kw.ep)
sig_aldex2_genus_group_exp3_result_top10 <- sig_aldex2_genus_group_exp3_result %>% top_n(-10, kw.ep)
write.csv(sig_aldex2_genus_group_exp3_result, "aldex/sig_aldex2_genus_group_exp3_result.csv")

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20, all three experiments
sig_aldex2_phy_group_result_count <- left_join(sig_aldex2_phy_group_result, otu_count)
write.csv(sig_aldex2_phy_group_result_count, "aldex/sig_aldex2_phy_group_result_count.csv")
#sig_aldex2_class_group_result_count <- left_join(sig_aldex2_class_group_result, otu_count)
#write.csv(sig_aldex2_class_group_result_count, "aldex/sig_aldex2_class_group_result_count.csv")
#sig_aldex2_order_group_result_count <- left_join(sig_aldex2_order_group_result, otu_count)
#write.csv(sig_aldex2_order_group_result_count, "aldex/sig_aldex2_order_group_result_count.csv")
sig_aldex2_fam_group_result_count <- left_join(sig_aldex2_fam_group_result, otu_count)
write.csv(sig_aldex2_fam_group_result_count, "aldex/sig_aldex2_fam_group_result_count.csv")
sig_aldex2_fam_group_result_top30_count <- left_join(sig_aldex2_fam_group_result_top30, otu_count)
write.csv(sig_aldex2_fam_group_result_top30_count, "aldex/sig_aldex2_fam_group_result_top30_count.csv")
sig_aldex2_genus_group_result_count <- left_join(sig_aldex2_genus_group_result, otu_count)
write.csv(sig_aldex2_genus_group_result_count, "aldex/sig_aldex2_genus_group_result_count.csv")
sig_aldex2_genus_group_result_top30_count <- left_join(sig_aldex2_genus_group_result_top30, otu_count)
write.csv(sig_aldex2_genus_group_result_top30_count, "aldex/sig_aldex2_genus_group_result_top30_count.csv")

clr_phy_group <- sig_aldex2_phy_group_result_count[, -(2:10)]
rownames(clr_phy_group) <- clr_phy_group$OTU
clr_phy_group <- clr_phy_group[, -1]
#clr_class_group <- sig_aldex2_class_group_result_count[, -(2:10)]
#rownames(clr_class_group) <- clr_class_group$OTU
#clr_class_group <- clr_class_group[, -1]
#clr_order_group <- sig_aldex2_order_group_result_count[, -(2:10)]
#rownames(clr_order_group) <- clr_order_group$OTU
#clr_order_group <- clr_order_group[, -1]
clr_fam_group <- sig_aldex2_fam_group_result_count[, -(2:10)]
rownames(clr_fam_group) <- clr_fam_group$OTU
clr_fam_group <- clr_fam_group[, -1]
clr_fam_group_top30 <- sig_aldex2_fam_group_result_top30_count[, -(2:10)]
rownames(clr_fam_group_top30) <- clr_fam_group_top30$OTU
clr_fam_group_top30 <- clr_fam_group_top30[, -1]
clr_genus_group <- sig_aldex2_genus_group_result_count[, -(2:10)]
rownames(clr_genus_group) <- clr_genus_group$OTU
clr_genus_group <- clr_genus_group[, -1]
clr_genus_group_top30 <- sig_aldex2_genus_group_result_top30_count[, -(2:10)]
rownames(clr_genus_group_top30) <- clr_genus_group_top30$OTU
clr_genus_group_top30 <- clr_genus_group_top30[, -1]

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20, experiment 1
sig_aldex2_phy_group_exp1_result_count <- left_join(sig_aldex2_phy_group_exp1_result, otu_info_rank1_exp1)
write.csv(sig_aldex2_phy_group_exp1_result_count, "aldex/sig_aldex2_phy_group_exp1_result_count.csv")
clr_phy_group_exp1 <- sig_aldex2_phy_group_exp1_result_count[, -(2:10)]
rownames(clr_phy_group_exp1) <- clr_phy_group_exp1$OTU
clr_phy_group_exp1 <- clr_phy_group_exp1[, -1]

sig_aldex2_phy_group_exp1_result_top10_count <- left_join(sig_aldex2_phy_group_exp1_result_top10, otu_info_rank1_exp1)
clr_phy_top10_group_exp1 <- sig_aldex2_phy_group_exp1_result_top10_count[, -(2:10)]
rownames(clr_phy_top10_group_exp1) <- clr_phy_top10_group_exp1$OTU
clr_phy_top10_group_exp1 <- clr_phy_top10_group_exp1[, -1]


sig_aldex2_class_group_exp1_result_count <- left_join(sig_aldex2_class_group_exp1_result, otu_info_rank2_exp1)
write.csv(sig_aldex2_class_group_exp1_result_count, "aldex/sig_aldex2_class_group_exp1_result_count.csv")
clr_class_group_exp1 <- sig_aldex2_class_group_exp1_result_count[, -(2:10)]
rownames(clr_class_group_exp1) <- clr_class_group_exp1$OTU
clr_class_group_exp1 <- clr_class_group_exp1[, -1]

sig_aldex2_order_group_exp1_result_count <- left_join(sig_aldex2_order_group_exp1_result, otu_info_rank3_exp1)
write.csv(sig_aldex2_order_group_exp1_result_count, "aldex/sig_aldex2_order_group_exp1_result_count.csv")
clr_order_group_exp1 <- sig_aldex2_order_group_exp1_result_count[, -(2:10)]
rownames(clr_order_group_exp1) <- clr_order_group_exp1$OTU
clr_order_group_exp1 <- clr_order_group_exp1[, -1]

sig_aldex2_order_group_exp1_result_top20_count <- left_join(sig_aldex2_order_group_exp1_result_top20, otu_info_rank3_exp1)
clr_order_top20_group_exp1 <- sig_aldex2_order_group_exp1_result_top20_count[, -(2:10)]
rownames(clr_order_top20_group_exp1) <- clr_order_top20_group_exp1$OTU
clr_order_top20_group_exp1 <- clr_order_top20_group_exp1[, -1]

sig_aldex2_fam_group_exp1_result_count <- left_join(sig_aldex2_fam_group_exp1_result, otu_info_rank4_exp1)
write.csv(sig_aldex2_fam_group_exp1_result_count, "aldex/sig_aldex2_fam_group_exp1_result_count.csv")
clr_fam_group_exp1 <- sig_aldex2_fam_group_exp1_result_count[, -(2:10)]
rownames(clr_fam_group_exp1) <- clr_fam_group_exp1$OTU
clr_fam_group_exp1 <- clr_fam_group_exp1[, -1]

sig_aldex2_fam_group_exp1_result_top20_count <- left_join(sig_aldex2_fam_group_exp1_result_top20, otu_info_rank4_exp1)
clr_fam_top20_group_exp1 <- sig_aldex2_fam_group_exp1_result_top20_count[, -(2:10)]
rownames(clr_fam_top20_group_exp1) <- clr_fam_top20_group_exp1$OTU
clr_fam_top20_group_exp1 <- clr_fam_top20_group_exp1[, -1]

sig_aldex2_genus_group_exp1_result_count <- left_join(sig_aldex2_genus_group_exp1_result, otu_info_rank5_exp1)
write.csv(sig_aldex2_genus_group_exp1_result_count, "aldex/sig_aldex2_genus_group_exp1_result_count.csv")
clr_genus_group_exp1 <- sig_aldex2_genus_group_exp1_result_count[, -(2:10)]
rownames(clr_genus_group_exp1) <- clr_genus_group_exp1$OTU
clr_genus_group_exp1 <- clr_genus_group_exp1[, -1]

sig_aldex2_genus_group_exp1_result_top20_count <- left_join(sig_aldex2_genus_group_exp1_result_top20, otu_info_rank5_exp1)
clr_genus_top20_group_exp1 <- sig_aldex2_genus_group_exp1_result_top20_count[, -(2:10)]
rownames(clr_genus_top20_group_exp1) <- clr_genus_top20_group_exp1$OTU
clr_genus_top20_group_exp1 <- clr_genus_top20_group_exp1[, -1]

sig_aldex2_genus_group_exp1_result_top10_count <- left_join(sig_aldex2_genus_group_exp1_result_top10, otu_info_rank5_exp1)
clr_genus_top10_group_exp1 <- sig_aldex2_genus_group_exp1_result_top10_count[, -(2:10)]
rownames(clr_genus_top10_group_exp1) <- clr_genus_top10_group_exp1$OTU
clr_genus_top10_group_exp1 <- clr_genus_top10_group_exp1[, -1]

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20, experiment 2
sig_aldex2_phy_group_exp2_result_count <- left_join(sig_aldex2_phy_group_exp2_result, otu_info_rank1_exp2)
write.csv(sig_aldex2_phy_group_exp2_result_count, "aldex/sig_aldex2_phy_group_exp2_result_count.csv")
clr_phy_group_exp2 <- sig_aldex2_phy_group_exp2_result_count[, -(2:10)]
rownames(clr_phy_group_exp2) <- clr_phy_group_exp2$OTU
clr_phy_group_exp2 <- clr_phy_group_exp2[, -1]

sig_aldex2_phy_group_exp2_result_top10_count <- left_join(sig_aldex2_phy_group_exp2_result_top10, otu_info_rank1_exp2)
clr_phy_top10_group_exp2 <- sig_aldex2_phy_group_exp2_result_top10_count[, -(2:10)]
rownames(clr_phy_top10_group_exp2) <- clr_phy_top10_group_exp2$OTU
clr_phy_top10_group_exp2 <- clr_phy_top10_group_exp2[, -1]

sig_aldex2_class_group_exp2_result_count <- left_join(sig_aldex2_class_group_exp2_result, otu_info_rank2_exp2)
write.csv(sig_aldex2_class_group_exp2_result_count, "aldex/sig_aldex2_class_group_exp2_result_count.csv")
clr_class_group_exp2 <- sig_aldex2_class_group_exp2_result_count[, -(2:10)]
rownames(clr_class_group_exp2) <- clr_class_group_exp2$OTU
clr_class_group_exp2 <- clr_class_group_exp2[, -1]

sig_aldex2_class_group_exp2_result_top20_count <- left_join(sig_aldex2_class_group_exp2_result_top20, otu_info_rank2_exp2)
clr_class_top20_group_exp2 <- sig_aldex2_class_group_exp2_result_top20_count[, -(2:10)]
rownames(clr_class_top20_group_exp2) <- clr_class_top20_group_exp2$OTU
clr_class_top20_group_exp2 <- clr_class_top20_group_exp2[, -1]

sig_aldex2_order_group_exp2_result_count <- left_join(sig_aldex2_order_group_exp2_result, otu_info_rank3_exp2)
write.csv(sig_aldex2_order_group_exp2_result_count, "aldex/sig_aldex2_order_group_exp2_result_count.csv")
clr_order_group_exp2 <- sig_aldex2_order_group_exp2_result_count[, -(2:10)]
rownames(clr_order_group_exp2) <- clr_order_group_exp2$OTU
clr_order_group_exp2 <- clr_order_group_exp2[, -1]

sig_aldex2_order_group_exp2_result_top20_count <- left_join(sig_aldex2_order_group_exp2_result_top20, otu_info_rank3_exp2)
clr_order_top20_group_exp2 <- sig_aldex2_order_group_exp2_result_top20_count[, -(2:10)]
rownames(clr_order_top20_group_exp2) <- clr_order_top20_group_exp2$OTU
clr_order_top20_group_exp2 <- clr_order_top20_group_exp2[, -1]

sig_aldex2_fam_group_exp2_result_count <- left_join(sig_aldex2_fam_group_exp2_result, otu_info_rank4_exp2)
write.csv(sig_aldex2_fam_group_exp2_result_count, "aldex/sig_aldex2_fam_group_exp2_result_count.csv")
clr_fam_group_exp2 <- sig_aldex2_fam_group_exp2_result_count[, -(2:10)]
rownames(clr_fam_group_exp2) <- clr_fam_group_exp2$OTU
clr_fam_group_exp2 <- clr_fam_group_exp2[, -1]

sig_aldex2_fam_group_exp2_result_top20_count <- left_join(sig_aldex2_fam_group_exp2_result_top20, otu_info_rank4_exp2)
clr_fam_top20_group_exp2 <- sig_aldex2_fam_group_exp2_result_top20_count[, -(2:10)]
rownames(clr_fam_top20_group_exp2) <- clr_fam_top20_group_exp2$OTU
clr_fam_top20_group_exp2 <- clr_fam_top20_group_exp2[, -1]

sig_aldex2_genus_group_exp2_result_count <- left_join(sig_aldex2_genus_group_exp2_result, otu_info_rank5_exp2)
write.csv(sig_aldex2_genus_group_exp2_result_count, "aldex/sig_aldex2_genus_group_exp2_result_count.csv")
clr_genus_group_exp2 <- sig_aldex2_genus_group_exp2_result_count[, -(2:10)]
rownames(clr_genus_group_exp2) <- clr_genus_group_exp2$OTU
clr_genus_group_exp2 <- clr_genus_group_exp2[, -1]

sig_aldex2_genus_group_exp2_result_top20_count <- left_join(sig_aldex2_genus_group_exp2_result_top20, otu_info_rank5_exp2)
clr_genus_top20_group_exp2 <- sig_aldex2_genus_group_exp2_result_top20_count[, -(2:10)]
rownames(clr_genus_top20_group_exp2) <- clr_genus_top20_group_exp2$OTU
clr_genus_top20_group_exp2 <- clr_genus_top20_group_exp2[, -1]

sig_aldex2_genus_group_exp2_result_top10_count <- left_join(sig_aldex2_genus_group_exp2_result_top10, otu_info_rank5_exp2)
clr_genus_top10_group_exp2 <- sig_aldex2_genus_group_exp2_result_top10_count[, -(2:10)]
rownames(clr_genus_top10_group_exp2) <- clr_genus_top10_group_exp2$OTU
clr_genus_top10_group_exp2 <- clr_genus_top10_group_exp2[, -1]

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20, experiment 3
sig_aldex2_phy_group_exp3_result_count <- left_join(sig_aldex2_phy_group_exp3_result, otu_info_rank1_exp3)
write.csv(sig_aldex2_phy_group_exp3_result_count, "aldex/sig_aldex2_phy_group_exp3_result_count.csv")
clr_phy_group_exp3 <- sig_aldex2_phy_group_exp3_result_count[, -(2:10)]
rownames(clr_phy_group_exp3) <- clr_phy_group_exp3$OTU
clr_phy_group_exp3 <- clr_phy_group_exp3[, -1]

sig_aldex2_phy_group_exp3_result_top10_count <- left_join(sig_aldex2_phy_group_exp3_result_top10, otu_info_rank1_exp3)
clr_phy_top10_group_exp3 <- sig_aldex2_phy_group_exp3_result_top10_count[, -(2:10)]
rownames(clr_phy_top10_group_exp3) <- clr_phy_top10_group_exp3$OTU
clr_phy_top10_group_exp3 <- clr_phy_top10_group_exp3[, -1]

sig_aldex2_class_group_exp3_result_count <- left_join(sig_aldex2_class_group_exp3_result, otu_info_rank2_exp3)
write.csv(sig_aldex2_class_group_exp3_result_count, "aldex/sig_aldex2_class_group_exp3_result_count.csv")
clr_class_group_exp3 <- sig_aldex2_class_group_exp3_result_count[, -(2:10)]
rownames(clr_class_group_exp3) <- clr_class_group_exp3$OTU
clr_class_group_exp3 <- clr_class_group_exp3[, -1]

sig_aldex2_class_group_exp3_result_top20_count <- left_join(sig_aldex2_class_group_exp3_result_top20, otu_info_rank2_exp3)
clr_class_top20_group_exp3 <- sig_aldex2_class_group_exp3_result_top20_count[, -(2:10)]
rownames(clr_class_top20_group_exp3) <- clr_class_top20_group_exp3$OTU
clr_class_top20_group_exp3 <- clr_class_top20_group_exp3[, -1]

sig_aldex2_order_group_exp3_result_count <- left_join(sig_aldex2_order_group_exp3_result, otu_info_rank3_exp3)
write.csv(sig_aldex2_order_group_exp3_result_count, "aldex/sig_aldex2_order_group_exp3_result_count.csv")
clr_order_group_exp3 <- sig_aldex2_order_group_exp3_result_count[, -(2:10)]
rownames(clr_order_group_exp3) <- clr_order_group_exp3$OTU
clr_order_group_exp3 <- clr_order_group_exp3[, -1]

sig_aldex2_order_group_exp3_result_top20_count <- left_join(sig_aldex2_order_group_exp3_result_top20, otu_info_rank3_exp3)
clr_order_top20_group_exp3 <- sig_aldex2_order_group_exp3_result_top20_count[, -(2:10)]
rownames(clr_order_top20_group_exp3) <- clr_order_top20_group_exp3$OTU
clr_order_top20_group_exp3 <- clr_order_top20_group_exp3[, -1]

sig_aldex2_fam_group_exp3_result_count <- left_join(sig_aldex2_fam_group_exp3_result, otu_info_rank4_exp3)
write.csv(sig_aldex2_fam_group_exp3_result_count, "aldex/sig_aldex2_fam_group_exp3_result_count.csv")
clr_fam_group_exp3 <- sig_aldex2_fam_group_exp3_result_count[, -(2:10)]
rownames(clr_fam_group_exp3) <- clr_fam_group_exp3$OTU
clr_fam_group_exp3 <- clr_fam_group_exp3[, -1]

sig_aldex2_fam_group_exp3_result_top20_count <- left_join(sig_aldex2_fam_group_exp3_result_top20, otu_info_rank4_exp3)
clr_fam_top20_group_exp3 <- sig_aldex2_fam_group_exp3_result_top20_count[, -(2:10)]
rownames(clr_fam_top20_group_exp3) <- clr_fam_top20_group_exp3$OTU
clr_fam_top20_group_exp3 <- clr_fam_top20_group_exp3[, -1]

sig_aldex2_genus_group_exp3_result_count <- left_join(sig_aldex2_genus_group_exp3_result, otu_info_rank5_exp3)
write.csv(sig_aldex2_genus_group_exp3_result_count, "aldex/sig_aldex2_genus_group_exp3_result_count.csv")
clr_genus_group_exp3 <- sig_aldex2_genus_group_exp3_result_count[, -(2:10)]
rownames(clr_genus_group_exp3) <- clr_genus_group_exp3$OTU
clr_genus_group_exp3 <- clr_genus_group_exp3[, -1]

sig_aldex2_genus_group_exp3_result_top20_count <- left_join(sig_aldex2_genus_group_exp3_result_top20, otu_info_rank5_exp3)
clr_genus_top20_group_exp3 <- sig_aldex2_genus_group_exp3_result_top20_count[, -(2:10)]
rownames(clr_genus_top20_group_exp3) <- clr_genus_top20_group_exp3$OTU
clr_genus_top20_group_exp3 <- clr_genus_top20_group_exp3[, -1]

sig_aldex2_genus_group_exp3_result_top10_count <- left_join(sig_aldex2_genus_group_exp3_result_top10, otu_info_rank5_exp3)
clr_genus_top10_group_exp3 <- sig_aldex2_genus_group_exp3_result_top10_count[, -(2:10)]
rownames(clr_genus_top10_group_exp3) <- clr_genus_top10_group_exp3$OTU
clr_genus_top10_group_exp3 <- clr_genus_top10_group_exp3[, -1]

# Adjusting zeros on the matrix, all three experiments
clr_phy_group_czm <- cmultRepl(t(clr_phy_group),  label=0, method="CZM")
clr_phy_group_czm_trans <- t(apply(clr_phy_group_czm, 1, function(x){log(x) - mean(log(x))}))
#clr_class_group_czm <- cmultRepl(t(clr_class_group),  label=0, method="CZM")
#clr_class_group_czm_trans <- t(apply(clr_class_group_czm, 1, function(x){log(x) - mean(log(x))}))
#clr_order_group_czm <- cmultRepl(t(clr_order_group),  label=0, method="CZM")
#clr_order_group_czm_trans <- t(apply(clr_order_group_czm, 1, function(x){log(x) - mean(log(x))}))
clr_fam_group_czm <- cmultRepl(t(clr_fam_group),  label=0, method="CZM")
clr_fam_group_czm_trans <- t(apply(clr_fam_group_czm, 1, function(x){log(x) - mean(log(x))}))
clr_fam_group_top30_czm <- cmultRepl(t(clr_fam_group_top30),  label=0, method="CZM")
clr_fam_group_top30_czm_trans <- t(apply(clr_fam_group_top30_czm, 1, function(x){log(x) - mean(log(x))}))
clr_genus_group_czm <- cmultRepl(t(clr_genus_group),  label=0, method="CZM")
clr_genus_group_czm_trans <- t(apply(clr_genus_group_czm, 1, function(x){log(x) - mean(log(x))}))
clr_genus_group_top30_czm <- cmultRepl(t(clr_genus_group_top30),  label=0, method="CZM")
clr_genus_group_top30_czm_trans <- t(apply(clr_genus_group_top30_czm, 1, function(x){log(x) - mean(log(x))}))

# Adjusting zeros on the matrix, experiment 1
clr_phy_group_exp1_czm <- cmultRepl(t(clr_phy_group_exp1),  label=0, method="CZM")
clr_phy_group_exp1_czm_trans <- t(apply(clr_phy_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_phy_top10_group_exp1_czm <- cmultRepl(t(clr_phy_top10_group_exp1),  label=0, method="CZM")
clr_phy_top10_group_exp1_czm_trans <- t(apply(clr_phy_top10_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_class_group_exp1_czm <- cmultRepl(t(clr_class_group_exp1),  label=0, method="CZM")
clr_class_group_exp1_czm_trans <- t(apply(clr_class_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_group_exp1_czm <- cmultRepl(t(clr_order_group_exp1),  label=0, method="CZM")
clr_order_group_exp1_czm_trans <- t(apply(clr_order_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_top20_group_exp1_czm <- cmultRepl(t(clr_order_top20_group_exp1),  label=0, method="CZM")
clr_order_top20_group_exp1_czm_trans <- t(apply(clr_order_top20_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_group_exp1_czm <- cmultRepl(t(clr_fam_group_exp1),  label=0, method="CZM")
clr_fam_group_exp1_czm_trans <- t(apply(clr_fam_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_top20_group_exp1_czm <- cmultRepl(t(clr_fam_top20_group_exp1),  label=0, method="CZM")
clr_fam_top20_group_exp1_czm_trans <- t(apply(clr_fam_top20_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_group_exp1_czm <- cmultRepl(t(clr_genus_group_exp1),  label=0, method="CZM")
clr_genus_group_exp1_czm_trans <- t(apply(clr_genus_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top20_group_exp1_czm <- cmultRepl(t(clr_genus_top20_group_exp1),  label=0, method="CZM")
clr_genus_top20_group_exp1_czm_trans <- t(apply(clr_genus_top20_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top10_group_exp1_czm <- cmultRepl(t(clr_genus_top10_group_exp1),  label=0, method="CZM")
clr_genus_top10_group_exp1_czm_trans <- t(apply(clr_genus_top10_group_exp1_czm, 1, function(x){log(x) - mean(log(x))}))

# Adjusting zeros on the matrix, experiment 2
clr_phy_group_exp2_czm <- cmultRepl(t(clr_phy_group_exp2),  label=0, method="CZM")
clr_phy_group_exp2_czm_trans <- t(apply(clr_phy_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_phy_top10_group_exp2_czm <- cmultRepl(t(clr_phy_top10_group_exp2),  label=0, method="CZM")
clr_phy_top10_group_exp2_czm_trans <- t(apply(clr_phy_top10_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_class_group_exp2_czm <- cmultRepl(t(clr_class_group_exp2),  label=0, method="CZM")
clr_class_group_exp2_czm_trans <- t(apply(clr_class_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_class_top20_group_exp2_czm <- cmultRepl(t(clr_class_top20_group_exp2),  label=0, method="CZM")
clr_class_top20_group_exp2_czm_trans <- t(apply(clr_class_top20_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_group_exp2_czm <- cmultRepl(t(clr_order_group_exp2),  label=0, method="CZM")
clr_order_group_exp2_czm_trans <- t(apply(clr_order_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_top20_group_exp2_czm <- cmultRepl(t(clr_order_top20_group_exp2),  label=0, method="CZM")
clr_order_top20_group_exp2_czm_trans <- t(apply(clr_order_top20_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_group_exp2_czm <- cmultRepl(t(clr_fam_group_exp2),  label=0, method="CZM")
clr_fam_group_exp2_czm_trans <- t(apply(clr_fam_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_top20_group_exp2_czm <- cmultRepl(t(clr_fam_top20_group_exp2),  label=0, method="CZM")
clr_fam_top20_group_exp2_czm_trans <- t(apply(clr_fam_top20_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_group_exp2_czm <- cmultRepl(t(clr_genus_group_exp2),  label=0, method="CZM")
clr_genus_group_exp2_czm_trans <- t(apply(clr_genus_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top20_group_exp2_czm <- cmultRepl(t(clr_genus_top20_group_exp2),  label=0, method="CZM")
clr_genus_top20_group_exp2_czm_trans <- t(apply(clr_genus_top20_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top10_group_exp2_czm <- cmultRepl(t(clr_genus_top10_group_exp2),  label=0, method="CZM")
clr_genus_top10_group_exp2_czm_trans <- t(apply(clr_genus_top10_group_exp2_czm, 1, function(x){log(x) - mean(log(x))}))

# Adjusting zeros on the matrix, experiment 3
clr_phy_group_exp3_czm <- cmultRepl(t(clr_phy_group_exp3),  label=0, method="CZM")
clr_phy_group_exp3_czm_trans <- t(apply(clr_phy_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_phy_top10_group_exp3_czm <- cmultRepl(t(clr_phy_top10_group_exp3),  label=0, method="CZM")
clr_phy_top10_group_exp3_czm_trans <- t(apply(clr_phy_top10_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_class_group_exp3_czm <- cmultRepl(t(clr_class_group_exp3),  label=0, method="CZM")
clr_class_group_exp3_czm_trans <- t(apply(clr_class_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_class_top20_group_exp3_czm <- cmultRepl(t(clr_class_top20_group_exp3),  label=0, method="CZM")
clr_class_top20_group_exp3_czm_trans <- t(apply(clr_class_top20_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_group_exp3_czm <- cmultRepl(t(clr_order_group_exp3),  label=0, method="CZM")
clr_order_group_exp3_czm_trans <- t(apply(clr_order_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_order_top20_group_exp3_czm <- cmultRepl(t(clr_order_top20_group_exp3),  label=0, method="CZM")
clr_order_top20_group_exp3_czm_trans <- t(apply(clr_order_top20_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_group_exp3_czm <- cmultRepl(t(clr_fam_group_exp3),  label=0, method="CZM")
clr_fam_group_exp3_czm_trans <- t(apply(clr_fam_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_fam_top20_group_exp3_czm <- cmultRepl(t(clr_fam_top20_group_exp3),  label=0, method="CZM")
clr_fam_top20_group_exp3_czm_trans <- t(apply(clr_fam_top20_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_group_exp3_czm <- cmultRepl(t(clr_genus_group_exp3),  label=0, method="CZM")
clr_genus_group_exp3_czm_trans <- t(apply(clr_genus_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top20_group_exp3_czm <- cmultRepl(t(clr_genus_top20_group_exp3),  label=0, method="CZM")
clr_genus_top20_group_exp3_czm_trans <- t(apply(clr_genus_top20_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

clr_genus_top10_group_exp3_czm <- cmultRepl(t(clr_genus_top10_group_exp3),  label=0, method="CZM")
clr_genus_top10_group_exp3_czm_trans <- t(apply(clr_genus_top10_group_exp3_czm, 1, function(x){log(x) - mean(log(x))}))

############ Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(clr_phy_group_czm_trans)
clr_phy_group_czm_trans_trav <- t(clr_phy_group_czm_trans)
clr_phy_group_czm_trans[, order(colnames(clr_phy_group_czm_trans))]
heatmap(clr_phy_group_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)
head(clr_fam_group_czm_trans)
clr_fam_group_czm_trans_trav <- t(clr_fam_group_czm_trans)
clr_fam_group_czm_trans[, order(colnames(clr_fam_group_czm_trans))]
heatmap(clr_fam_group_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)
head(clr_fam_group_top30_czm_trans)
clr_fam_group_top30_czm_trans_trav <- t(clr_fam_group_top30_czm_trans)
clr_fam_group_top30_czm_trans[, order(colnames(clr_fam_group_top30_czm_trans))]
heatmap(clr_fam_group_top30_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)

head(clr_genus_group_czm_trans)
clr_genus_group_czm_trans_trav <- t(clr_genus_group_czm_trans)
clr_genus_group_czm_trans[, order(colnames(clr_genus_group_czm_trans))]
heatmap(clr_genus_group_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)
head(clr_genus_group_top30_czm_trans)
clr_genus_group_top30_czm_trans_trav <- t(clr_genus_group_top30_czm_trans)
clr_genus_group_top30_czm_trans[, order(colnames(clr_genus_group_top30_czm_trans))]
heatmap(clr_genus_group_top30_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
z_score_phy_group <- scale(clr_phy_group_czm_trans_trav)
#z_score_class_group <- scale(clr_class_group_czm_trans_trav)
#z_score_order_group <- scale(clr_order_group_czm_trans_trav)
z_score_fam_group <- scale(clr_fam_group_czm_trans_trav)
z_score_fam_group_top30 <- scale(clr_fam_group_top30_czm_trans_trav)
z_score_genus_group <- scale(clr_genus_group_czm_trans_trav)
z_score_genus_group_top30 <- scale(clr_genus_group_top30_czm_trans_trav)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top = HeatmapAnnotation(border = TRUE,
  Experiment = as.vector(sample_tab_group_factors$Experiment),
  Sample = as.vector(sample_tab_group_factors$SampleType),
  col = list(
    Experiment = c("Experiment1" = "#e0ecf4", "Experiment2" = "#9ebcda", "Experiment3" = "#8856a7"),
    Sample = c("EmptyGut" = "#1B9E77", "FullGut" = "#D59D08", "FecalPellets24" = "#7570B3",
                "FecalPellets2" = "#D95F02", "Water" = "#666666")
  ),
  annotation_legend_param = list(
    Experiment = list(
      title = "Experiment",
      at = c("Experiment1", "Experiment2", "Experiment3"),
      labels = c ("Experiment 1", "Experiment 2", "Experiment 3")
    ),
    Sample = list(
      title = "Sample Type",
      at = c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water"),
      labels = c ("EG", "FG", "FP24Hrs", "FP2Hrs", "SW")
    )
  ))

# Organize total abundance and taxa name
z_score_phy_group_name <- as.data.frame(z_score_phy_group)
str(z_score_phy_group_name)
z_score_phy_group_name <- rownames_to_column(z_score_phy_group_name, var = "OTU")
z_score_phy_group_name <- left_join(z_score_phy_group_name, taxa_info)
head(z_score_phy_group_name)
z_score_phy_group_name <- z_score_phy_group_name[, -(2:112)]

z_score_fam_group_name <- as.data.frame(z_score_fam_group)
str(z_score_fam_group_name)
z_score_fam_group_name <- rownames_to_column(z_score_fam_group_name, var = "OTU")
z_score_fam_group_name <- left_join(z_score_fam_group_name, taxa_info)
head(z_score_fam_group_name)
z_score_fam_group_name <- z_score_fam_group_name[, -(2:112)]

z_score_fam_group_top30_name <- as.data.frame(z_score_fam_group_top30)
str(z_score_fam_group_top30_name)
z_score_fam_group_top30_name <- rownames_to_column(z_score_fam_group_top30_name, var = "OTU")
z_score_fam_group_top30_name <- left_join(z_score_fam_group_top30_name, taxa_info)
head(z_score_fam_group_top30_name)
z_score_fam_group_top30_name <- z_score_fam_group_top30_name[, -(2:112)]

z_score_genus_group_name <- as.data.frame(z_score_genus_group)
str(z_score_genus_group_name)
z_score_genus_group_name <- rownames_to_column(z_score_genus_group_name, var = "OTU")
z_score_genus_group_name <- left_join(z_score_genus_group_name, taxa_info)
head(z_score_genus_group_name)
z_score_genus_group_name <- z_score_genus_group_name[, -(2:112)]

z_score_genus_group_top30_name <- as.data.frame(z_score_genus_group_top30)
str(z_score_genus_group_top30_name)
z_score_genus_group_top30_name <- rownames_to_column(z_score_genus_group_top30_name, var = "OTU")
z_score_genus_group_top30_name <- left_join(z_score_genus_group_top30_name, taxa_info)
head(z_score_genus_group_top30_name)
z_score_genus_group_top30_name <- z_score_genus_group_top30_name[, -(2:112)]

otu_count_total <- as.data.frame(otu_count)
otu_count_total$Total <- rowSums(otu_count_total[, -1])
head(otu_count_total)
otu_count_total <- otu_count_total[, -(2:112)]

z_score_phy_group_count_total <- left_join(z_score_phy_group_name, otu_count_total)
head(z_score_phy_group_count_total)
z_score_fam_group_count_total <- left_join(z_score_fam_group_name, otu_count_total)
head(z_score_fam_group_count_total)
z_score_fam_group_top30_count_total <- left_join(z_score_fam_group_top30_name, otu_count_total)
head(z_score_fam_group_top30_count_total)
z_score_genus_group_count_total <- left_join(z_score_genus_group_name, otu_count_total)
head(z_score_genus_group_count_total)
z_score_genus_group_top30_count_total <- left_join(z_score_genus_group_top30_name, otu_count_total)
head(z_score_genus_group_top30_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 300, 3000, 30000, 300000),
                               c("#f5f5f5",
                                 "#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))
ha_right_phy = rowAnnotation(
  Abundance = z_score_phy_group_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_phy = z_score_phy_group_count_total$Phylum

ha_right_fam = rowAnnotation(
  Abundance = z_score_fam_group_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_fam = z_score_fam_group_count_total$Family

ha_right_fam_top30 = rowAnnotation(
  Abundance = z_score_fam_group_top30_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_fam_top30 = z_score_fam_group_top30_count_total$Family

ha_right_genus = rowAnnotation(
  Abundance = z_score_genus_group_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_genus = z_score_genus_group_count_total$Genus

ha_right_genus_top30 = rowAnnotation(
  Abundance = z_score_genus_group_top30_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_genus_top30 = z_score_genus_group_top30_count_total$Genus

# Plot heatmap at the phylum level
hm_phylum <- Heatmap(z_score_phy_group, name = "Z-score, CLR", col = col_matrix3,
                     #column_title = "Experiments & Sample Types", 
                     #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     cluster_columns = cluster_within_group(as.vector(sample_tab_group_factors$SampleType)),
                     column_split = as.vector(sample_tab_group_factors$Experiment),
                     border = TRUE,
                     top_annotation = ha_top,
                     right_annotation = ha_right_phy,
                     row_title = "Phylum",
                     row_labels = row_labels_phy,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     row_order = order(row_labels_phy),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_phylum
ggsave("results/16S_heatmap_phy_SampleType_Experiment_new.pdf", width = 12, height = 6, dpi = 150) # save graphic

hm_fam <- Heatmap(z_score_fam_group, name = "Z-score, CLR", col = col_matrix3,
                     #column_title = "Experiments & Sample Types", 
                     #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     #cluster_columns = FALSE,
                     column_split = as.vector(sample_tab_group_factors$Experiment),
                     border = TRUE,
                     top_annotation = ha_top,
                     right_annotation = ha_right_fam,
                     row_title = "Family",
                     row_labels = row_labels_fam,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     row_order = order(row_labels_fam),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_fam
ggsave("results/16S_heatmap_fam_SampleType_Experiment_new.pdf", width = 12, height = 16, dpi = 150) # save graphic

hm_fam_top30 <- Heatmap(z_score_fam_group_top30, name = "Z-score, CLR", col = col_matrix3,
                  #column_title = "Experiments & Sample Types", 
                  #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                  #cluster_columns = FALSE,
                  column_split = as.vector(sample_tab_group_factors$Experiment),
                  border = TRUE,
                  top_annotation = ha_top,
                  right_annotation = ha_right_fam_top30,
                  row_title = "Family",
                  row_labels = row_labels_fam_top30,
                  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8),
                  row_order = order(row_labels_fam_top30),
                  rect_gp = gpar(col = "white", lwd = 1),
                  show_column_names = FALSE,
                  show_heatmap_legend = TRUE)
hm_fam_top30
ggsave("results/16S_heatmap_fam_top30_SampleType_Experiment_new.pdf", width = 12, height = 16, dpi = 150) # save graphic

hm_genus <- Heatmap(z_score_genus_group, name = "Z-score, CLR", col = col_matrix3,
                  #column_title = "Experiments & Sample Types", 
                  #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                  #cluster_columns = FALSE,
                  column_split = as.vector(sample_tab_group_factors$Experiment),
                  border = TRUE,
                  top_annotation = ha_top,
                  right_annotation = ha_right_genus,
                  row_title = "Genus",
                  row_labels = row_labels_genus,
                  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8),
                  row_order = order(row_labels_genus),
                  rect_gp = gpar(col = "white", lwd = 1),
                  show_column_names = FALSE,
                  show_heatmap_legend = TRUE)
hm_genus
ggsave("results/16S_heatmap_genus_SampleType_Experiment_new.pdf", width = 12, height = 18, dpi = 150) # save graphic

hm_genus_top30 <- Heatmap(z_score_genus_group_top30, name = "Z-score, CLR", col = col_matrix3,
                        #column_title = "Experiments & Sample Types", 
                        #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                        #cluster_columns = FALSE,
                        column_split = as.vector(sample_tab_group_factors$Experiment),
                        border = TRUE,
                        top_annotation = ha_top,
                        right_annotation = ha_right_genus_top30,
                        row_title = "Genus",
                        row_labels = row_labels_genus_top30,
                        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        row_order = order(row_labels_genus_top30),
                        rect_gp = gpar(col = "white", lwd = 1),
                        show_column_names = FALSE,
                        show_heatmap_legend = TRUE)
hm_genus_top30
ggsave("results/16S_heatmap_genus_top30_SampleType_Experiment_new.pdf", width = 12, height = 6, dpi = 150) # save graphic
#################### Plot heatmap with ggplot2 ####################
otu_count_total <- as.data.frame(otu_count)
otu_count_total$Total <- rowSums(otu_count_total[, -1])
head(otu_count_total)
otu_count_total <- otu_count_total[, -(2:112)]

# Apply scale to calculate Z-score on the clr_transf zero adjusted objects, then tranform from matrix to data.frame, all three experiments
z_score_phy_group_df <- as.data.frame(z_score_phy_group)
z_score_phy_group_df <- left_join(z_score_phy_group_df, otu_count_total)
rownames(z_score_phy_group_df) <- z_score_phy_group_df$OTU
z_score_phy_group_df <- z_score_phy_group_df[, -1]
z_score_phy_group_df <-as.data.frame(t(z_score_phy_group_df))

# Use rownames to create another variable called sample
z_score_phy_group_df <- tibble::rownames_to_column(z_score_phy_group_df, "Sample")

#Add sample factors (Experiment, SampleType) to the Z.score object
sample_tab_group_factors <- rownames_to_column(sample_tab_group_factors, var = "Sample")
z_score_phy_group_df_var <- data.frame(full_join(z_score_phy_group_df, sample_tab_group_factors))

#Add sample factors (Experiment, SampleType) to the Z.score object
#Remove all categorical variables
z_score_phy_group_df_long <- pivot_longer(data = z_score_phy_group_df_var,
                                       cols = -c(Sample, SampleType,Experiment),
                                       names_to = "ASV",
                                       values_to = "Z")
class(z_score_phy_group_df_long)
z_score_phy_group_df_long <- as.data.frame(z_score_phy_group_df_long)
class(z_score_phy_group_df_long)
str(z_score_phy_group_df_long)
z_score_phy_group_df_long_merged <- merge(z_score_phy_group_df_long, taxa_info, by.x = 4, by.y = 1, all.x = TRUE)

# Change facet SampleType labels for plotting on heatmaps
sample.type.labs <- c("EG", "FG", "FP24Hrs", "FP2Hrs", "SW")
names(sample.type.labs) <- c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water")

experiment.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(experiment.labs) <- c("Experiment1", "Experiment2", "Experiment3")

#Heatmaps all three experiments
theme_set(theme_bw())
heatmap_phy <- ggplot(z_score_phy_group_df_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile(data = filter(z_score_phy_group_df_long_merged, Phylum != "Total", Sample != "Total"), aes(fill = Z)) +
  geom_point(data = filter(z_score_phy_group_df_long_merged, Phylum == "Total" | Sample == "Total"), aes(color = Z), size = 3, shape = 19) +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "Sample") + # Add a nicer x-axis title
  facet_nested(. ~Experiment+SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs, Experiment = experiment.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy
ggsave("results/16S_heatmap_phy_SampleType_Experiment.pdf", width = 11, height = 4, dpi = 150) # save graphic
ggsave("results/16S_heatmap_phy_SampleType.pdf", width = 6, height = 4, dpi = 150) # save graphic

# Apply scale to calculate Z-score on the clr_transf zero adjusted objects, then tranform from matrix to data.frame, experiment 1
z_score_phy_group_exp1 <- as.data.frame(scale(clr_phy_group_exp1_czm_trans))
z_score_phy_top10_group_exp1 <- as.data.frame(scale(clr_phy_top10_group_exp1_czm_trans))
z_score_class_group_exp1 <- as.data.frame(scale(clr_class_group_exp1_czm_trans))
z_score_order_group_exp1 <- as.data.frame(scale(clr_order_group_exp1_czm_trans))
z_score_order_top20_group_exp1 <- as.data.frame(scale(clr_order_top20_group_exp1_czm_trans))
z_score_fam_group_exp1 <- as.data.frame(scale(clr_fam_group_exp1_czm_trans))
z_score_fam_top20_group_exp1 <- as.data.frame(scale(clr_fam_top20_group_exp1_czm_trans))
z_score_genus_group_exp1 <- as.data.frame(scale(clr_genus_group_exp1_czm_trans))
z_score_genus_top20_group_exp1 <- as.data.frame(scale(clr_genus_top20_group_exp1_czm_trans))
z_score_genus_top10_group_exp1 <- as.data.frame(scale(clr_genus_top10_group_exp1_czm_trans))

# Apply scale to calculate Z-score on the clr_transf zero adjusted objects, then tranform from matrix to data.frame, experiment 2
z_score_phy_group_exp2 <- as.data.frame(scale(clr_phy_group_exp2_czm_trans))
z_score_phy_top10_group_exp2 <- as.data.frame(scale(clr_phy_top10_group_exp2_czm_trans))
z_score_class_group_exp2 <- as.data.frame(scale(clr_class_group_exp2_czm_trans))
z_score_class_top20_group_exp2 <- as.data.frame(scale(clr_class_top20_group_exp2_czm_trans))
z_score_order_group_exp2 <- as.data.frame(scale(clr_order_group_exp2_czm_trans))
z_score_order_top20_group_exp2 <- as.data.frame(scale(clr_order_top20_group_exp2_czm_trans))
z_score_fam_group_exp2 <- as.data.frame(scale(clr_fam_group_exp2_czm_trans))
z_score_fam_top20_group_exp2 <- as.data.frame(scale(clr_fam_top20_group_exp2_czm_trans))
z_score_genus_group_exp2 <- as.data.frame(scale(clr_genus_group_exp2_czm_trans))
z_score_genus_top20_group_exp2 <- as.data.frame(scale(clr_genus_top20_group_exp2_czm_trans))
z_score_genus_top10_group_exp2 <- as.data.frame(scale(clr_genus_top10_group_exp2_czm_trans))

# Apply scale to calculate Z-score on the clr_transf zero adjusted objects, then transform from matrix to data.frame, experiment 3
z_score_phy_group_exp3 <- as.data.frame(scale(clr_phy_group_exp3_czm_trans))
z_score_phy_top10_group_exp3 <- as.data.frame(scale(clr_phy_top10_group_exp3_czm_trans))
z_score_class_group_exp3 <- as.data.frame(scale(clr_class_group_exp3_czm_trans))
z_score_class_top20_group_exp3 <- as.data.frame(scale(clr_class_top20_group_exp3_czm_trans))
z_score_order_group_exp3 <- as.data.frame(scale(clr_order_group_exp3_czm_trans))
z_score_order_top20_group_exp3 <- as.data.frame(scale(clr_order_top20_group_exp3_czm_trans))
z_score_fam_group_exp3 <- as.data.frame(scale(clr_fam_group_exp3_czm_trans))
z_score_fam_top20_group_exp3 <- as.data.frame(scale(clr_fam_top20_group_exp3_czm_trans))
z_score_genus_group_exp3 <- as.data.frame(scale(clr_genus_group_exp3_czm_trans))
z_score_genus_top20_group_exp3 <- as.data.frame(scale(clr_genus_top20_group_exp3_czm_trans))
z_score_genus_top10_group_exp3 <- as.data.frame(scale(clr_genus_top10_group_exp3_czm_trans))

# Use rownames to create another variable called Sample, all three experiments
z_score_phy_group_df <- tibble::rownames_to_column(z_score_phy_group_df, "OTU")

# Use rownames to create another variable called Sample, experiment 1
z_score_phy_group_exp1 <- tibble::rownames_to_column(z_score_phy_group_exp1, "Sample")
z_score_phy_top10_group_exp1 <- tibble::rownames_to_column(z_score_phy_top10_group_exp1, "Sample")
z_score_class_group_exp1 <- tibble::rownames_to_column(z_score_class_group_exp1, "Sample")
z_score_order_group_exp1 <- tibble::rownames_to_column(z_score_order_group_exp1, "Sample")
z_score_order_top20_group_exp1 <- tibble::rownames_to_column(z_score_order_top20_group_exp1, "Sample")
z_score_fam_group_exp1 <- tibble::rownames_to_column(z_score_fam_group_exp1, "Sample")
z_score_fam_top20_group_exp1 <- tibble::rownames_to_column(z_score_fam_top20_group_exp1, "Sample")
z_score_genus_group_exp1 <- tibble::rownames_to_column(z_score_genus_group_exp1, "Sample")
z_score_genus_top20_group_exp1 <- tibble::rownames_to_column(z_score_genus_top20_group_exp1, "Sample")
z_score_genus_top10_group_exp1 <- tibble::rownames_to_column(z_score_genus_top10_group_exp1, "Sample")

# Use rownames to create another variable called Sample, experiment 2
z_score_phy_group_exp2 <- tibble::rownames_to_column(z_score_phy_group_exp2, "Sample")
z_score_phy_top10_group_exp2 <- tibble::rownames_to_column(z_score_phy_top10_group_exp2, "Sample")
z_score_class_group_exp2 <- tibble::rownames_to_column(z_score_class_group_exp2, "Sample")
z_score_class_top20_group_exp2 <- tibble::rownames_to_column(z_score_class_top20_group_exp2, "Sample")
z_score_order_group_exp2 <- tibble::rownames_to_column(z_score_order_group_exp2, "Sample")
z_score_order_top20_group_exp2 <- tibble::rownames_to_column(z_score_order_top20_group_exp2, "Sample")
z_score_fam_group_exp2 <- tibble::rownames_to_column(z_score_fam_group_exp2, "Sample")
z_score_fam_top20_group_exp2 <- tibble::rownames_to_column(z_score_fam_top20_group_exp2, "Sample")
z_score_genus_group_exp2 <- tibble::rownames_to_column(z_score_genus_group_exp2, "Sample")
z_score_genus_top20_group_exp2 <- tibble::rownames_to_column(z_score_genus_top20_group_exp2, "Sample")
z_score_genus_top10_group_exp2 <- tibble::rownames_to_column(z_score_genus_top10_group_exp2, "Sample")

# Use rownames to create another variable called Sample, experiment 3
z_score_phy_group_exp3 <- tibble::rownames_to_column(z_score_phy_group_exp3, "Sample")
z_score_phy_top10_group_exp3 <- tibble::rownames_to_column(z_score_phy_top10_group_exp3, "Sample")
z_score_class_group_exp3 <- tibble::rownames_to_column(z_score_class_group_exp3, "Sample")
z_score_class_top20_group_exp3 <- tibble::rownames_to_column(z_score_class_top20_group_exp3, "Sample")
z_score_order_group_exp3 <- tibble::rownames_to_column(z_score_order_group_exp3, "Sample")
z_score_order_top20_group_exp3 <- tibble::rownames_to_column(z_score_order_top20_group_exp3, "Sample")
z_score_fam_group_exp3 <- tibble::rownames_to_column(z_score_fam_group_exp3, "Sample")
z_score_fam_top20_group_exp3 <- tibble::rownames_to_column(z_score_fam_top20_group_exp3, "Sample")
z_score_genus_group_exp3 <- tibble::rownames_to_column(z_score_genus_group_exp3, "Sample")
z_score_genus_top20_group_exp3 <- tibble::rownames_to_column(z_score_genus_top20_group_exp3, "Sample")
z_score_genus_top10_group_exp3 <- tibble::rownames_to_column(z_score_genus_top10_group_exp3, "Sample")

#Add sample factors (SampleType, Experiment) to the Z.score object, all three experiments
z_score_phy_group <- data.frame(cbind(z_score_phy_group, sample_tab_group_factors))

#Add sample factors (SampleType, Experiment) to the Z.score object, experiment 1
z_score_phy_group_exp1 <- data.frame(cbind(z_score_phy_group_exp1, sample_tab_group_exp1_factors))
z_score_phy_top10_group_exp1 <- data.frame(cbind(z_score_phy_top10_group_exp1, sample_tab_group_exp1_factors))
z_score_class_group_exp1 <- data.frame(cbind(z_score_class_group_exp1, sample_tab_group_exp1_factors))
z_score_order_group_exp1 <- data.frame(cbind(z_score_order_group_exp1, sample_tab_group_exp1_factors))
z_score_order_top20_group_exp1 <- data.frame(cbind(z_score_order_top20_group_exp1, sample_tab_group_exp1_factors))
z_score_fam_group_exp1 <- data.frame(cbind(z_score_fam_group_exp1, sample_tab_group_exp1_factors))
z_score_fam_top20_group_exp1 <- data.frame(cbind(z_score_fam_top20_group_exp1, sample_tab_group_exp1_factors))
z_score_genus_group_exp1 <- data.frame(cbind(z_score_genus_group_exp1, sample_tab_group_exp1_factors))
z_score_genus_top20_group_exp1 <- data.frame(cbind(z_score_genus_top20_group_exp1, sample_tab_group_exp1_factors))
z_score_genus_top10_group_exp1 <- data.frame(cbind(z_score_genus_top10_group_exp1, sample_tab_group_exp1_factors))

#Add sample factors (SampleType, Experiment) to the Z.score object, experiment 2
z_score_phy_group_exp2 <- data.frame(cbind(z_score_phy_group_exp2, sample_tab_group_exp2_factors))
z_score_phy_top10_group_exp2 <- data.frame(cbind(z_score_phy_top10_group_exp2, sample_tab_group_exp2_factors))
z_score_class_group_exp2 <- data.frame(cbind(z_score_class_group_exp2, sample_tab_group_exp2_factors))
z_score_class_top20_group_exp2 <- data.frame(cbind(z_score_class_top20_group_exp2, sample_tab_group_exp2_factors))
z_score_order_group_exp2 <- data.frame(cbind(z_score_order_group_exp2, sample_tab_group_exp2_factors))
z_score_order_top20_group_exp2 <- data.frame(cbind(z_score_order_top20_group_exp2, sample_tab_group_exp2_factors))
z_score_fam_group_exp2 <- data.frame(cbind(z_score_fam_group_exp2, sample_tab_group_exp2_factors))
z_score_fam_top20_group_exp2 <- data.frame(cbind(z_score_fam_top20_group_exp2, sample_tab_group_exp2_factors))
z_score_genus_group_exp2 <- data.frame(cbind(z_score_genus_group_exp2, sample_tab_group_exp2_factors))
z_score_genus_top20_group_exp2 <- data.frame(cbind(z_score_genus_top20_group_exp2, sample_tab_group_exp2_factors))
z_score_genus_top10_group_exp2 <- data.frame(cbind(z_score_genus_top10_group_exp2, sample_tab_group_exp2_factors))

#Add sample factors (SampleType, Experiment) to the Z.score object, experiment 3
z_score_phy_group_exp3 <- data.frame(cbind(z_score_phy_group_exp3, sample_tab_group_exp3_factors))
z_score_phy_top10_group_exp3 <- data.frame(cbind(z_score_phy_top10_group_exp3, sample_tab_group_exp3_factors))
z_score_class_group_exp3 <- data.frame(cbind(z_score_class_group_exp3, sample_tab_group_exp3_factors))
z_score_class_top20_group_exp3 <- data.frame(cbind(z_score_class_top20_group_exp3, sample_tab_group_exp3_factors))
z_score_order_group_exp3 <- data.frame(cbind(z_score_order_group_exp3, sample_tab_group_exp3_factors))
z_score_order_top20_group_exp3 <- data.frame(cbind(z_score_order_top20_group_exp3, sample_tab_group_exp3_factors))
z_score_fam_group_exp3 <- data.frame(cbind(z_score_fam_group_exp3, sample_tab_group_exp3_factors))
z_score_fam_top20_group_exp3 <- data.frame(cbind(z_score_fam_top20_group_exp3, sample_tab_group_exp3_factors))
z_score_genus_group_exp3 <- data.frame(cbind(z_score_genus_group_exp3, sample_tab_group_exp3_factors))
z_score_genus_top20_group_exp3 <- data.frame(cbind(z_score_genus_top20_group_exp3, sample_tab_group_exp3_factors))
z_score_genus_top10_group_exp3 <- data.frame(cbind(z_score_genus_top10_group_exp3, sample_tab_group_exp3_factors))

#Add sample factors (Sample, SampleType, Experiment) to the Z.score object, experiment 1
z_score_phy_group_long <- pivot_longer(data = z_score_phy_group,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_phy_group_long <- as.data.frame(z_score_phy_group_long)
z_score_phy_group_long_merged <- merge(z_score_phy_group_long, sig_aldex2_phy_group_result, by.x = 4, by.y = 1, all.x = TRUE)
z_score_phy_group_long_merged <- merge(z_score_phy_group_long_merged, otu_count_total, by.x = 1, by.y = 1, all.x = TRUE)

#Add sample factors (Sample, SampleType, Experiment) to the Z.score object, experiment 1
z_score_phy_group_exp1_long <- pivot_longer(data = z_score_phy_group_exp1,
                                       cols = -c(Sample, SampleType, Experiment),
                                       names_to = "ASV",
                                       values_to = "Z")
z_score_phy_group_exp1_long <- as.data.frame(z_score_phy_group_exp1_long)
z_score_phy_group_exp1_long_merged <- merge(z_score_phy_group_exp1_long, sig_aldex2_phy_group_exp1_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_phy_top10_group_exp1_long <- pivot_longer(data = z_score_phy_top10_group_exp1,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_phy_top10_group_exp1_long <- as.data.frame(z_score_phy_top10_group_exp1_long)
z_score_phy_top10_group_exp1_long_merged <- merge(z_score_phy_top10_group_exp1_long, sig_aldex2_phy_group_exp1_result_top10, by.x = 4, by.y = 1, all.x = TRUE)


z_score_class_group_exp1_long <- pivot_longer(data = z_score_class_group_exp1,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_class_group_exp1_long <- as.data.frame(z_score_class_group_exp1_long)
z_score_class_group_exp1_long_merged <- merge(z_score_class_group_exp1_long, sig_aldex2_class_group_exp1_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_group_exp1_long <- pivot_longer(data = z_score_order_group_exp1,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_order_group_exp1_long <- as.data.frame(z_score_order_group_exp1_long)
z_score_order_group_exp1_long_merged <- merge(z_score_order_group_exp1_long, sig_aldex2_order_group_exp1_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_top20_group_exp1_long <- pivot_longer(data = z_score_order_top20_group_exp1,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_order_top20_group_exp1_long <- as.data.frame(z_score_order_top20_group_exp1_long)
z_score_order_top20_group_exp1_long_merged <- merge(z_score_order_top20_group_exp1_long, sig_aldex2_order_group_exp1_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_group_exp1_long <- pivot_longer(data = z_score_fam_group_exp1,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_fam_group_exp1_long <- as.data.frame(z_score_fam_group_exp1_long)
z_score_fam_group_exp1_long_merged <- merge(z_score_fam_group_exp1_long, sig_aldex2_fam_group_exp1_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_top20_group_exp1_long <- pivot_longer(data = z_score_fam_top20_group_exp1,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_fam_top20_group_exp1_long <- as.data.frame(z_score_fam_top20_group_exp1_long)
z_score_fam_top20_group_exp1_long_merged <- merge(z_score_fam_top20_group_exp1_long, sig_aldex2_fam_group_exp1_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_group_exp1_long <- pivot_longer(data = z_score_genus_group_exp1,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_genus_group_exp1_long <- as.data.frame(z_score_genus_group_exp1_long)
z_score_genus_group_exp1_long_merged <- merge(z_score_genus_group_exp1_long, sig_aldex2_genus_group_exp1_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top20_group_exp1_long <- pivot_longer(data = z_score_genus_top20_group_exp1,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_genus_top20_group_exp1_long <- as.data.frame(z_score_genus_top20_group_exp1_long)
z_score_genus_top20_group_exp1_long_merged <- merge(z_score_genus_top20_group_exp1_long, sig_aldex2_genus_group_exp1_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top10_group_exp1_long <- pivot_longer(data = z_score_genus_top10_group_exp1,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_genus_top10_group_exp1_long <- as.data.frame(z_score_genus_top10_group_exp1_long)
z_score_genus_top10_group_exp1_long_merged <- merge(z_score_genus_top10_group_exp1_long, sig_aldex2_genus_group_exp1_result_top10, by.x = 4, by.y = 1, all.x = TRUE)

#Add sample factors (Sample, SampleType, Experiment) to the Z.score object, experiment 2
z_score_phy_group_exp2_long <- pivot_longer(data = z_score_phy_group_exp2,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_phy_group_exp2_long <- as.data.frame(z_score_phy_group_exp2_long)
z_score_phy_group_exp2_long_merged <- merge(z_score_phy_group_exp2_long, sig_aldex2_phy_group_exp2_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_phy_top10_group_exp2_long <- pivot_longer(data = z_score_phy_top10_group_exp2,
                                                  cols = -c(Sample, SampleType, Experiment),
                                                  names_to = "ASV",
                                                  values_to = "Z")
z_score_phy_top10_group_exp2_long <- as.data.frame(z_score_phy_top10_group_exp2_long)
z_score_phy_top10_group_exp2_long_merged <- merge(z_score_phy_top10_group_exp2_long, sig_aldex2_phy_group_exp2_result_top10, by.x = 4, by.y = 1, all.x = TRUE)

z_score_class_group_exp2_long <- pivot_longer(data = z_score_class_group_exp2,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_class_group_exp2_long <- as.data.frame(z_score_class_group_exp2_long)
z_score_class_group_exp2_long_merged <- merge(z_score_class_group_exp2_long, sig_aldex2_class_group_exp2_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_class_top20_group_exp2_long <- pivot_longer(data = z_score_class_top20_group_exp2,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_class_top20_group_exp2_long <- as.data.frame(z_score_class_top20_group_exp2_long)
z_score_class_top20_group_exp2_long_merged <- merge(z_score_class_top20_group_exp2_long, sig_aldex2_class_group_exp2_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_group_exp2_long <- pivot_longer(data = z_score_order_group_exp2,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_order_group_exp2_long <- as.data.frame(z_score_order_group_exp2_long)
z_score_order_group_exp2_long_merged <- merge(z_score_order_group_exp2_long, sig_aldex2_order_group_exp2_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_top20_group_exp2_long <- pivot_longer(data = z_score_order_top20_group_exp2,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_order_top20_group_exp2_long <- as.data.frame(z_score_order_top20_group_exp2_long)
z_score_order_top20_group_exp2_long_merged <- merge(z_score_order_top20_group_exp2_long, sig_aldex2_order_group_exp2_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_group_exp2_long <- pivot_longer(data = z_score_fam_group_exp2,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_fam_group_exp2_long <- as.data.frame(z_score_fam_group_exp2_long)
z_score_fam_group_exp2_long_merged <- merge(z_score_fam_group_exp2_long, sig_aldex2_fam_group_exp2_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_top20_group_exp2_long <- pivot_longer(data = z_score_fam_top20_group_exp2,
                                                  cols = -c(Sample, SampleType, Experiment),
                                                  names_to = "ASV",
                                                  values_to = "Z")
z_score_fam_top20_group_exp2_long <- as.data.frame(z_score_fam_top20_group_exp2_long)
z_score_fam_top20_group_exp2_long_merged <- merge(z_score_fam_top20_group_exp2_long, sig_aldex2_fam_group_exp2_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_group_exp2_long <- pivot_longer(data = z_score_genus_group_exp2,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_genus_group_exp2_long <- as.data.frame(z_score_genus_group_exp2_long)
z_score_genus_group_exp2_long_merged <- merge(z_score_genus_group_exp2_long, sig_aldex2_genus_group_exp2_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top20_group_exp2_long <- pivot_longer(data = z_score_genus_top20_group_exp2,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_genus_top20_group_exp2_long <- as.data.frame(z_score_genus_top20_group_exp2_long)
z_score_genus_top20_group_exp2_long_merged <- merge(z_score_genus_top20_group_exp2_long, sig_aldex2_genus_group_exp2_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top10_group_exp2_long <- pivot_longer(data = z_score_genus_top10_group_exp2,
                                                  cols = -c(Sample, SampleType, Experiment),
                                                  names_to = "ASV",
                                                  values_to = "Z")
z_score_genus_top10_group_exp2_long <- as.data.frame(z_score_genus_top10_group_exp2_long)
z_score_genus_top10_group_exp2_long_merged <- merge(z_score_genus_top10_group_exp2_long, sig_aldex2_genus_group_exp2_result_top10, by.x = 4, by.y = 1, all.x = TRUE)

#Add sample factors (Sample, SampleType, Experiment) to the Z.score object, experiment 3
z_score_phy_group_exp3_long <- pivot_longer(data = z_score_phy_group_exp3,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_phy_group_exp3_long <- as.data.frame(z_score_phy_group_exp3_long)
z_score_phy_group_exp3_long_merged <- merge(z_score_phy_group_exp3_long, sig_aldex2_phy_group_exp3_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_phy_top10_group_exp3_long <- pivot_longer(data = z_score_phy_top10_group_exp3,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_phy_top10_group_exp3_long <- as.data.frame(z_score_phy_top10_group_exp3_long)
z_score_phy_top10_group_exp3_long_merged <- merge(z_score_phy_top10_group_exp3_long, sig_aldex2_phy_group_exp3_result_top10, by.x = 4, by.y = 1, all.x = TRUE)


z_score_class_group_exp3_long <- pivot_longer(data = z_score_class_group_exp3,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_class_group_exp3_long <- as.data.frame(z_score_class_group_exp3_long)
z_score_class_group_exp3_long_merged <- merge(z_score_class_group_exp3_long, sig_aldex2_class_group_exp3_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_class_top20_group_exp3_long <- pivot_longer(data = z_score_class_top20_group_exp3,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_class_top20_group_exp3_long <- as.data.frame(z_score_class_top20_group_exp3_long)
z_score_class_top20_group_exp3_long_merged <- merge(z_score_class_top20_group_exp3_long, sig_aldex2_class_group_exp3_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_group_exp3_long <- pivot_longer(data = z_score_order_group_exp3,
                                              cols = -c(Sample, SampleType, Experiment),
                                              names_to = "ASV",
                                              values_to = "Z")
z_score_order_group_exp3_long <- as.data.frame(z_score_order_group_exp3_long)
z_score_order_group_exp3_long_merged <- merge(z_score_order_group_exp3_long, sig_aldex2_order_group_exp3_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_order_top20_group_exp3_long <- pivot_longer(data = z_score_order_top20_group_exp3,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_order_top20_group_exp3_long <- as.data.frame(z_score_order_top20_group_exp3_long)
z_score_order_top20_group_exp3_long_merged <- merge(z_score_order_top20_group_exp3_long, sig_aldex2_order_group_exp3_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_group_exp3_long <- pivot_longer(data = z_score_fam_group_exp3,
                                            cols = -c(Sample, SampleType, Experiment),
                                            names_to = "ASV",
                                            values_to = "Z")
z_score_fam_group_exp3_long <- as.data.frame(z_score_fam_group_exp3_long)
z_score_fam_group_exp3_long_merged <- merge(z_score_fam_group_exp3_long, sig_aldex2_fam_group_exp3_result, by.x = 4, by.y = 1, all.x = TRUE)

z_score_fam_top20_group_exp3_long <- pivot_longer(data = z_score_fam_top20_group_exp3,
                                                  cols = -c(Sample, SampleType, Experiment),
                                                  names_to = "ASV",
                                                  values_to = "Z")
z_score_fam_top20_group_exp3_long <- as.data.frame(z_score_fam_top20_group_exp3_long)
z_score_fam_top20_group_exp3_long_merged <- merge(z_score_fam_top20_group_exp3_long, sig_aldex2_fam_group_exp3_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top20_group_exp3_long <- pivot_longer(data = z_score_genus_top20_group_exp3,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_genus_top20_group_exp3_long <- as.data.frame(z_score_genus_top20_group_exp3_long)
z_score_genus_top20_group_exp3_long_merged <- merge(z_score_genus_top20_group_exp3_long, sig_aldex2_genus_group_exp3_result_top20, by.x = 4, by.y = 1, all.x = TRUE)

z_score_genus_top10_group_exp3_long <- pivot_longer(data = z_score_genus_top10_group_exp3,
                                                    cols = -c(Sample, SampleType, Experiment),
                                                    names_to = "ASV",
                                                    values_to = "Z")
z_score_genus_top10_group_exp3_long <- as.data.frame(z_score_genus_top10_group_exp3_long)
z_score_genus_top10_group_exp3_long_merged <- merge(z_score_genus_top10_group_exp3_long, sig_aldex2_genus_group_exp3_result_top10, by.x = 4, by.y = 1, all.x = TRUE)

################################################# Heatmaps for experiment 1 #################################################

#Heatmaps experiment 1
theme_set(theme_bw())
heatmap_phy_exp1 <- ggplot(z_score_phy_group_exp1_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp1

#Heatmaps
theme_set(theme_bw())
heatmap_phy_exp1_top10 <- ggplot(z_score_phy_top10_group_exp1_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp1_top10

theme_set(theme_bw())
heatmap_class_exp1 <- ggplot(z_score_class_group_exp1_long_merged, aes(x = Sample, y = Class, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_class_exp1

theme_set(theme_bw())
heatmap_order_top20_exp1 <- ggplot(z_score_order_top20_group_exp1_long_merged, aes(x = Sample, y = Order, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_order_top20_exp1

theme_set(theme_bw())
heatmap_fam_top20_exp1 <- ggplot(z_score_fam_top20_group_exp1_long_merged, aes(x = Sample, y = Family, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_fam_top20_exp1

theme_set(theme_bw())
heatmap_genus_top20_exp1 <- ggplot(z_score_genus_top20_group_exp1_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "Sample") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 9)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_top20_exp1

#Heatmaps
theme_set(theme_bw())
heatmap_genus_exp1_top10 <- ggplot(z_score_genus_top10_group_exp1_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_exp1_top10

Figure_heatmap_taxa_exp1 <- ggarrange(heatmap_phy_exp1,  heatmap_class_exp1, heatmap_order_top20_exp1, heatmap_fam_top20_exp1, heatmap_genus_top20_exp1,
                                      heights = c(13, 18, 20, 20, 20), ncol = 1, nrow = 5, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa_exp1
ggsave("results/16S_heatmap_FEPE_all_levels_SampleType_exp1.tiff", width = 10, height = 20, dpi = 150) # save graphic
ggsave("results/16S_heatmap_FEPE_all_levels_SampleType_exp1.pdf", width = 10, height = 20, dpi = 150) # save graphic
################################################# Heatmaps for experiment 2 #################################################
# Change facet SampleType labels
sample.type.labs <- c("EG", "FG", "FP24Hrs", "FP2Hrs", "SW")
names(sample.type.labs) <- c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water")

#z_score_phy_group_exp3_long_merged$SampleType <- factor(z_score_phy_group_exp3_long_merged$SampleType,
#                                                levels = c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water"))

# Heatmaps
theme_set(theme_bw())
heatmap_phy_exp2 <- ggplot(z_score_phy_group_exp2_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  #theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp2

#Heatmaps
theme_set(theme_bw())
heatmap_phy_exp2_top10 <- ggplot(z_score_phy_top10_group_exp2_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp2_top10

theme_set(theme_bw())
heatmap_class_exp2 <- ggplot(z_score_class_top20_group_exp2_long_merged, aes(x = Sample, y = Class, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_class_exp2

theme_set(theme_bw())
heatmap_order_top20_exp2 <- ggplot(z_score_order_top20_group_exp2_long_merged, aes(x = Sample, y = Order, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_order_top20_exp2

theme_set(theme_bw())
heatmap_fam_top20_exp2 <- ggplot(z_score_fam_top20_group_exp2_long_merged, aes(x = Sample, y = Family, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_fam_top20_exp2

theme_set(theme_bw())
heatmap_genus_top20_exp2 <- ggplot(z_score_genus_top20_group_exp2_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "Sample") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 9)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_top20_exp2

theme_set(theme_bw())
heatmap_genus_exp2_top10 <- ggplot(z_score_genus_top10_group_exp2_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_exp2_top10

Figure_heatmap_taxa_exp2 <- ggarrange(heatmap_phy_exp2_top10,  heatmap_class_exp2, heatmap_order_top20_exp2, heatmap_fam_top20_exp2, heatmap_genus_top20_exp2,
                                      heights = c(17, 20, 20, 20, 20), ncol = 1, nrow = 5, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa_exp2
ggsave("results/16S_heatmap_FEPE_all_levels_SampleType_exp2.tiff", width = 14, height = 21, dpi = 150) # save graphic
ggsave("results/16S_heatmap_FEPE_all_levels_SampleType_exp2.pdf", width = 14, height = 21, dpi = 150) # save graphic
################################################# Heatmaps for experiment 3 #################################################
# Change facet SampleType labels
sample.type.labs <- c("EG", "FG", "FP24Hrs", "FP2Hrs", "SW")
names(sample.type.labs) <- c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water")

# Heatmaps
theme_set(theme_bw())
heatmap_phy_exp3 <- ggplot(z_score_phy_group_exp3_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  #theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp3

theme_set(theme_bw())
heatmap_phy_exp3_top10 <- ggplot(z_score_phy_top10_group_exp3_long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy_exp3_top10

theme_set(theme_bw())
heatmap_class_exp3 <- ggplot(z_score_class_top20_group_exp3_long_merged, aes(x = Sample, y = Class, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_class_exp3

theme_set(theme_bw())
heatmap_order_top20_exp3 <- ggplot(z_score_order_top20_group_exp3_long_merged, aes(x = Sample, y = Order, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_order_top20_exp3

theme_set(theme_bw())
heatmap_fam_top20_exp3 <- ggplot(z_score_fam_top20_group_exp3_long_merged, aes(x = Sample, y = Family, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_fam_top20_exp3

theme_set(theme_bw())
heatmap_genus_top20_exp3 <- ggplot(z_score_genus_top20_group_exp3_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 9)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_top20_exp3

theme_set(theme_bw())
heatmap_genus_exp3_top10 <- ggplot(z_score_genus_top10_group_exp3_long_merged, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus_exp3_top10

Figure_heatmap_taxa_exp3 <- ggarrange(heatmap_phy_exp3_top10,  heatmap_class_exp3, heatmap_order_top20_exp3, heatmap_fam_top20_exp3, heatmap_genus_top20_exp3,
                                      heights = c(17, 20, 20, 20, 20), ncol = 1, nrow = 5, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa_exp3
ggsave("results/16S_heatmap_FEPE_all_levels_SampleType_exp3.pdf", width = 13, height = 21, dpi = 150) # save graphic

Figure_heatmap_taxa_all_exp_top20 <- ggarrange(heatmap_phy_exp1, heatmap_genus_top20_exp1, heatmap_phy_exp2,
                                         heatmap_genus_top20_exp2,  heatmap_phy_exp3, heatmap_genus_top20_exp3,labels = c("A", " ", "B", " ", "C", " "), 
                                      heights = c(13, 21, 19, 20, 20, 20), ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa_all_exp_top20
ggsave("results/16S_heatmap_phylum_genus_taxa_all_exp_top20.pdf", width = 8, height = 14, dpi = 150) # save graphic

Figure_heatmap_taxa_all_exp_top10 <- ggarrange(heatmap_phy_exp1_top10, heatmap_genus_exp1_top10, heatmap_phy_exp2_top10,
                                               heatmap_genus_exp2_top10,  heatmap_phy_exp3_top10, heatmap_genus_exp3_top10,labels = c("A", " ", "B", " ", "C", " "), 
                                               heights = c(10, 10, 10, 10, 10, 10), ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa_all_exp_top10
ggsave("results/16S_heatmap_phylum_genus_taxa_all_exp_top10.pdf", width = 7, height = 10, dpi = 200) # save graphic
