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
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("metagMisc")    # Miscellaneous functions for metagenomic analysis
library("vegan")
library("remotes")
library("mctoolsr")
library("RColorBrewer")
library("viridis")
library("decontam")     # Simple statistical identification of contaminated sequence features
library("plyr")
library("ALDEx2")
library("tibble")
library("gplots")
library("dplyr")        # filter and reformat data frames
library("ComplexHeatmap") # Package for Complex heatmaps reveal patterns and correlations in multidimensional genomic data
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

# Load sources
source("scr/miseqR.R")
source("scr/summary.R")

###################### Set plotting theme to blank ######################
theme_set(theme_bw())

###################### Import all input files ######################
# Import ASV table, taxonomy table, and sample table in excel format
ec_tab <- read_excel("data/EC_data.xlsx", sheet = "ec_matrix") # 47 observations (ECs) and 112 varibles (111 samples + EC IDs)
ec_tax <- read_excel("data/EC_data.xlsx", sheet = "ec_tax") # 47 observations (EC) and 7 variables (EC description)
ec_meta <- read_excel("data/EC_data.xlsx", sheet = "ec_meta") # 111 observations (samples) and 6 variables (sample factors)

ec_tab_all <- read_excel("data/EC_data_all.xlsx", sheet = "ec_matrix") # 1681 observations (ECs) and 112 varibles (111 samples + EC IDs)
ec_tax_all <- read_excel("data/EC_data_all.xlsx", sheet = "ec_tax") # 1681 observations (EC) and 5 variables (EC description)
ec_meta_all <- read_excel("data/EC_data_all.xlsx", sheet = "ec_meta") # 1681 observations (samples) and 6 variables (sample factors)

# Convert ASV and taxonomy tables to matrix format; sample table as data frame
ec_tab <- as.matrix(ec_tab)
class(ec_tab) # check object class/type
ec_tax <- as.matrix(ec_tax)
ec_meta <- as.data.frame(ec_meta)

ec_tab_all <- as.matrix(ec_tab_all)
class(ec_tab_all) # check object class/type
ec_tax_all<- as.matrix(ec_tax_all)
ec_meta_all <- as.data.frame(ec_meta_all)

# Rename row names using ASV column for ASV table and taxonomy table; for sample table rename using Sample column
# Takes all rows "," and uses first column as row names
rownames(ec_tab) <- ec_tab[,1]
rownames(ec_tax) <- ec_tax[,1]
rownames(ec_meta) <- ec_meta[,1]

rownames(ec_tab_all) <- ec_tab_all[,1]
rownames(ec_tax_all) <- ec_tax_all[,1]
rownames(ec_meta_all) <- ec_meta_all[,1]

# Exclude first column from all three tables; set ASV table as numeric
ec_tab <- ec_tab[,-1]
class(ec_tab) <- "numeric"
ec_tax <- ec_tax[,-1]
ec_meta <- ec_meta[,-1]
class(ec_meta)

ec_tab_all <- ec_tab_all[,-1]
class(ec_tab_all) <- "numeric"
ec_tax_all <- ec_tax_all[,-1]
ec_meta_all <- ec_meta_all[,-1]
class(ec_meta_all)

# Transform each matrix to a phyloseq component
EC = otu_table(ec_tab, taxa_are_rows = TRUE) # notes taxa (AVSs) are as rows
TAX = tax_table(ec_tax)
samples = sample_data(ec_meta)

EC_all = otu_table(ec_tab_all, taxa_are_rows = TRUE) # notes taxa (AVSs) are as rows
TAX_all = tax_table(ec_tax_all)
samples_all = sample_data(ec_meta_all)

# Create a phyloseq object by merging all three phyloseq components.
# Other components such as a phylogentic tree and representative sequences can also be added.
ec_phy <- phyloseq(EC, TAX, samples)
ec_phy # 47 taxa and 111 samples

ec_phy_all <- phyloseq(EC_all, TAX_all, samples_all)
ec_phy_all # 1681 taxa and 111 samples

# Remove low abundance samples if needed. Only two low abundance samples, both negative controls. These will be removed later.
# Check what is the minimum total abundance across all samples in the dataset
smin <- min(sample_sums(ec_phy)) # 27272 abundance
smin <- min(sample_sums(ec_phy_all)) # 1469305 abundance

# Ordination with nMDS without filtering
# All samples included
# All three experiments together
set.seed(1)
ec_phy_nmds <- ordinate(
  physeq = ec_phy,
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1487094

ec_phy_all_comp <- microbiome::transform(ec_phy_all, "compositional")
set.seed(1)

ec_phy_all_nmds <- ordinate(
  physeq = ec_phy_all,
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1346773

ec_phy_all_comp_nmds <- ordinate(
  physeq = ec_phy_all_comp,
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.1056146

# Set colors manually for the different sample types including controls
colors_all <- c("#1B9E77", "#D95F02", "#7570B3", "#D59D08", "#66A61E", "#E7298A", "#666666")

# Plot NMDS without any filtering
# All three experiments together
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = ec_phy_all,
  ordination = ec_phy_all_nmds,
  type = "samples",
  color = "SampleType",
  shape = factor("PcrReplicate")) +
  facet_nested(. ~ Experiment, scales = "free", labeller = labeller(Experiment = exp.labs.ec)) +
  scale_color_manual(values = colors_all,
                     name = "Sample",
                     labels = c("EG", "FP2Hrs", "FP24Hrs", "FG", "SW")) +
  scale_shape_manual(values = c(1, 8, 17, 3),
                     name = "PCR Replicate",
                     breaks = c("1", "2", "3"),
                     labels = c("PCR1", "PCR2", "PCR3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = SampleType, shape = factor(PcrReplicate)), size = 3) +#geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("nMDS2") + # add the title on y axis
  xlab("nMDS1") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # adjusts the title of the legend
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold"))  # remove the background of titles  
#annotate("text", x = 2, y = 2.5, label ="2D Stress: 0.18")
ggsave("results/ec_phy_nmds_all_samples_experiments.pdf", width = 8, height = 3, dpi = 150)

exp.labs.ec <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(exp.labs.ec) <- c("Experiment1", "Experiment2", "Experiment3")
  
theme_set(theme_bw())
phyloseq::plot_ordination(
    physeq = ec_phy_all_comp,
    ordination = ec_phy_all_comp_nmds,
    type = "samples",
    color = "SampleType",
    shape = factor("PcrReplicate")) +
    facet_nested(. ~ Experiment, scales = "free", labeller = labeller(Experiment = exp.labs.ec)) +
    scale_color_manual(values = colors_all,
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
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
    guides(fill = guide_legend(title = "Sample", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("nMDS2") + # add the title on y axis
    xlab("nMDS1") + # add the title on x axis
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
    theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold"))  # remove the background of titles
    #annotate("text", x = 2, y = 2.5, label ="2D Stress: 0.18")
    ggsave("results/ec_phy_comp_nmds_all_samples_experiments.pdf", width = 8, height = 3, dpi = 150)
  
# Import Aldex2 results and rename the X variable by OTU, all three experiments
aldex2_ec_result <- read_excel("data/EC_data.xlsx", sheet = "ec_matrix") # 47 observations (EC) and 7 variables (EC description)
aldex2_ec_result <- as.data.frame(aldex2_ec_result)
rownames(aldex2_ec_result) <- aldex2_ec_result[, 1]
aldex2_ec_result <- aldex2_ec_result[, -(1)]

# Adjusting zeros on the matrix, all three experiments
aldex2_ec_result_czm <- cmultRepl(t(aldex2_ec_result),  label=0, method="CZM")
aldex2_ec_result_czm_trans <- t(apply(aldex2_ec_result_czm, 1, function(x){log(x) - mean(log(x))}))

############ Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans_trav <- t(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans[, order(colnames(aldex2_ec_result_czm_trans))]
heatmap(aldex2_ec_result_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)


#Clean up presentation, all three experiments
ec_info <- data.frame(tax_table(ec_phy))
ec_info <- ec_info %>% rownames_to_column(var = "EC_number")

ec_count <- data.frame(otu_table(ec_phy))
ec_count <- ec_count %>% rownames_to_column(var = "EC_number")

ec_sample_tab <- data.frame(sample_data(ec_phy))
ec_sample_tab <- ec_sample_tab %>% rownames_to_column(var = "SampleID")

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
aldex2_ec_result_czm_trans_trav_group <- scale(aldex2_ec_result_czm_trans_trav)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top_ec = HeatmapAnnotation(border = TRUE,
                           Experiment = as.vector(ec_sample_tab$Experiment),
                           Sample = as.vector(ec_sample_tab$SampleType),
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
aldex2_ec_result_name <- as.data.frame(aldex2_ec_result)
str(aldex2_ec_result_name)
aldex2_ec_result_name <- rownames_to_column(aldex2_ec_result_name, var = "EC_number")
aldex2_ec_result_name <- left_join(aldex2_ec_result_name, ec_info)
head(aldex2_ec_result_name)
aldex2_ec_result_name <- aldex2_ec_result_name[, -(2:112)]

ec_count_total <- as.data.frame(ec_count)
ec_count_total$Total <- rowSums(ec_count_total[, -1])
head(ec_count_total)
ec_count_total <- ec_count_total[, -(2:112)]
head(ec_count_total)

z_score_aldex2_ec_result_name_count_total <- left_join(aldex2_ec_result_name, ec_count_total)
head(z_score_aldex2_ec_result_name_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 300, 3000, 30000, 300000),
                               c("#f5f5f5",
                                 "#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))
process_col_fun = colorRamp2(c("Carbon", "Nitrogen", "Sulfur"),
                             c("#af8dc3", "#f7f7f7", "#7fbf7b"))

process_col_discr = gpar(c("Carbon", "Nitrogen", "Sulfur"),
                         c("#af8dc3", "#f7f7f7", "#7fbf7b"))
colors = structure(6:8, names = c("Carbon", "Nitrogen", "Sulfur"))

ha_right_ec1 = rowAnnotation(
  Abundance = z_score_aldex2_ec_result_name_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))

ha_right_ec2 = rowAnnotation(
  Process = z_score_aldex2_ec_result_name_count_total$Element, border = TRUE, col  = list(Process = colors))
row_labels_ec = z_score_aldex2_ec_result_name_count_total$description

hm_ec_sel <- Heatmap(aldex2_ec_result_czm_trans_trav_group, name = "Z-score, CLR", col = col_matrix,
                     #column_title = "Experiments & Sample Types", 
                     #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     #cluster_columns = cluster_within_group(as.vector(ec_sample_tab$SampleType)),
                     column_split = as.vector(ec_sample_tab$Experiment),
                     border = TRUE,
                     top_annotation = ha_top_ec,
                     right_annotation = ha_right_ec1,
                     row_title = "Description",
                     row_labels = row_labels_ec,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     row_order = order(row_labels_ec),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     #show_row_names = TRUE,
                     show_heatmap_legend = TRUE) + ha_right_ec2
hm_ec_sel

ggsave("results/16S_heatmap_phy_SampleType_Experiment_new.pdf", width = 12, height = 6, dpi = 150) # save graphic
