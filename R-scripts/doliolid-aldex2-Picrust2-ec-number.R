# Set working directory

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
library("DataCombine")
library("FD")
library("circlize")

# Load tables and make adjustments on each one
fepe_exp1_ec_number <- read.table("picrust/EC_exp1_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", row.names = 1, check.names=F, comment.char="")
fepe_exp2_ec_number <- read.table("picrust/EC_exp2_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", row.names = 1, check.names=F, comment.char="")
fepe_exp3_ec_number <- read.table("picrust/EC_exp3_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", row.names = 1, check.names=F, comment.char="")
fepe_all_ec_number <- read.table("picrust/EC_all_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", row.names = 1, check.names=F, comment.char="")

fepe_exp1_ec_number <- fepe_exp1_ec_number[, -(1)]
fepe_exp2_ec_number <- fepe_exp2_ec_number[, -(1)]
fepe_exp3_ec_number <- fepe_exp3_ec_number[, -(1)]
fepe_all_ec_number <- fepe_all_ec_number[, -(1)]

# Create tables with gene function ids and description only
fepe_exp1_ec_number_descrip <- read.table("picrust/EC_exp1_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_exp2_ec_number_descrip <- read.table("picrust/EC_exp2_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_exp3_ec_number_descrip <- read.table("picrust/EC_exp3_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_all_ec_number_descrip <- read.table("picrust/EC_all_pred_metagenome_unstrat_descrip.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")

fepe_exp1_ec_number_descrip <- fepe_exp1_ec_number_descrip[, -(3:26)]
fepe_exp2_ec_number_descrip <- fepe_exp2_ec_number_descrip[, -(3:47)]
fepe_exp3_ec_number_descrip <- fepe_exp3_ec_number_descrip[, -(3:44)]
fepe_all_ec_number_descrip <- fepe_all_ec_number_descrip [, -(3:113)]

# Load metadata files for each experiment
fepe_exp1_metadata <- read.table("picrust/metadata_exp1_picrust.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_exp2_metadata <- read.table("picrust/metadata_exp2_picrust.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_exp3_metadata <- read.table("picrust/metadata_exp3_picrust.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
fepe_all_metadata <- read.table("picrust/metadata_all_picrust.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")

# Create a vector by using originals metadata files
conds_sample_type_exp1 <- c(fepe_exp1_metadata$SampleType)
conds_sample_type_exp2 <- c(fepe_exp2_metadata$SampleType)
conds_sample_type_exp3 <- c(fepe_exp3_metadata$SampleType)
conds_sample_type_all <- c(fepe_all_metadata$SampleType)

# Perform aldex2 analysis using the transverse matrix and conds (factors)
# Aldex2 performs the entire analysis, i.e., calculates the clr, p-values, and FDR for each sample group.
# Write Aldex2 results to file as csv
exp1_ec_number_aldex <- aldex(fepe_exp1_ec_number, conds_sample_type_exp1, mc.samples=128, test="kw", effect=TRUE,
                                   include.sample.summary=FALSE, denom="all", verbose=FALSE)
exp1_ec_number_aldex_res <- rownames_to_column(exp1_ec_number_aldex, var ="function")
#write.csv(exp1_ec_number_aldex_res, "/home/tjp99569/fepe/picrust2/aldex/exp1_ec_number_aldex_res.csv")

exp2_ec_number_aldex <- aldex(fepe_exp2_ec_number, conds_sample_type_exp2, mc.samples=128, test="kw", effect=TRUE,
                                   include.sample.summary=FALSE, denom="all", verbose=FALSE)
exp2_ec_number_aldex_res <- rownames_to_column(exp2_ec_number_aldex, var ="function")
#write.csv(exp2_ec_number_aldex_res, "/home/tjp99569/fepe/picrust2/aldex/exp2_ec_number_aldex_res.csv")

exp3_ec_number_aldex <- aldex(fepe_exp3_ec_number, conds_sample_type_exp3, mc.samples=128, test="kw", effect=TRUE,
                                   include.sample.summary=FALSE, denom="all", verbose=FALSE)
exp3_ec_number_aldex_res <- rownames_to_column(exp3_ec_number_aldex, var ="function")
#write.csv(exp3_ec_number_aldex_res, "/home/tjp99569/fepe/picrust2/aldex/exp3_ec_number_aldex_res.csv")

all_ec_number_aldex <- aldex(fepe_all_ec_number, conds_sample_type_exp3, mc.samples=128, test="kw", effect=TRUE,
                              include.sample.summary=FALSE, denom="all", verbose=FALSE)
all_ec_number_aldex_res <- rownames_to_column(all_ec_number_aldex, var ="function")
#write.csv(all_ec_number_aldex_res, "/home/tjp99569/fepe/picrust2/aldex/all_ec_number_aldex_res.csv")

# Filter Aldex2 results by sig kw.ep and join the function description
exp1_ec_number_aldex_res <- read.csv("picrust/exp1_ec_number_aldex_res.csv")
exp2_ec_number_aldex_res <- read.csv("picrust/exp2_ec_number_aldex_res.csv")
exp3_ec_number_aldex_res <- read.csv("picrust/exp3_ec_number_aldex_res.csv")
all_ec_number_aldex_res <- read.csv("picrust/all_ec_number_aldex_res.csv")

exp1_ec_number_aldex_sig <- exp1_ec_number_aldex_res %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(EC_number, kw.ep, kw.eBH)
exp1_ec_number_aldex_sig <- left_join(exp1_ec_number_aldex_sig, fepe_exp1_ec_number_descrip)
write.csv(exp1_ec_number_aldex_sig, "picrust/exp1_ec_number_aldex_sig.csv")
exp1_ec_number_aldex_sig_top20 <- exp1_ec_number_aldex_sig %>% top_n(-20, kw.ep)
exp1_ec_number_aldex_sig_top20 <- left_join(exp1_ec_number_aldex_sig_top20, fepe_exp1_ec_number_descrip)

#exp2_ec_number_aldex_sig <- exp2_ec_number_aldex_res %>%
exp2_ec_number_aldex_sig <- exp2_ec_number_aldex_res %>%
    filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(EC_number, kw.ep, kw.eBH)
exp2_ec_number_aldex_sig <- left_join(exp2_ec_number_aldex_sig, fepe_exp2_ec_number_descrip)
write.csv(exp2_ec_number_aldex_sig, "picrust/exp2_ec_number_aldex_sig.csv")
exp2_ec_number_aldex_sig_top20 <- exp2_ec_number_aldex_sig %>% top_n(-20, kw.ep)
exp2_ec_number_aldex_sig_top20 <- left_join(exp2_ec_number_aldex_sig_top20, fepe_exp2_ec_number_descrip)

#exp3_ec_number_aldex_sig <- exp3_ec_number_aldex_res %>%
exp3_ec_number_aldex_sig <- exp3_ec_number_aldex_res %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(EC_number, kw.ep, kw.eBH)
exp3_ec_number_aldex_sig <- left_join(exp3_ec_number_aldex_sig, fepe_exp3_ec_number_descrip)
write.csv(exp3_ec_number_aldex_sig, "picrust/exp3_ec_number_aldex_sig.csv")
exp3_ec_number_aldex_sig_top20 <- exp3_ec_number_aldex_sig %>% top_n(-20, kw.ep)
exp3_ec_number_aldex_sig_top20 <- left_join(exp3_ec_number_aldex_sig_top20, fepe_exp3_ec_number_descrip)

# All three experiments
all_ec_number_aldex_sig <- all_ec_number_aldex_res %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(EC_number, kw.ep, kw.eBH)
all_ec_number_aldex_sig <- left_join(all_ec_number_aldex_sig, fepe_all_ec_number_descrip)
write.csv(all_ec_number_aldex_sig, "picrust/all_ec_number_aldex_sig.csv")

all_ec_number_aldex_sig_top20 <- all_ec_number_aldex_sig %>% top_n(-20, kw.ep)
all_ec_number_aldex_sig_top20 <- left_join(all_ec_number_aldex_sig_top20, fepe_all_ec_number_descrip)
all_ec_number_aldex_sig_top100 <- all_ec_number_aldex_sig %>% top_n(-100, kw.ep)
all_ec_number_aldex_sig_top100 <- left_join(all_ec_number_aldex_sig_top100, fepe_all_ec_number_descrip)

# Create clr objects by using EC_number ids and original otu table, all significant OTUs and only top 20, experiment 1
exp1_ec_number_aldex_sig_count <- rownames_to_column(fepe_exp1_ec_number, var ="EC_number")
exp1_ec_number_aldex_sig_count <- left_join(exp1_ec_number_aldex_sig, exp1_ec_number_aldex_sig_count)
write.csv(exp1_ec_number_aldex_sig_count, "picrust/exp1_ec_number_aldex_sig_count.csv")

exp1_ec_number_aldex_clr <- exp1_ec_number_aldex_sig_count[, -(2:4)]
rownames(exp1_ec_number_aldex_clr) <- exp1_ec_number_aldex_clr$EC_number
exp1_ec_number_aldex_clr <- exp1_ec_number_aldex_clr[, -1]
class(exp1_ec_number_aldex_clr)

exp1_ec_number_aldex_sig_top20_count <- rownames_to_column(fepe_exp1_ec_number, var ="EC_number")
exp1_ec_number_aldex_sig_top20_count <- left_join(exp1_ec_number_aldex_sig_top20, exp1_ec_number_aldex_sig_top20_count)
exp1_ec_number_aldex_clr_top20 <- exp1_ec_number_aldex_sig_top20_count[, -(2:4)]
rownames(exp1_ec_number_aldex_clr_top20) <- exp1_ec_number_aldex_clr_top20$EC_number
exp1_ec_number_aldex_clr_top20 <- exp1_ec_number_aldex_clr_top20[, -1]

# Create clr objects by using EC_number ids and original otu table, all significant OTUs and only top 20, experiment 2
exp2_ec_number_aldex_sig_count <- rownames_to_column(fepe_exp2_ec_number, var ="EC_number")
exp2_ec_number_aldex_sig_count <- left_join(exp2_ec_number_aldex_sig, exp2_ec_number_aldex_sig_count)
write.csv(exp2_ec_number_aldex_sig_count, "picrust/exp2_ec_number_aldex_sig_count.csv")
exp2_ec_number_aldex_clr <- exp2_ec_number_aldex_sig_count[, -(2:4)]
rownames(exp2_ec_number_aldex_clr) <- exp2_ec_number_aldex_clr$EC_number
exp2_ec_number_aldex_clr <- exp2_ec_number_aldex_clr[, -1]

exp2_ec_number_aldex_sig_top20_count <- rownames_to_column(fepe_exp2_ec_number, var ="EC_number")
exp2_ec_number_aldex_sig_top20_count <- left_join(exp2_ec_number_aldex_sig_top20, exp2_ec_number_aldex_sig_top20_count)
exp2_ec_number_aldex_clr_top20 <- exp2_ec_number_aldex_sig_top20_count[, -(2:4)]
rownames(exp2_ec_number_aldex_clr_top20) <- exp2_ec_number_aldex_clr_top20$EC_number
exp2_ec_number_aldex_clr_top20 <- exp2_ec_number_aldex_clr_top20[, -1]

# Create clr objects by using EC_number ids and original otu table, all significant OTUs and only top 20, experiment 3
exp3_ec_number_aldex_sig_count <- rownames_to_column(fepe_exp3_ec_number, var ="EC_number")
exp3_ec_number_aldex_sig_count <- left_join(exp3_ec_number_aldex_sig, exp3_ec_number_aldex_sig_count)
write.csv(exp3_ec_number_aldex_sig_count, "picrust/exp3_ec_number_aldex_sig_count.csv")
exp3_ec_number_aldex_clr <- exp3_ec_number_aldex_sig_count[, -(2:4)]
rownames(exp3_ec_number_aldex_clr) <- exp3_ec_number_aldex_clr$EC_number
exp3_ec_number_aldex_clr <- exp3_ec_number_aldex_clr[, -1]

exp3_ec_number_aldex_sig_top20_count <- rownames_to_column(fepe_exp3_ec_number, var ="EC_number")
exp3_ec_number_aldex_sig_top20_count <- left_join(exp3_ec_number_aldex_sig_top20, exp3_ec_number_aldex_sig_top20_count)
exp3_ec_number_aldex_clr_top20 <- exp3_ec_number_aldex_sig_top20_count[, -(2:4)]
rownames(exp3_ec_number_aldex_clr_top20) <- exp3_ec_number_aldex_clr_top20$EC_number
exp3_ec_number_aldex_clr_top20 <- exp3_ec_number_aldex_clr_top20[, -1]

# Create clr objects by using EC_number ids and original otu table, all significant OTUs and only top 20, experiment 3
all_ec_number_aldex_sig_count <- rownames_to_column(fepe_all_ec_number, var ="EC_number")
all_ec_number_aldex_sig_count <- left_join(all_ec_number_aldex_sig, all_ec_number_aldex_sig_count)
write.csv(all_ec_number_aldex_sig_count, "picrust/all_ec_number_aldex_sig_count.csv")
all_ec_number_aldex_clr <- all_ec_number_aldex_sig_count[, -(2:4)]
rownames(all_ec_number_aldex_clr) <- all_ec_number_aldex_clr$EC_number
all_ec_number_aldex_clr <- all_ec_number_aldex_clr[, -1]

all_ec_number_aldex_sig_top20_count <- rownames_to_column(fepe_all_ec_number, var ="EC_number")
all_ec_number_aldex_sig_top20_count <- left_join(all_ec_number_aldex_sig_top20, all_ec_number_aldex_sig_top20_count)
all_ec_number_aldex_clr_top20 <- all_ec_number_aldex_sig_top20_count[, -(2:4)]
rownames(all_ec_number_aldex_clr_top20) <- all_ec_number_aldex_clr_top20$EC_number
all_ec_number_aldex_clr_top20 <- all_ec_number_aldex_clr_top20[, -1]

all_ec_number_aldex_sig_top100_count <- rownames_to_column(fepe_all_ec_number, var ="EC_number")
all_ec_number_aldex_sig_top100_count <- left_join(all_ec_number_aldex_sig_top100, all_ec_number_aldex_sig_top100_count)
all_ec_number_aldex_clr_top100 <- all_ec_number_aldex_sig_top100_count[, -(2:4)]
rownames(all_ec_number_aldex_clr_top100) <- all_ec_number_aldex_clr_top100$EC_number
all_ec_number_aldex_clr_top100 <- all_ec_number_aldex_clr_top100[, -1]

# Adjusting zeros on the matrix, experiment 1
sum(exp1_ec_number_aldex_clr == 0)
exp1_ec_number_aldex_clr_czm <- cmultRepl(t(exp1_ec_number_aldex_clr),  label=0, method="CZM")
exp1_ec_number_aldex_clr_trans <- t(apply(exp1_ec_number_aldex_clr_czm, 1, function(x){log(x) - mean(log(x))}))
#exp1_ec_number_aldex_clr_trans <- t(exp1_ec_number_aldex_clr_trans)

sum(exp1_ec_number_aldex_clr_top20 == 0) # no need to use cmultRepl function
exp1_ec_number_aldex_clr_top20_trans <- apply(exp1_ec_number_aldex_clr_top20, 1, function(x){log(x) - mean(log(x))})

# Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(exp1_ec_number_aldex_clr_top20_trans)
exp1_ec_number_aldex_clr_top20_trans_trav <- t(exp1_ec_number_aldex_clr_top20_trans)
exp1_ec_number_aldex_clr_top20_trans[, order(colnames(exp1_ec_number_aldex_clr_top20_trans))]
heatmap(exp1_ec_number_aldex_clr_top20_trans_trav, scale = "none", col = bluered(100))

head(all_ec_number_aldex_clr_top20_trans)
all_ec_number_aldex_clr_top20_trans_trav <- t(all_ec_number_aldex_clr_top20_trans)
all_ec_number_aldex_clr_top20_trans[, order(colnames(all_ec_number_aldex_clr_top20_trans))]
heatmap(all_ec_number_aldex_clr_top20_trans_trav, scale = "none", col = bluered(100))

sum(exp2_ec_number_aldex_clr_top20 == 0) # no need to use cmultRepl function
exp2_ec_number_aldex_clr_czm <- cmultRepl(t(exp2_ec_number_aldex_clr),  label=0, method="CZM")
exp2_ec_number_aldex_clr_trans <- t(apply(exp2_ec_number_aldex_clr_czm, 1, function(x){log(x) - mean(log(x))}))
#exp2_ec_number_aldex_clr_trans <- t(exp2_ec_number_aldex_clr_trans)

sum(exp2_ec_number_aldex_clr_top20 == 0)
exp2_ec_number_aldex_clr_top20_czm <- cmultRepl(t(exp2_ec_number_aldex_clr_top20),  label=0, method="CZM")
exp2_ec_number_aldex_clr_top20_trans <- t(apply(exp2_ec_number_aldex_clr_top20_czm, 1, function(x){log(x) - mean(log(x))}))
#exp2_ec_number_aldex_clr_top20_trans <- t(exp2_ec_number_aldex_clr_top20_trans)

sum(exp3_ec_number_aldex_clr == 0)
exp3_ec_number_aldex_clr_czm <- cmultRepl(t(exp3_ec_number_aldex_clr),  label=0, method="CZM")
exp3_ec_number_aldex_clr_trans <- t(apply(exp3_ec_number_aldex_clr_czm, 1, function(x){log(x) - mean(log(x))}))
#exp3_ec_number_aldex_clr_trans <- t(exp3_ec_number_aldex_clr_trans)

sum(exp3_ec_number_aldex_clr_top20 == 0) # no need to use cmultRepl function
exp3_ec_number_aldex_clr_top20_trans <- apply(exp3_ec_number_aldex_clr_top20, 1, function(x){log(x) - mean(log(x))})

sum(all_ec_number_aldex_clr == 0)
all_ec_number_aldex_clr_czm <- cmultRepl(t(all_ec_number_aldex_clr),  label=0, method="CZM")
all_ec_number_aldex_clr_trans <- t(apply(all_ec_number_aldex_clr_czm, 1, function(x){log(x) - mean(log(x))}))
#all_ec_number_aldex_clr_trans <- t(all_ec_number_aldex_clr_trans)

sum(all_ec_number_aldex_clr_top20 == 0) 
all_ec_number_aldex_clr_top20_czm <- cmultRepl(t(all_ec_number_aldex_clr_top20),  label=0, method="CZM")
all_ec_number_aldex_clr_top20_trans <- t(apply(all_ec_number_aldex_clr_top20_czm, 1, function(x){log(x) - mean(log(x))}))
#all_ec_number_aldex_clr_top20_trans <- t(all_ec_number_aldex_clr_top20_trans)

sum(all_ec_number_aldex_clr_top100 == 0) 
all_ec_number_aldex_clr_top100_czm <- cmultRepl(t(all_ec_number_aldex_clr_top100),  label=0, method="CZM")
all_ec_number_aldex_clr_top100_trans <- t(apply(all_ec_number_aldex_clr_top100_czm, 1, function(x){log(x) - mean(log(x))}))
#all_ec_number_aldex_clr_top100_trans <- t(all_ec_number_aldex_clr_top100_trans)

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
z_exp1_ec_number_aldex_clr_top20_trans_trav_group <- scale(exp1_ec_number_aldex_clr_top20_trans_trav)
z_all_ec_number_aldex_clr_top20_trans_trav_group <- scale(all_ec_number_aldex_clr_top20_trans_trav)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top_exp1 = HeatmapAnnotation(border = TRUE,
                           Experiment = as.vector(z_score_exp1_ec_number_factors$Experiment),
                           Sample = as.vector(z_score_exp1_ec_number_factors$SampleType),
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

ha_top_all = HeatmapAnnotation(border = TRUE,
                           Experiment = as.vector(z_score_all_ec_number_factors$Experiment),
                           Sample = as.vector(z_score_all_ec_number_factors$SampleType),
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
z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name <- as.data.frame(z_exp1_ec_number_aldex_clr_top20_trans_trav_group)
str(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name)
z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name <- rownames_to_column(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name, var = "EC_number")
z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name <- left_join(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name, fepe_exp1_ec_number_descrip)
head(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name)
z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name <- z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name[, -(2:25)]

z_all_ec_number_aldex_clr_top20_trans_trav_group_name <- as.data.frame(z_all_ec_number_aldex_clr_top20_trans_trav_group)
str(z_all_ec_number_aldex_clr_top20_trans_trav_group_name)
z_all_ec_number_aldex_clr_top20_trans_trav_group_name <- rownames_to_column(z_all_ec_number_aldex_clr_top20_trans_trav_group_name, var = "EC_number")
z_all_ec_number_aldex_clr_top20_trans_trav_group_name <- left_join(z_all_ec_number_aldex_clr_top20_trans_trav_group_name, fepe_all_ec_number_descrip)
head(z_all_ec_number_aldex_clr_top20_trans_trav_group_name)
z_all_ec_number_aldex_clr_top20_trans_trav_group_name <- z_all_ec_number_aldex_clr_top20_trans_trav_group_name[, -(2:25)]


ec_count_total_exp1 <- as.data.frame(fepe_exp1_ec_number)
ec_count_total_exp1 <- rownames_to_column(ec_count_total_exp1, var = "EC_number")
ec_count_total_exp1$Total <- rowSums(ec_count_total_exp1[, -1])
head(ec_count_total_exp1)
ec_count_total_exp1 <- ec_count_total_exp1[, -(2:25)]

ec_count_total_all <- as.data.frame(fepe_all_ec_number)
ec_count_total_all <- rownames_to_column(ec_count_total_all, var = "EC_number")
ec_count_total_all$Total <- rowSums(ec_count_total_all[, -1])
head(ec_count_total_all)
ec_count_total_all <- ec_count_total_all[, -(2:112)]


z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name_total <- left_join(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name, ec_count_total_exp1)
head(z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name_total)

z_all_ec_number_aldex_clr_top20_trans_trav_group_name_total <- left_join(z_all_ec_number_aldex_clr_top20_trans_trav_group_name, ec_count_total_all)
head(z_all_ec_number_aldex_clr_top20_trans_trav_group_name_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 300, 3000, 30000, 300000),
                               c("#f5f5f5",
                                 "#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))
ha_right_exp1 = rowAnnotation(
  Abundance = z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_exp1 = z_exp1_ec_number_aldex_clr_top20_trans_trav_group_name_total$description

ha_right_all = rowAnnotation(
  Abundance = z_all_ec_number_aldex_clr_top20_trans_trav_group_name_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))
row_labels_all = z_all_ec_number_aldex_clr_top20_trans_trav_group_name_total$description


# Plot heatmap at the phylum level
hm_z_exp1_ec_number_aldex_clr_top20 <- Heatmap(z_exp1_ec_number_aldex_clr_top20_trans_trav_group, name = "Z-score, CLR", col = col_matrix3,
                     #column_title = "Experiments & Sample Types", 
                     #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     #cluster_columns = cluster_within_group(as.vector(z_score_exp1_ec_number_factors$SampleType)),
                     column_split = as.vector(z_score_exp1_ec_number_factors$Experiment),
                     border = TRUE,
                     top_annotation = ha_top_exp1,
                     right_annotation = ha_right_exp1,
                     row_title = "EC description",
                     row_labels = row_labels_exp1,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     row_order = order(row_labels_exp1),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_z_exp1_ec_number_aldex_clr_top20

hm_z_all_ec_number_aldex_clr_top20 <- Heatmap(z_all_ec_number_aldex_clr_top20_trans_trav_group, name = "Z-score, CLR", col = col_matrix3,
                                               #column_title = "Experiments & Sample Types", 
                                               #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                                               #cluster_columns = cluster_within_group(as.vector(z_score_all_ec_number_factors$SampleType)),
                                               column_split = as.vector(z_score_all_ec_number_factors$Experiment),
                                               border = TRUE,
                                               top_annotation = ha_top_all,
                                               right_annotation = ha_right_all,
                                               row_title = "EC description",
                                               row_labels = row_labels_all,
                                               row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                               row_names_gp = gpar(fontsize = 8),
                                               column_names_gp = gpar(fontsize = 8),
                                               row_order = order(row_labels_all),
                                               rect_gp = gpar(col = "white", lwd = 1),
                                               show_column_names = FALSE,
                                               show_heatmap_legend = TRUE)
hm_z_all_ec_number_aldex_clr_top20

# Apply scale to calculate Z-score on the clr_transf zero adjusted objects, then transform from matrix to data.frame, experiment 1
z_score_exp1_ec_number <- as.data.frame(scale(exp1_ec_number_aldex_clr_trans))
z_score_exp2_ec_number <- as.data.frame(scale(exp2_ec_number_aldex_clr_trans))
z_score_exp3_ec_number <- as.data.frame(scale(exp3_ec_number_aldex_clr_trans))
z_score_all_ec_number <- as.data.frame(scale(all_ec_number_aldex_clr_trans))

z_score_exp1_ec_number_top20 <- as.data.frame(scale(exp1_ec_number_aldex_clr_top20_trans))
z_score_exp2_ec_number_top20 <- as.data.frame(scale(exp2_ec_number_aldex_clr_top20_trans))
z_score_exp3_ec_number_top20 <- as.data.frame(scale(exp3_ec_number_aldex_clr_top20_trans))
z_score_all_ec_number_top20 <- as.data.frame(scale(all_ec_number_aldex_clr_top20_trans))
z_score_all_ec_number_top100 <- as.data.frame(scale(all_ec_number_aldex_clr_top100_trans))

# Use rownames to create another variable called SampleID, experiment 1
z_score_exp1_ec_number_group <- tibble::rownames_to_column(z_score_exp1_ec_number, "SampleID")
z_score_exp2_ec_number_group <- tibble::rownames_to_column(z_score_exp2_ec_number, "SampleID")
z_score_exp3_ec_number_group <- tibble::rownames_to_column(z_score_exp3_ec_number, "SampleID")
z_score_all_ec_number_group <- tibble::rownames_to_column(z_score_all_ec_number, "SampleID")

z_score_exp1_ec_number_group_top20 <- tibble::rownames_to_column(z_score_exp1_ec_number_top20, "SampleID")
z_score_exp2_ec_number_group_top20 <- tibble::rownames_to_column(z_score_exp2_ec_number_top20, "SampleID")
z_score_exp3_ec_number_group_top20 <- tibble::rownames_to_column(z_score_exp3_ec_number_top20, "SampleID")
z_score_all_ec_number_group_top20 <- tibble::rownames_to_column(z_score_all_ec_number_top20, "SampleID")
z_score_all_ec_number_group_top100 <- tibble::rownames_to_column(z_score_all_ec_number_top100, "SampleID")

# Add sample factors (Sample, SampleType, Experiment) to the Z.score object, experiment 1
z_score_exp1_ec_number_factors <- fepe_exp1_metadata[, -(5:18)]
z_score_exp1_ec_number_factors <- z_score_exp1_ec_number_factors[, -(1)]
z_score_exp1_ec_number_group <- data.frame(cbind(z_score_exp1_ec_number_factors, z_score_exp1_ec_number_group))
z_score_exp1_ec_number_group_top20 <- data.frame(cbind(z_score_exp1_ec_number_factors, z_score_exp1_ec_number_group_top20))

z_score_exp2_ec_number_factors <- fepe_exp2_metadata[, -(5:18)]
z_score_exp2_ec_number_factors <- z_score_exp2_ec_number_factors[, -(1)]
z_score_exp2_ec_number_group <- data.frame(cbind(z_score_exp2_ec_number_factors, z_score_exp2_ec_number_group))
z_score_exp2_ec_number_group_top20 <- data.frame(cbind(z_score_exp2_ec_number_factors, z_score_exp2_ec_number_group_top20))

z_score_exp3_ec_number_factors <- fepe_exp3_metadata[, -(5:18)]
z_score_exp3_ec_number_factors <- z_score_exp3_ec_number_factors[, -(1)]
z_score_exp3_ec_number_group <- data.frame(cbind(z_score_exp3_ec_number_factors, z_score_exp3_ec_number_group))
z_score_exp3_ec_number_group_top20 <- data.frame(cbind(z_score_exp3_ec_number_factors, z_score_exp3_ec_number_group_top20))

z_score_all_ec_number_factors <- fepe_all_metadata[, -(5:18)]
z_score_all_ec_number_factors <- z_score_all_ec_number_factors[, -(1)]
z_score_all_ec_number_group <- data.frame(cbind(z_score_all_ec_number_factors, z_score_all_ec_number_group))
z_score_all_ec_number_group_top20 <- data.frame(cbind(z_score_all_ec_number_factors, z_score_all_ec_number_group_top20))
z_score_all_ec_number_group_top100 <- data.frame(cbind(z_score_all_ec_number_factors, z_score_all_ec_number_group_top100))

#Add sample factors (SampleType) to the Z.score object, experiment 1
z_score_all_ec_number_group_long <- pivot_longer(data = z_score_all_ec_number_group,
                                                  cols = -c(Sample, Experiment, SampleType, SampleID),
                                                  names_to = "EC_number",
                                                  values_to = "Z")
z_score_all_ec_number_group_long <- as.data.frame(z_score_all_ec_number_group_long)

z_score_all_ec_number_group_long_top20 <- pivot_longer(data = z_score_all_ec_number_group_top20,
                                                        cols = -c(Sample, Experiment, SampleType, SampleID),
                                                        names_to = "EC_number",
                                                        values_to = "Z")
z_score_all_ec_number_group_long_top20 <- as.data.frame(z_score_all_ec_number_group_long_top20)

z_score_all_ec_number_group_long_top100 <- pivot_longer(data = z_score_all_ec_number_group_top100,
                                                       cols = -c(Sample, Experiment, SampleType, SampleID),
                                                       names_to = "EC_number",
                                                       values_to = "Z")
z_score_all_ec_number_group_long_top100 <- as.data.frame(z_score_all_ec_number_group_long_top100)

z_score_exp1_ec_number_group_long <- pivot_longer(data = z_score_exp1_ec_number_group,
                                           cols = -c(Sample, Experiment, SampleType, SampleID),
                                           names_to = "EC_number",
                                           values_to = "Z")
z_score_exp1_ec_number_group_long <- as.data.frame(z_score_exp1_ec_number_group_long)

z_score_exp1_ec_number_group_long_top20 <- pivot_longer(data = z_score_exp1_ec_number_group_top20,
                                                  cols = -c(Sample, Experiment, SampleType, SampleID),
                                                  names_to = "EC_number",
                                                  values_to = "Z")
z_score_exp1_ec_number_group_long_top20 <- as.data.frame(z_score_exp1_ec_number_group_long_top20)


z_score_exp2_ec_number_group_long <- pivot_longer(data = z_score_exp2_ec_number_group,
                                                  cols = -c(Sample, SampleType),
                                                  names_to = "EC_number",
                                                  values_to = "Z")
z_score_exp2_ec_number_group_long <- as.data.frame(z_score_exp2_ec_number_group_long)

z_score_exp2_ec_number_group_long_top20 <- pivot_longer(data = z_score_exp2_ec_number_group_top20,
                                                        cols = -c(Sample, SampleType),
                                                        names_to = "EC_number",
                                                        values_to = "Z")
z_score_exp2_ec_number_group_long_top20 <- as.data.frame(z_score_exp2_ec_number_group_long_top20)

z_score_exp3_ec_number_group_long <- pivot_longer(data = z_score_exp3_ec_number_group,
                                                  cols = -c(Sample, SampleType),
                                                  names_to = "EC_number",
                                                  values_to = "Z")
z_score_exp3_ec_number_group_long <- as.data.frame(z_score_exp3_ec_number_group_long)

z_score_exp3_ec_number_group_long_top20 <- pivot_longer(data = z_score_exp3_ec_number_group_top20,
                                                        cols = -c(Sample, SampleType),
                                                        names_to = "EC_number",
                                                        values_to = "Z")
z_score_exp3_ec_number_group_long_top20 <- as.data.frame(z_score_exp3_ec_number_group_long_top20)


# Replace "." by "-" so it matches the description of each EC number
z_score_all_ec_number_group_long$EC_number <- gsub("\\EC.", "EC:", z_score_all_ec_number_group_long$EC_number)
z_score_all_ec_number_group_long_merged <- merge(z_score_all_ec_number_group_long, all_ec_number_aldex_sig, by.x = 5, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_all_ec_number_group_long_top20$EC_number <- gsub("\\EC.", "EC:", z_score_all_ec_number_group_long_top20$EC_number)
z_score_all_ec_number_group_long_merged_top20 <- merge(z_score_all_ec_number_group_long_top20, all_ec_number_aldex_sig, by.x = 5, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_all_ec_number_group_long_top100$EC_number <- gsub("\\EC.", "EC:", z_score_all_ec_number_group_long_top100$EC_number)
z_score_all_ec_number_group_long_merged_top100 <- merge(z_score_all_ec_number_group_long_top100, all_ec_number_aldex_sig, by.x = 5, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp1_ec_number_group_long$EC_number <- gsub("\\EC.", "EC:", z_score_exp1_ec_number_group_long$EC_number)
z_score_exp1_ec_number_group_long_merged <- merge(z_score_exp1_ec_number_group_long, exp1_ec_number_aldex_sig, by.x = 5, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp1_ec_number_group_long_top20$EC_number <- gsub("\\EC.", "EC:", z_score_exp1_ec_number_group_long_top20$EC_number)
z_score_exp1_ec_number_group_long_merged_top20 <- merge(z_score_exp1_ec_number_group_long_top20, exp1_ec_number_aldex_sig, by.x = 5, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp2_ec_number_group_long$EC_number <- gsub("\\EC.", "EC:", z_score_exp2_ec_number_group_long$EC_number)
z_score_exp2_ec_number_group_long_merged <- merge(z_score_exp2_ec_number_group_long, exp2_ec_number_aldex_sig, by.x = 3, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp2_ec_number_group_long_top20$EC_number <- gsub("\\EC.", "EC:", z_score_exp2_ec_number_group_long_top20$EC_number)
z_score_exp2_ec_number_group_long_merged_top20 <- merge(z_score_exp2_ec_number_group_long_top20, exp2_ec_number_aldex_sig, by.x = 3, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp3_ec_number_group_long$EC_number <- gsub("\\EC.", "EC:", z_score_exp3_ec_number_group_long$EC_number)
z_score_exp3_ec_number_group_long_merged <- merge(z_score_exp3_ec_number_group_long, exp3_ec_number_aldex_sig, by.x = 3, by.y = 1, all.x = TRUE, sort = FALSE)

z_score_exp3_ec_number_group_long_top20$EC_number <- gsub("\\EC.", "EC:", z_score_exp3_ec_number_group_long_top20$EC_number)
z_score_exp3_ec_number_group_long_merged_top20 <- merge(z_score_exp3_ec_number_group_long_top20, exp3_ec_number_aldex_sig, by.x = 3, by.y = 1, all.x = TRUE, sort = FALSE)

################################################# Heatmaps for experiment 1 #################################################
# Change facet SampleType labels for plotting on heatmaps
sample.type.labs <- c("EG", "FG", "FP24H", "FP2H", "SW")
names(sample.type.labs) <- c("EmptyGut", "FullGut", "FecalPellets24", "FecalPellets2", "Water")

exp.type.labs <- c("Experiment 1", "Experiment 2", "Experiment 3")
names(exp.type.labs) <- c("Experiment1", "Experiment2", "Experiment3")

#Heatmaps
theme_set(theme_bw())
heatmap_z_score_all_ec_number_top100 <- ggplot(z_score_all_ec_number_group_long_merged_top100, aes(x = Sample, y = description, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "Sample Type", y= "EC description") + # Add a nicer x-axis title
  facet_nested(. ~ Experiment+SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs, Experiment = exp.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5,size = 8, color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  #scale_x_discrete(labels = sample.type.labs) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_z_score_all_ec_number_top100


heatmap_z_score_exp1_ec_number <- ggplot(z_score_exp1_ec_number_group_long_merged, aes(x = SampleType, y = EC_number, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL", y= "EC description") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5,size = 8, color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_z_score_exp1_ec_number

heatmap_z_score_exp1_ec_number_top20 <- ggplot(z_score_exp1_ec_number_group_long_merged_top20, aes(x = SampleType, y = description, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL", y= "EC description") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5,size = 8, color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_z_score_exp1_ec_number_top20

#Heatmaps
heatmap_z_score_exp2_ec_number_top20 <- ggplot(z_score_exp2_ec_number_group_long_merged_top20, aes(x = Sample, y = description, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL", y= "EC description") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_z_score_exp2_ec_number_top20
ggsave("results/heatmap_z_score_exp2_ec_number_top20.tiff", width = 8, height = 4, dpi = 150) # save graphic

#Heatmaps
heatmap_z_score_exp3_ec_number_top20 <- ggplot(z_score_exp3_ec_number_group_long_merged_top20, aes(x = Sample, y = description, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL", y= "EC description") + # Add a nicer x-axis title
  facet_nested(~ SampleType, scales = "free", labeller = labeller(SampleType = sample.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 9, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_z_score_exp3_ec_number_top20
ggsave("results/heatmap_z_score_exp3_ec_number_top20.tiff", width = 8, height = 4, dpi = 150) # save graphic

heatmap_ec_number_top20_all_exp <- ggarrange(heatmap_z_score_exp1_ec_number_top20, heatmap_z_score_exp2_ec_number_top20, heatmap_z_score_exp3_ec_number_top20,
                                             heights = c(20, 20, 20), ncol = 1, nrow = 3, align = "v", common.legend = TRUE,
                                             legend = "right", labels = c("A", "B", "C"))
heatmap_ec_number_top20_all_exp
ggsave("results/heatmap_ec_number_top20_all_exp.tiff", width = 8, height = 10, dpi = 150) # save graphic
ggsave("results/heatmap_ec_number_top20_all_exp.pdf", width = 8, height = 10, dpi = 150) # save graphic
