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

library("readxl")
library("ggplot2")
library("ggthemes")
library("extrafont")
library("plyr")
library("scales")
library("ggpubr")
library("RColorBrewer")
library("gtable")
library("grid")
library("cowplot")
#library("dplyr")
library("vegan")
library("stringr")
library("tidyr")
library("ggord")

source("scr/summary.R")

###################### Set plotting theme to blank ######################
theme_set(theme_bw())

###################### Import all input files ######################
# Import qPCR table in excel format and convert to a data frame
q_pcr <- read_excel("data/q-PCR.xlsx")
class(q_pcr)
head(q_pcr)

core <- read_excel("data/eg_core.xlsx")
class(core)
head(core)

# Summary +/-SE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
q_pcr_sum <- summarySE(q_pcr, measurevar=c("Abundance"), groupvars=c("SampleType"))

core_sum <- summarySE(core, measurevar=c("Abundance"), groupvars=c("Experiment", "Genus"))
log_core_sum <- summarySE(core, measurevar=c("Log_Abundance"), groupvars=c("Experiment", "Genus"))
ra_core_sum <- summarySE(core, measurevar=c("RA"), groupvars=c("Experiment", "Genus"))

my_comparisons <- list(c("FP24Hrs", "EG"), c("FP24Hrs", "SW"),
                       c("EG", "SW")) 
# Barplots of environmental variables by basin
q_pcr_sum_plot <- ggplot(q_pcr_sum, aes(x= SampleType, y= Abundance, fill = SampleType)) + 
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Abundance-se, ymax=Abundance+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  stat_compare_means(comparisons = my_comparisons, p.adjust.method = "none", aes(label = ..p.signif..), size = 3, hide.ns = TRUE)+
  stat_compare_means(size = 3, label.y = 0)+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#666666")) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     expand = c(0.05, 0.05)) +
  coord_cartesian(ylim = c(10, 10000000)) +
  theme(axis.title.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axiselement_blank()) + # adjusts text of x axis
  xlab("Sample") + # add the title on y axis
  ylab(expression(bold("16S rRNA gene copies"~(mm^3)))) +
  theme(legend.position="none") +
  #theme(legend.title = element_text(face = "bold")) +
  #theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
q_pcr_sum_plot
ggsave("results/16S_q_pcr_sum_plot.pdf", width = 4, height = 3, dpi = 200) # save graphic

# KW analysis 
kw_q_pcr_abundance <- kruskal.test(Abundance ~ SampleType, data = q_pcr)
pwt_q_pcr_abundance_no_adjust <- pairwise.t.test(q_pcr$Abundance, q_pcr$SampleType, p.adjust.method = "none")
pwt_q_pcr_abundance_adjust <- pairwise.t.test(q_pcr$Abundance, q_pcr$SampleType, p.adjust.method = "BH")

#### plotting core microbiome
core_sum_plot <- factor(core_sum$Experiment,
                        levels = c("Experiment 1", "Experiment 2", "Experiment 3"))

log_core_sum_plot <- ggplot(log_core_sum, aes(x= Experiment, y= Log_Abundance, fill = Genus)) + 
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Log_Abundance-se, ymax=Log_Abundance+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  facet_wrap(Genus ~ ., scales = "free", ncol = 3) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02",
                               "#A6761D", "#666666", "lightblue")) +
  scale_x_discrete(
    labels = c("Exp1", "Exp2", "Exp3"),
    drop = TRUE) +
  #scale_y_continuous(trans = "log10",
   #                  breaks = trans_breaks("log10", function(x) 10^x),
    #                 labels = trans_format("log10", math_format(10^.x)),
     #                expand = c(0.05, 0.05)) +
  #coord_cartesian(ylim = c(0, 1000000)) +
  theme(axis.title.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axiselement_blank()) + # adjusts text of x axis
  xlab("Experiment") + # add the title on y axis
  ylab("Abundance") +
  theme(legend.position="none") +
  #theme(legend.title = element_text(face = "bold")) +
  #theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
theme(strip.background =element_rect(
  color = "black",
  fill = "white",
  size = 1,
  linetype = "solid"),
  strip.text = element_text(
    size = 12, color = "black", face = "bold.italic")) 
log_core_sum_plot
ggsave("results/16S_log_core_sum_plot.pdf", width = 7, height = 6, dpi = 200) # save graphic

ra_core_sum_plot <- ggplot(ra_core_sum, aes(x= Experiment, y= RA, fill = Genus)) + 
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=RA-se, ymax=RA+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  facet_wrap(Genus ~ ., scales = "free", ncol = 3) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02",
                               "#A6761D", "#666666", "lightblue")) +
  scale_x_discrete(
    labels = c("Exp1", "Exp2", "Exp3"),
    drop = TRUE) +
  scale_y_continuous(labels= scales::percent_format(accuracy = 0.1L)) + # plot as % and removes the internal margins
  #scale_y_continuous(trans = "log10",
  #                  breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                expand = c(0.05, 0.05)) +
  #coord_cartesian(ylim = c(0, 1000000)) +
  theme(axis.title.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axiselement_blank()) + # adjusts text of x axis
  xlab("Experiment") + # add the title on y axis
  ylab("Relative Abundance") +
  theme(legend.position="none") +
  #theme(legend.title = element_text(face = "bold")) +
  #theme(legend.title.align = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold.italic")) 
ra_core_sum_plot
ggsave("results/16S_ra_core_sum_plot.pdf", width = 7, height = 6, dpi = 200) # save graphic
