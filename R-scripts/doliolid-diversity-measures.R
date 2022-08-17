###################### Setting working directory ######################
setwd("/Users/tiagopereira/Documents/doliolid")

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

theme_set(theme_bw())

# Import matrices using excel files - diversity measures
alpha_div_all <- read_excel("exported_tables/alpha_div.xlsx")

# Summary +/-SE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
# Summary for environmental variables across basins
alpha_div_all_sum_ob <- summarySE(alpha_div_all, measurevar=c("Observed"), groupvars=c("SampleType", "Experiment"))
alpha_div_all_sum_sh <- summarySE(alpha_div_all, measurevar=c("Shannon"), groupvars=c("SampleType", "Experiment"))
alpha_div_all_sum_sim <- summarySE(alpha_div_all, measurevar=c("Simpson"), groupvars=c("SampleType", "Experiment"))

alpha_div_all_sum_ob$SampleType <- factor(alpha_div_all_sum_ob$SampleType,
                               level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"),
                               labels = c("SW", "FG", "FP2Hrs", "FP24Hrs", "EG"))

alpha_div_all_sum_sh$SampleType <- factor(alpha_div_all_sum_sh$SampleType,
                                          level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"),
                                          labels = c("SW", "FG", "FP2Hrs", "FP24Hrs", "EG"))

alpha_div_all_sum_sim$SampleType <- factor(alpha_div_all_sum_sim$SampleType,
                                          level = c("Water", "FullGut", "FecalPellets2", "FecalPellets24", "EmptyGut"),
                                          labels = c("SW", "FG", "FP2Hrs", "FP24Hrs", "EG"))

alpha_color <- c("#666666", "#D59D08","#D95F02", "#7570B3", "#1B9E77")

# Barplots of diversity indices
  alpha_div_all_sum_ob_plot <- ggplot(alpha_div_all_sum_ob, aes(x= SampleType, y= Observed, fill = SampleType)) + 
  facet_grid(. ~ factor(Experiment, levels = c("Experiment1", "Experiment2", "Experiment3"), labels = c("Experiment 1", "Experiment 2", "Experiment 3")), scales = "free") +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("#666666", "#D59D08","#D95F02", "#7570B3", "#1B9E77")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) + # adjusts text of x axis
  theme(axis.ticks.x = element_blank())+
  scale_y_continuous(expand = c(0.0, 0.5),
                     limits = c(0, 400)) +
  xlab("Sample Type") + # add the title on y axis
  ylab("ASV Richness") +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ # removes the gridlines
  theme(legend.position="none") +
    theme(strip.background =element_rect(
      color = "black",
      fill = "white",
      size = 1,
      linetype = "solid"),
      strip.text = element_text(
        size = 12, color = "black", face = "bold"))
alpha_div_all_sum_ob_plot

alpha_div_all_sum_sh_plot <- ggplot(alpha_div_all_sum_sh, aes(x= SampleType, y= Shannon, fill = SampleType)) + 
  facet_grid(. ~ factor(Experiment, levels = c("Experiment1", "Experiment2", "Experiment3"), labels = c("Experiment 1", "Experiment 2", "Experiment 3")), scales = "free") +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("#666666", "#D59D08","#D95F02", "#7570B3", "#1B9E77")) +
  theme(strip.text.x = element_blank()) +
  theme(strip.text.y = element_blank()) +
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) + # adjusts text of x axis
  theme(axis.ticks.x = element_blank())+
  scale_y_continuous(expand = c(0.0, 0.0)) +
  xlab("Sample Type") + # add the title on y axis
  ylab("Shannon") +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ # removes the gridlines
  theme(legend.position="none") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
alpha_div_all_sum_sh_plot

alpha_div_all_sum_sim_plot <- ggplot(alpha_div_all_sum_sim, aes(x= SampleType, y= Simpson, fill = SampleType)) + 
  facet_grid(. ~ factor(Experiment, levels = c("Experiment1", "Experiment2", "Experiment3"), labels = c("Experiment 1", "Experiment 2", "Experiment 3")), scales = "free") +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Simpson-se, ymax=Simpson+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("#666666", "#D59D08","#D95F02", "#7570B3", "#1B9E77")) +
  theme(strip.text.x = element_blank()) +
  theme(strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 12, color = "black", face = "bold")) + # adjusts the title of x axis
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.0, 0.0),
                     limits = c(0, 40)) +
  xlab("Sample Type") + # add the title on y axis
  ylab("Simpson") +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ # removes the gridlines
  theme(legend.position="none") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
alpha_div_all_sum_sim_plot

alpha_div_all_figure <- ggarrange(alpha_div_all_sum_ob_plot, alpha_div_all_sum_sh_plot, alpha_div_all_sum_sim_plot,
                                  labels = c("A", "B", "C"), nrow = 3, ncol = 1, legend = FALSE,
                                  align = "v")
alpha_div_all_figure

ggsave("results/alpha_div_all_figure.pdf", width = 5, height = 7, dpi = 150) # save graphic
