#!/usr/bin/env Rscript
# Author: Jiqiu Wu
# Contact: jiqiuwuwhy@gmail.com
# This is an integrated R script to plot figures of the CRISPR study!
# Enjoy  ~

# Throw everything away and restart
rm(list=ls())
graphics.off()

# load libraries
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ape)
library(ggsci)
library(vegan)
library(Rmisc)
library(stringr)


# Fig 1 -------------------------------------------------------------------

# Fig 1a
setwd("~/Documents/ImperialCollege/CRISPR/Codes/")
contig_species <- read.csv("../Files/contig_species.csv", header = F, sep = ",")
colnames(contig_species) <- c("baby_id", "month", "contig_id", "acc_num", "phylum", 
                              "genus", "species")
genus_data <- contig_species %>% group_by(genus) %>% summarise(n()) %>% as.data.frame()
colnames(genus_data) <- c("genus","value")

phylum_data <- contig_species %>% group_by(phylum) %>% summarise(n()) %>% as.data.frame()
colnames(phylum_data) <- c("phylum","value")

phylum_data = phylum_data[order(phylum_data$value, decreasing = TRUE),]
myLabel = as.vector(phylum_data$phylum)   
myLabel = paste(myLabel, "(", round(phylum_data$value / sum(phylum_data$value) * 100, 2), "%)", sep = "")   

colnames(phylum_data)[1] <- "Phylum"

phylum_plot <-  ggplot(phylum_data, aes(x = "", y = value, fill = Phylum)) +
  geom_bar(stat = "identity") +    
  coord_polar(theta = "y") +  theme_classic(6) +
  labs(x = "", y = "", title = "") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586", "#F5DF3D", "#93CA76", "#9D9D9D")) + 
   theme(panel.border = element_rect( fill = NA))

pdf("../Figures/Fig1/a_phylum_plot.pdf", width = 7/2.54, height = 6/2.54)
phylum_plot
graphics.off()

# Fig 1b
genus_data = genus_data[order(genus_data$value, decreasing = TRUE),]

genus2plot = subset(genus_data, value > 5)

genus_plot <-  ggplot(genus2plot, mapping = aes(x = reorder(genus, -value), y = value)) +
  geom_bar(stat = "identity") + 
  theme_classic(6) +
  theme(axis.text.x = element_text(angle=60, vjust = 0.5))   + 
  theme(panel.border = element_rect( fill = NA)) +
  labs(x = "Genus", y = "Number of CRISPRs", title = "") 


pdf("../Figures/Fig1/b_genus_summary.pdf", width = 10.8/2.54, height = 6/2.54)
genus_plot
graphics.off()

# Fig 1c
data_bacteria <- read.table("../Files/Baby_merged.txt", header=T, sep="\t", 
                            comment.char="", row.names = 1)
data_bacteria$levels <- str_count(rownames(data_bacteria), "\\|")
data_bacteria_species <- subset(data_bacteria, levels == 6)
data_bacteria_species$levels <- NULL

# contain CRISPRs
data_bacteria_species$speices <- str_split_fixed(rownames(data_bacteria_species), 
                                                 "\\|", n = 7)[,7]
data_bacteria_species$speices <- substr(data_bacteria_species$speices, 4, 70)
contig_species$bacteria_name <- paste(contig_species$genus, contig_species$species, sep = "_")

data_bacteria_cris <- subset(data_bacteria_species, data_bacteria_species$speices %in% contig_species$bacteria_name)
data_bacteria_cris$speices <- NULL

bacteria_cris_abundance <- colSums(data_bacteria_cris)
data_bacteria_cris <- rbind(data_bacteria_cris, bacteria_cris_abundance)

data_bacteria_cris_abundance <- subset(data_bacteria_cris, row.names(data_bacteria_cris) == 64)
rownames(data_bacteria_cris_abundance)[1] <- "Total_abundance"

data_abundance_plot <- data.frame(t(data_bacteria_cris_abundance))
data_abundance_plot$month <- str_split_fixed(rownames(data_abundance_plot), 
                                             "\\_", n = 2)[,2]

mean <- data_abundance_plot %>% group_by(month) %>% summarise(mean(Total_abundance))
sd <- data_abundance_plot %>% group_by(month) %>% summarise(sd(Total_abundance))


# plot 
plot_bacteria_cris_abundance <- ggplot(data = data_abundance_plot, 
                                       aes(x = factor(month, levels = c("B", "4M", "12M")), 
                                           y = Total_abundance)) + 
  geom_violin(aes(fill = month)) + geom_boxplot(width = 0.1) + 
  ylab("Total abundance") + xlab("Age (month)") + 
  theme_classic(base_size = 6) + #theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig1/c_bacteria_cris_abundance.pdf", width = 6.2/2.54, height = 5.5/2.54)
plot_bacteria_cris_abundance
graphics.off()

# Statistical test
model_1c <- aov(Total_abundance ~ month, data = data_abundance_plot)
summary(model_1c)

TukeyHSD(model_1c, conf.level=.95)


# Fig 2 -------------------------------------------------------------------

# Fig 2a
fig2_data <- read.csv("../Files/Crispr_meta_quanlity.csv", header = T, sep = ",")
CRIS_data <- subset(fig2_data, month != "M")



plot_crispr <- ggplot(data = CRIS_data, aes(x = factor(month, levels = c("B", "4M", "12M")), 
                                            y = num_CRISPRs)) + 
  geom_violin(aes(fill = month)) + geom_boxplot(width = 0.1) + 
  ylab("Number of CRISPRs") + xlab("Age (month)") + ylim(0, 40) + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig2/a_crispr_num.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_crispr
graphics.off()

# Statistical test
model_2a <- aov(num_CRISPRs ~ month, data = CRIS_data)
summary(model_2a)

TukeyHSD(model_2a, conf.level=.95)

# Fig 2b
richness_data <- read.csv("../Files/bacterial_richness.csv", header = T, sep = ",")
richness2com <- subset(richness_data, Age != "M")
cris_data2com <- select(CRIS_data, baby_id, gender)
colnames(cris_data2com)[1] <- "BabyID"
cris_data2com <- unique(cris_data2com)
richness2plot <- left_join(richness2com, cris_data2com, by = "BabyID")


cris_model <- select(CRIS_data, baby_id, month, num_CRISPRs, spacers_sample)
colnames(cris_model)[1:2] <- c("BabyID", "Age")
data_model <- left_join(cris_model, richness2plot, by = c("BabyID", "Age"))

plot_model_cris <- ggplot(data = data_model, aes(x = Bacterial_richness, y = num_CRISPRs, colour = Age)) + 
  geom_point(size = 0.5) + geom_smooth(aes(group = 1) , method = 'lm')  + 
  ylab("The number of CRISPRs") + xlab("Bacterial richness") + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA)) + 
  scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none")

pdf("../Figures/Fig2/b_model_cris.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_model_cris 
graphics.off()

fit_2b <- lm(num_CRISPRs ~ Bacterial_richness, data = data_model)
summary(fit_2b)

# Fig 2c
plot_crispr_density <- ggplot(data = CRIS_data, aes(x = factor(month, levels = c("B", "4M", "12M")), 
                                                    y = num_CRISPRs/total_bp * 1000000)) + 
  geom_violin(aes(fill = month)) + geom_boxplot(width = 0.1) + 
  ylab("Density of CRISPRs") + xlab("Age (month)") + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none")

pdf("../Figures/Fig2/c_crispr_density.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_crispr_density
graphics.off()

# Statistical test
model_2c <- aov(num_CRISPRs/total_bp ~ month, data = CRIS_data)
summary(model_2c)

TukeyHSD(model_2c, conf.level=.95)

# Fig 2d
plot_spacer <- ggplot(data = CRIS_data, aes(x = factor(month, levels = c("B", "4M", "12M")), y = spacers_sample)) + 
  geom_violin(aes(fill = month)) + geom_boxplot(width = 0.1) + 
  ylab("Number of spacers") + xlab("Age (month)") + ylim(0, 500) +
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig2/d_spacer_num.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_spacer
graphics.off()

# Statistical test
model_2d <- aov(spacers_sample ~ month, data = CRIS_data)
summary(model_2d)

TukeyHSD(model_2d, conf.level=.95)

# Fig 2e
plot_model_spacer <- ggplot(data = data_model, aes(x = Bacterial_richness, y = spacers_sample, colour = Age)) + 
  geom_point(size = 0.5) + geom_smooth(aes(group = 1) , method = 'lm')  + 
  ylab("The number of spacers") + xlab("Bacterial richness") + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA)) + 
  scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "none")

pdf("../Figures/Fig2/e_model_spacer.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_model_spacer 
graphics.off()

fit_2e <- lm(spacers_sample ~ Bacterial_richness, data = data_model)
summary(fit_2e)

# Fig 2f
plot_spacer_density <- ggplot(data = CRIS_data, aes(x = factor(month, levels = c("B", "4M", "12M")), 
                                                    y = spacers_sample/num_CRISPRs)) + 
  geom_violin(aes(fill = month)) + geom_boxplot(width = 0.1) + 
  ylab("Number of spacers in one CRISPR") + xlab("Age (month)") + 
  theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#9EBC19")) +
  theme(legend.position = "bottom")

pdf("../Figures/Fig2/f_spacer_density.pdf", width = 4.3/2.54, height = 5.6/2.54)
plot_spacer_density
graphics.off()


# Fig 3 -------------------------------------------------------------------

# Fig 3a
# input 
sankey_data = read.table("../Files/sankey_input.txt", 
                         header=T, sep="\t", comment.char="", row.names = 1) 
colnames(sankey_data) <- str_sub(colnames(sankey_data), 4,20)

sankey_plot <- data.frame(t(sankey_data))
sankey_plot$genus <- rownames(sankey_plot)
sankey_plot <- melt(sankey_plot, id.vars =  "genus")
sankey_plot$Value <- sqrt(sankey_plot$value)
sankey_plot$value <- NULL
sankey_plot$variable <- as.numeric(substr(sankey_plot$variable, 2, 5))

set.seed(39) # for nice colours
cols <- hsv(h = sample(1:20/20), s = sample(3:12)/15, v = sample(3:12)/15)


pdf("../Figures/Fig3/a_sankey_genus_20.pdf", width = 26/2.54, height = 15/2.54)
plot_sankey <- alluvial_ts(sankey_plot, wave = .3, ygap = 5, col = cols, plotdir = 'centred', alpha=.9,
                           grid = TRUE, grid.lwd = 5, xmargin = 0.2, lab.cex = .7, xlab = '',
                           ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
                           leg.col='white', 
                           title = "")
graphics.off()

# Fig 3b
# Bacteria containing CRISPRs 
genus_gender <- read.csv("../Files/genus_gender.csv", sep = ",", header = T)
genus_gender_num <- genus_gender %>% group_by(genus, month) %>% summarise(frequency = sum(value)) %>% as.data.frame()
genus_gender_num$value <- sqrt(genus_gender_num$frequency)


genus_gender_num$month <- gsub("4M", as.numeric(4), genus_gender_num$month)
genus_gender_num$month <- gsub("12M", as.numeric(12), genus_gender_num$month)
genus_gender_num$month <- gsub("B", as.numeric(0), genus_gender_num$month)

genus_gender_num$month <- as.numeric(genus_gender_num$month)

cols_1 <- hsv(h = sample(1:48/48), s = sample(3:12)/15, v = sample(3:12)/15)
pdf("../Figures/Fig3/b_sankey_genus_crispr.pdf", width = 26/2.54, height = 25/2.54)
plot_sankey_crispr_bacteria <- alluvial_ts(genus_gender_num, wave = .3, ygap = 5, col = cols_1, plotdir = 'centred', alpha=.9,
                                           grid = TRUE, grid.lwd = 5, xmargin = 0.2, lab.cex = .7, xlab = '',
                                           ylab = '', border = NA, axis.cex = .8, leg.cex = .7,
                                           leg.col='white', 
                                           title = "")
graphics.off()

# Fig 3c

tax_sample <- read.table("../Files/Baby_genus.txt", header=T, sep="\t", 
                        comment.char="", row.names = 1) 
micro_data <- tax_sample
micro_data$bacteria <- rownames(micro_data)
micro_data$genus <- str_split_fixed(micro_data$bacteria, "\\.", n = 6)[,6]
micro_data$genus_pure <- substr(micro_data$genus, 4, 20)
micro_data_cris <- subset(micro_data, genus_pure %in% unique(genus_gender_num$genus))
micro_data_cris$bacteria <- NULL
micro_data_cris$genus <- NULL
rownames(micro_data_cris) <- micro_data_cris$genus_pure
micro_data_cris$genus_pure <- NULL
micro_data_cris$Taxonomy <- NULL
micro_data_cris_t <- data.frame(t(micro_data_cris))


micro_genus_12 <- subset(micro_data_cris_t, 
                         str_split_fixed(colnames(micro_data_cris), "\\_", n = 2)[,2] == "12M")
micro_genus_4 <- subset(micro_data_cris_t, 
                        str_split_fixed(colnames(micro_data_cris), "\\_", n = 2)[,2] == "4M")
micro_genus_b <- subset(micro_data_cris_t, 
                        str_split_fixed(colnames(micro_data_cris), "\\_", n = 2)[,2] == "B")
micro_genus_cris <- data.frame(mean_12 = colMeans(micro_genus_12), 
                               mean_4 = colMeans(micro_genus_4),
                               mean_0 = colMeans(micro_genus_b))
micro_genus_cris$genus <- rownames(micro_genus_cris)

micro_genus_plot <- melt(micro_genus_cris, id.vars = "genus")
micro_genus_plot$variable <- str_split_fixed(micro_genus_plot$variable, "\\_", n = 2)[,2]

# association


colnames(genus_gender_num)[3] <- "num_crispr"
colnames(micro_genus_plot)[3] <- "bact_relative"


genus_gender_num$month <- gsub("4M", as.numeric(4), genus_gender_num$month)
genus_gender_num$month <- gsub("12M", as.numeric(12), genus_gender_num$month)
genus_gender_num$month <- gsub("B", as.numeric(0), genus_gender_num$month)


genus_gender_num$by_id <- paste(genus_gender_num$genus, genus_gender_num$month, sep = "_")
micro_genus_plot$by_id <- paste(micro_genus_plot$genus, micro_genus_plot$variable, sep = "_")



associa_plot <- full_join(micro_genus_plot, genus_gender_num, by = "by_id")
associa_plot <- drop_na(associa_plot)


associa_plot <- subset(associa_plot, by_id != "Lactococcus_0")

associa_plot$month <- gsub(4, "4M", associa_plot$month)
associa_plot$month <- gsub(12, "12M", associa_plot$month)
associa_plot$month <- gsub(0, "B", associa_plot$month)

plot_ass <- ggplot(data = associa_plot, aes(x = bact_relative, y = num_crispr, color = as.factor(month))) +
  geom_point(size = 0.5)  + geom_smooth(aes(group = 1), method = "lm") + theme_classic(base_size = 6) +  
  theme(panel.border = element_rect(fill = NA)) + ylab("Number of CRISPRs") + 
  xlab("Bacterial relative abundance") + 
  scale_color_manual(values = c("#C35C6A", "#88ABDA","#9EBC19")) + 
  theme(legend.position = "bottom")


pdf("../Figures/Fig3/c_associa.pdf", width = 4.3/2.54, height = 4.6/2.54)
plot_ass
graphics.off()

fit_3c <- lm(num_crispr ~ log(bact_relative), data = associa_plot )
summary(fit_3c)

fit_3c_lm <- lm(num_crispr ~ bact_relative, data = associa_plot )
summary(fit_3c_lm)

# Fig 3e

indiv_shared <- read.csv("../Files/individual_shared.csv", header = F, sep = ",")
colnames(indiv_shared) <- c("baby_id", "B", "M4", "M12","B_4M", "B_12M", "b4M_12M",
                            "B_4M_12M")

indi_shared <- melt(indiv_shared, id = "baby_id")
colnames(indi_shared)[2:3] <- c("Month", "spacers")


indi_shared$Month <- gsub("B_4M_12M", "B & 4M & 12M",  indi_shared$Month)
indi_shared$Month <- gsub("b4M_12M", "4M & 12M",  indi_shared$Month)
indi_shared$Month <- gsub("B_12M", "B & 12M",  indi_shared$Month)
indi_shared$Month <- gsub("B_4M", "B & 4M",  indi_shared$Month)
indi_shared$Month <- gsub("M12", "12M",  indi_shared$Month)
indi_shared$Month <- gsub("M4", "4M", indi_shared$Month)
indi_shared$Month <- gsub("B", "B", indi_shared$Month)

indi_shared2plot <- subset(indi_shared, spacers != 0)

plot_share <- ggplot(data = indi_shared2plot) + 
  geom_violin(aes(x = Month, y = sqrt(spacers), fill = Month))   + 
  geom_boxplot(aes(x = Month, y = sqrt(spacers)), width = 0.05) + 
  ylab("The square root of numbers of spacers") + xlab("") + 
  scale_x_discrete(limits = c("B", "4M", "12M","B & 4M","B & 12M", "4M & 12M", "B & 4M & 12M")) +
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA), 
                                       axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + 
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#422256", "#9EBC19",
                               "#422256","#422256","#422256")) +
  theme(legend.position = "none") 
  
  
  ggplot(data = indi_shared2plot) +
  geom_violin(aes(x = Month, y = spacers, fill = Month)) +
  geom_boxplot(aes(x = Month, y = spacers), outlier.shape = NA, width = 0.05) + ylab("Intersection number") + 
  scale_x_discrete(limits = c("B", "4M", "12M",
                              "B & 4M","B & 12M",
                              "4M & 12M", "B & 4M & 12M")) + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect(fill = NA)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))  +  
  scale_fill_manual(values = c("#C35C6A", "#88ABDA", "#422256", "#9EBC19",
                               "#422256","#422256","#422256"))  +
  xlab("") + ylab("The square root of numbers of spacers") + theme(legend.position = "none")

pdf("../Figures/Fig3/e_indi_shared.pdf", width = 6.3/2.54, height = 4.6/2.54)
plot_share
graphics.off()

# Statistical test
model_3e <- aov(spacers ~ Month, data = indi_shared2plot)
summary(model_3e)

TukeyHSD(model_3e, conf.level=.95)

indiv_shared <- read.csv("../Files/individual_shared.csv", header = F, sep = ",")
colnames(indiv_shared) <- c("baby_id", "B", "M4", "M12","B_4M", "B_12M", "b4M_12M",
                            "B_4M_12M")

indi_shared <- melt(indiv_shared, id = "baby_id")
colnames(indi_shared)[2:3] <- c("Month", "spacers")


indi_shared$Month <- gsub("B_4M_12M", "At birth & 4 months & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("b4M_12M", "4 months & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("B_12M", "At birth & 12 months",  indi_shared$Month)
indi_shared$Month <- gsub("B_4M", "At birth & 4 months",  indi_shared$Month)
indi_shared$Month <- gsub("M12", "12 months",  indi_shared$Month)
indi_shared$Month <- gsub("M4", "4 months", indi_shared$Month)
indi_shared$Month <- gsub("B", "At birth", indi_shared$Month)

indi_shared2plot <- subset(indi_shared, spacers != 0)

plot_share <- ggplot(data = indi_shared2plot, aes(x = Month, y = spacers, colour = Month)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1) + ylab("Intersection number") + 
  scale_x_discrete(limits = c("At birth", "4 months", "12 months",
                              "At birth & 4 months","At birth & 12 months",
                              "4 months & 12 months", "At birth & 4 months & 12 months")) + theme_classic(base_size = 8) +
  theme(axis.text.x =element_text(angle = 60, hjust =1, vjust = 1))  +  
  scale_colour_manual(values = c("#C35C6A", "#88ABDA", "#422256", "#9EBC19","#422256","#422256","#422256"))  +
  xlab("") + ylab("The numbers of spacers")

pdf("../Figures/Fig3/b_indi_shared.pdf", width = 9.2/2.54, height = 7/2.54)
plot_share
graphics.off()


# Fig 3f
# jaccard index for bacteria
# species level
metaphlan_all <- read.table("../Files/Baby_merged.txt", sep = "\t", header = T, row.names = 1)
metaphlan_all$levels <- str_count(rownames(metaphlan_all), "\\|")

metaphlan_species <- subset(metaphlan_all, levels == 6)
metaphlan_species$levels <- NULL

species_num <- nrow(metaphlan_species)
sample_num <- ncol(metaphlan_species)

species_birth <- NULL
species_4M <- NULL
species_12M <- NULL


for (j in 1:sample_num){
  {for (i in 1:species_num)
    if (metaphlan_species[i,j] != 0)
      metaphlan_species[i,j] = rownames(metaphlan_species)[i]
  }
}

for (k in 1:sample_num){
  if (str_split_fixed(colnames(metaphlan_species)[k], "\\_", n = 2)[,2] == "B"){
    species_birth <- union(species_birth, unique(metaphlan_species[,k]))
  } else if (str_split_fixed(colnames(metaphlan_species)[k], "\\_", n = 2)[,2] == "4M"){
    species_4M <- union(species_4M, unique(metaphlan_species[,k]))
  } else {species_12M <- union(species_12M, unique(metaphlan_species[,k]))}
}

num2remove <- c("0")
species_birth = species_birth[!(species_birth %in% num2remove)]
species_4M = species_4M[!(species_4M %in% num2remove)]
species_12M = species_12M[!(species_12M %in% num2remove)]


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


# jaccard index for spacer
ja_b_4 <- (748+555) / (748+555+1399+193+3852+1074)
ja_4_12 <- (555+1074) / (555+1074+748+3852+193+9448)

# for individuals
individual_venn <- data.frame(matrix(ncol = 2, nrow = 82))
colnames(individual_venn) <- c("jac_b_4", "jac_4_12")
rownames(individual_venn) <- unique(str_split_fixed(colnames(metaphlan_species), 
                                                    "\\_", n = 2)[,1])

for (k in rownames(individual_venn)){
  for (l in colnames(metaphlan_species)){
    if (str_split_fixed(l, "\\_", n = 2)[,1] == k & 
        str_split_fixed(l, "\\_", n = 2)[,2] == "B"){
      species_birth <- unique(metaphlan_species[,l])
      species_b <- species_birth[!(species_birth %in% num2remove)]
    } else if (str_split_fixed(l, "\\_", n = 2)[,1] == k & 
               str_split_fixed(l, "\\_", n = 2)[,2] == "4M"){
      species_4m <- unique(metaphlan_species[,l])
      species_4 <- species_4m[!(species_4m %in% num2remove)]
    } else if (str_split_fixed(l, "\\_", n = 2)[,1] == k & 
               str_split_fixed(l, "\\_", n = 2)[,2] == "12M"){
      species_12m <- unique(metaphlan_species[,l])
      species_12 <- species_12m[!(species_12m %in% num2remove)]
    }
  }
  individual_venn[k, "jac_b_4"] <- jaccard(species_b, species_4)
  individual_venn[k, "jac_4_12"] <- jaccard(species_4, species_12)
}

individual_venn$baby_id <- substr(rownames(individual_venn),2, 100L)
individual_venn$Type <- "Species"

# spacer
indiv_shared <- read.csv("../Files/individual_shared.csv", 
                         header = F, sep = ",")
colnames(indiv_shared) <- c("baby_id", "B", "M4", "M12","B_4M", "B_12M", "b4M_12M", "B_4M_12M")

indiv_shared$jac_b_4 <- rowSums(indiv_shared[,c('B_4M', 'B_4M_12M')]) / 
  rowSums(indiv_shared[,c('B_4M', 'B_4M_12M', 'B', 'M4', 'B_12M', 'b4M_12M')])

indiv_shared$jac_4_12 <- rowSums(indiv_shared[,c('b4M_12M', 'B_4M_12M')]) / 
  rowSums(indiv_shared[,c('b4M_12M', 'B_4M_12M', 'M4', 'b4M_12M', "M12", 'B_12M')])

indiv_shared_plot <- select(indiv_shared, baby_id, jac_b_4, jac_4_12)
indiv_shared_plot$Type <- "Spacer"

indi_overlap <- rbind(indiv_shared_plot, individual_venn)

indi_overlap_plot <- melt(indi_overlap, id = c("baby_id", "Type"))
colnames(indi_overlap_plot)[3:4] <- c("Time", "Jaccard")

plot_indi_overlap <- ggplot(transform(indi_overlap_plot, 
                                      Time = factor(Time, levels = c("jac_b_4", "jac_4_12")))) + 
  geom_violin(aes(x = Type, y = Jaccard, fill = Type))   + facet_wrap( ~ Time) +
  geom_boxplot(aes(x = Type, y = Jaccard), width = 0.1) + ylim(0, 0.6) +
  ylab("Jaccard Index") + xlab("") + 
  theme_classic(base_size = 6) + theme(panel.border = element_rect( fill = NA),
                                       axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + 
  scale_fill_manual(values = c("#89c3eb", "#c97586")) +
  theme(legend.position = "none") 

pdf("../Figures/Fig3/f_overlap_spacer_species.pdf", width = 3.3/2.54, height = 4.6/2.54)
plot_indi_overlap
graphics.off()


# Statistical test
indi_overlap_plot_1 <- subset(indi_overlap_plot, Time == "jac_b_4")
indi_overlap_plot_2 <- subset(indi_overlap_plot, Time == "jac_4_12")


model_3f_1 <- aov(Jaccard ~ Type, data = indi_overlap_plot_1)
summary(model_3f_1)

TukeyHSD(model_3f_1, conf.level=.95)

model_3f_2 <- aov(Jaccard ~ Type, data = indi_overlap_plot_2)
summary(model_3f_2)

TukeyHSD(model_3f_2, conf.level=.95)


# Fig5 --------------------------------------------------------------------

# Fig 5a

spacer_cris <- read.csv("../Files/spacer_unique_seq.csv", 
                        sep = ",", header = F )
colnames(spacer_cris) <- c("spacerID", "spaceruqiID", "seq")
spacer_cris2com <- select(spacer_cris, "spacerID", "spaceruqiID")

# process
interaction <- spacer_cris2com
interaction$cris_spacer <- substr(interaction$spacerID, 1, 40)

interaction <- left_join(interaction, contig_species, by = "cris_spacer")
interaction <- left_join(interaction, spacer_viruse, by = "spaceruqiID")
colnames(interaction)[7] <- "bacteria_accession"
colnames(interaction)[11] <- "phage_accession"

interaction <- left_join(interaction, phage_database, by = "phage_accession")

# full interaction
interaction_full <- interaction %>% drop_na()

# generate the interaction strength
interaction_full$interac <- paste(interaction_full$genus, interaction_full$species, 
                                  interaction_full$phage_name, sep = "_")
interaction_strength <- interaction_full %>% group_by(interac) %>% summarise(n())
interaction_full_dedup <- left_join(interaction_full_dedup, interaction_strength, by = "interac")
colnames(interaction_full_dedup)[7] <- "strength"
interaction_full_dedup$node1 <- paste(interaction_full_dedup$genus, 
                                      interaction_full_dedup$species, sep = " ")

data_links_cris <- select(interaction_full_dedup, node1, phage_name, type, strength)
colnames(data_links_cris)[2] <- "node2"
data_links_cris$strength_2 <- sqrt(data_links_cris$strength)

# generate a nodes table
data_phage_com <- select(interaction_full_dedup, phage_name)
data_phage_com <- unique(data_phage_com)
data_phage_com$feature <- "Phage"
colnames(data_phage_com)[1] <- "Nodes"

data_bacte_com <- select(interaction_full_dedup, node1)
data_bacte_com <- unique(data_bacte_com)
data_bacte_com$feature <- "Bacteria"
colnames(data_bacte_com)[1] <- "Nodes"

data_nodes_cris <- rbind(data_phage_com, data_bacte_com)
data_nodes_cris$name <- data_nodes_cris$Nodes

# network
network <- graph_from_data_frame(d = data_links_cris, vertices = data_nodes_cris, directed = F)


set.seed(42)
l <- layout_with_kk(network)

# node color（V）
V(network)[feature == "Phage"]$color <- "#BFC096"
V(network)[feature == "Bacteria"]$color <- "#f39b7fff"

# node size（E）
V(network)$size <- 10
V(network)[feature == "Phage"]$size <- 5

# node label（E）
V(network)$label <- ""
V(network)[feature == "Bacteria"]$label <- V(network)[feature == "Bacteria"]$name

# node shape（V）
V(network)$shape <- "circle"


# edge width（E）
E(network)$width <- E(network)$strength_2

# edge color（E）
E(network)[type == 1]$color <- "#DC143C"
E(network)[type == 2]$color <- "#9932CC"
E(network)[type == 3]$color <- "#708090"
E(network)[type == 4]$color <- "#00CED1"


pdf("~/Documents/ImperialCollege/Revision/bacteria/interaction_network.pdf", width = 6/2.54, height = 6/2.54)
plot(network, layout = l)
graphics.off()

Betweenness <- betweenness(network)
Eig <- evcent(network)$vector

# centrality
Degree <- data.frame(degree(network))
Degree$Nodes <- rownames(Degree)

Degree <- left_join(Degree, data_nodes_cris, by = "Nodes")

Degree_bact <- subset(Degree, feature == "Bacteria")
Degree_phag <- subset(Degree, feature == "Phage")

Degree_bact <- Degree_bact[order(Degree_bact$`degree.network.`, decreasing = T),]
Degree_bact <- head(Degree_bact, 5)

# Fig 5c
plot_cen_degree <- ggplot(data = Degree_bact, 
                          aes(x = factor(Nodes, levels = c("Escherichia coli",
                                                           "Streptococcus thermophilus",
                                                           "Lactobacillus paracasei",
                                                           "Lactobacillus rhamnosus",
                                                           "Lactobacillus paracase")), 
                              y = `degree.network.`)) + geom_bar(stat = "identity", width = 0.5) +
  theme_classic(6) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))   + 
  theme(panel.border = element_rect( fill = NA)) +
  labs(x = "", y = "Degree", title = "") 

pdf("../Figures/Fig5/c_bacteria_high_degree.pdf", width = 4.3/2.54, height = 4.3/2.54)
plot_cen_degree
graphics.off()

# Fig 5d
# phage
Degree_phag <- Degree_phag[order(Degree_phag$`degree.network.`, decreasing = T),]
Degree_phag <- head(Degree_phag, 5)

plot_cen_degree_phag <- ggplot(data = Degree_phag, 
                               aes(x = factor(Nodes, levels = c("Moineauvirus",
                                                                "Lactobacillus phage phiAT3",
                                                                "Lactobacillus phage Lc-Nu",
                                                                "Lactobacillus phage Lrm1",
                                                                "Streptococcus virus 9872")), 
                                   y = `degree.network.`)) + geom_bar(stat = "identity", width = 0.5) +
  theme_classic(6) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))   + 
  theme(panel.border = element_rect( fill = NA)) +
  labs(x = "", y = "Degree", title = "") 

pdf("../Figures/Fig5/d_phage_high_degree.pdf", width = 4.3/2.54, height = 4.3/2.54)
plot_cen_degree_phag
graphics.off()
