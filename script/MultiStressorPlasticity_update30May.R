# surface design predation & copper / SHF-DFG project
# update DB 24 May 2022

# packages
library(data.table)
library(tidyverse)
library(gridExtra)
library(lme4)
library(car)
library(patchwork)
library(viridis)
library(cowplot)
library(sjPlot)
library(Hmisc)
library(rms)
library(survival)
library(survminer)
library(coxme)


###################################
#### CLONE ID AS RANDOM EFFECT ####
###################################

## DATA IMPORT ---------------------------------------------------
AllData <- read.csv("data/response_surface_final_April2021.csv")
AllData <- as.data.table(AllData)
head(AllData)
names(AllData)
AllData[, cloneID := paste0(clone_type,"_",clone)]

# induction data
Ind <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
                 ind_1_total = induction_1, 
                 ind_2_total = induction_2, 
                 ind_3_total = induction_3,
                 ind_4_total = induction_4,
                 
                 ind_1_nteeth = induction_1_nteeth, 
                 ind_2_nteeth = induction_2_nteeth, 
                 ind_3_nteeth = induction_3_nteeth,
                 ind_4_nteeth = induction_4_nteeth,
                 
                 ind_1_pedestal = induction_1_pedestal, 
                 ind_2_pedestal = induction_2_pedestal, 
                 ind_3_pedestal = induction_3_pedestal,
                 ind_4_pedestal = induction_4_pedestal,
                 
                 # max induction across all days, including instar 1 
                 maxInd_total = max_induction_total_I1,
                 maxInd_nteeth = max_induction_nteeth_I1,
                 maxInd_pedestal = max_induction_pedestal_I1,
                
                 # max induction instars 2 and later 
                 maxInd_total_I2 = max_induction_total_I2,
                 maxInd_nteeth_I2 = max_induction_nteeth_I2,
                 maxInd_pedestal_I2 = max_induction_pedestal_I2)

IndDat <- Ind[!maxInd_total=="NA"]

# induction across all instars
IndDat_use <- IndDat %>% 
  select(clone, clone_type, cloneID, juju, copper, maxInd_total) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

IndDat.ag <- na.omit(IndDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanInd = mean(maxInd_total), 
    seInd = sd(maxInd_total)/sqrt(sum(!is.na(maxInd_total))))

# induction in instars 2 and older
IndDat2_use <- IndDat %>% 
  select(clone, clone_type, cloneID, juju, copper, maxInd_total_I2) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

IndDat2.ag <- na.omit(IndDat2_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanInd = mean(maxInd_total_I2), 
    seInd = sd(maxInd_total_I2)/sqrt(sum(!is.na(maxInd_total_I2))))


# size data
Size <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
                  size_mat = size_maturity)

SizeDat <- Size[!size_mat=="NA"]

SizeDat_use <- SizeDat %>% 
  select(clone, clone_type, cloneID, juju, copper, size_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

SizeDat.ag <- na.omit(SizeDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanSize = mean(size_mat), 
    seSize = sd(size_mat)/sqrt(sum(!is.na(size_mat))))


# growth data
Growth <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
                    growth_mat = somaticGrowth_maturity, 
                    growth_3d = somaticGrowth_3d, 
                    growth_4d = somaticGrowth_4d)

GrowthDat <- Growth[!growth_mat=="NA"]

GrowthDat_use <- GrowthDat %>% 
  select(clone, clone_type, cloneID, juju, copper, growth_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

GrowthDat.ag <- na.omit(GrowthDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanGrowth = mean(growth_mat), 
    seGrowth = sd(growth_mat)/sqrt(sum(!is.na(growth_mat))))


# age data
Age <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
                 age_mat = age_maturity)

AgeDat <- Age[!age_mat=="NA"]

AgeDat_use <- AgeDat %>% 
  select(clone, clone_type, cloneID, juju, copper, age_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

AgeDat.ag <- na.omit(AgeDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanAge = mean(age_mat), 
    seAge = sd(age_mat)/sqrt(sum(!is.na(age_mat))))


## INDUCTION ---------------------------------------------------

# total induction (i.e., pedestal and neckteeth)

indSum_1 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_total),
    seInd_1 = sd(ind_1_total)/sqrt(sum(!is.na(ind_1_total))))

indSum_2 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_total),
    seInd_2 = sd(ind_2_total)/sqrt(sum(!is.na(ind_2_total))))

indSum_3 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_total),
    seInd_3 = sd(ind_3_total)/sqrt(sum(!is.na(ind_3_total))))

indSum_4 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_total),
    seInd_4 = sd(ind_4_total)/sqrt(sum(!is.na(ind_4_total))))


plot_ind1 <- ggplot(transform(indSum_1, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_1, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
                geom_point(size = 2)+
                geom_line(size = 1)+
                scale_color_viridis(discrete = TRUE, option = "D")+
                geom_line()+
                geom_errorbar(width = 0, size = 0.5)+
                facet_wrap(~clone)+
                scale_x_continuous(breaks=c(0,0.25,0.5)) +
                labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                ylim(0,110) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_blank(),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm")) 

plot_ind2 <- ggplot(transform(indSum_2, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_2, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
                geom_point(size = 2)+
                geom_line(size = 1)+
                scale_color_viridis(discrete = TRUE, option = "D")+
                geom_line()+
                geom_errorbar(width = 0, size = 0.5)+
                facet_wrap(~clone)+
                scale_x_continuous(breaks=c(0,0.25,0.5)) +
                labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                ylim(0,110) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_blank(),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm")) 
              
plot_ind3 <- ggplot(transform(indSum_3, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_3, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
                geom_point(size = 2)+
                geom_line(size = 1)+
                scale_color_viridis(discrete = TRUE, option = "D")+
                geom_line()+
                geom_errorbar(width = 0, size = 0.5)+
                facet_wrap(~clone)+
                scale_x_continuous(breaks=c(0,0.25,0.5)) +
                labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                ylim(0,110) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_blank(),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm")) 

plot_ind4 <- ggplot(transform(indSum_4, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_4, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
                geom_point(size = 2)+
                geom_line(size = 1)+
                scale_color_viridis(discrete = TRUE, option = "D")+
                geom_line()+
                geom_errorbar(width = 0, size = 0.5)+
                facet_wrap(~clone)+
                scale_x_continuous(breaks=c(0,0.25,0.5)) +
                labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                ylim(0,110) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_blank(),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm")) 


# pedestal only

indSum_1_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_pedestal),
    seInd_1 = sd(ind_1_pedestal)/sqrt(sum(!is.na(ind_1_pedestal))))

indSum_2_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_pedestal),
    seInd_2 = sd(ind_2_pedestal)/sqrt(sum(!is.na(ind_2_pedestal))))

indSum_3_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_pedestal),
    seInd_3 = sd(ind_3_pedestal)/sqrt(sum(!is.na(ind_3_pedestal))))

indSum_4_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_pedestal),
    seInd_4 = sd(ind_4_pedestal)/sqrt(sum(!is.na(ind_4_pedestal))))


plot_ind1_ped <- ggplot(transform(indSum_1_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_1, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
                    geom_point(size = 2)+
                    geom_line(size = 1)+
                    scale_color_viridis(discrete = TRUE, option = "D")+
                    geom_line()+
                    geom_errorbar(width = 0, size = 0.5)+
                    facet_wrap(~clone)+
                    scale_x_continuous(breaks=c(0,0.25,0.5)) +
                    labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                    ylim(0,60) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_blank(),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm")) 

plot_ind2_ped <- ggplot(transform(indSum_2_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_2, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
                    geom_point(size = 2)+
                    geom_line(size = 1)+
                    scale_color_viridis(discrete = TRUE, option = "D")+
                    geom_line()+
                    geom_errorbar(width = 0, size = 0.5)+
                    facet_wrap(~clone)+
                    scale_x_continuous(breaks=c(0,0.25,0.5)) +
                    labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                    ylim(0,60) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_blank(),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm")) 

plot_ind3_ped <- ggplot(transform(indSum_3_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_3, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
                    geom_point(size = 2)+
                    geom_line(size = 1)+
                    scale_color_viridis(discrete = TRUE, option = "D")+
                    geom_line()+
                    geom_errorbar(width = 0, size = 0.5)+
                    facet_wrap(~clone)+
                    scale_x_continuous(breaks=c(0,0.25,0.5)) +
                    labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                    ylim(0,60) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_blank(),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm")) 

plot_ind4_ped <- ggplot(transform(indSum_4_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_4, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
                    geom_point(size = 2)+
                    geom_line(size = 1)+
                    scale_color_viridis(discrete = TRUE, option = "D")+
                    geom_line()+
                    geom_errorbar(width = 0, size = 0.5)+
                    facet_wrap(~clone)+
                    scale_x_continuous(breaks=c(0,0.25,0.5)) +
                    labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) + 
                    ylim(0,60) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_blank(),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm")) 
                  

# neckteeth only

indSum_1_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_nteeth),
    seInd_1 = sd(ind_1_nteeth)/sqrt(sum(!is.na(ind_1_nteeth))))

indSum_2_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_nteeth),
    seInd_2 = sd(ind_2_nteeth)/sqrt(sum(!is.na(ind_2_nteeth))))

indSum_3_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_nteeth),
    seInd_3 = sd(ind_3_nteeth)/sqrt(sum(!is.na(ind_3_nteeth))))

indSum_4_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_nteeth),
    seInd_4 = sd(ind_4_nteeth)/sqrt(sum(!is.na(ind_4_nteeth))))

plot_ind1_nteeth <- ggplot(transform(indSum_1_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_1, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
                      geom_point(size = 2)+
                      geom_line(size = 1)+
                      scale_color_viridis(discrete = TRUE, option = "D")+
                      geom_line()+
                      geom_errorbar(width = 0, size = 0.5)+
                      facet_wrap(~clone)+
                      scale_x_continuous(breaks=c(0,0.25,0.5)) +
                      labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                      ylim(0,60) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm")) 

plot_ind2_nteeth <- ggplot(transform(indSum_2_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_2, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
                      geom_point(size = 2)+
                      geom_line(size = 1)+
                      scale_color_viridis(discrete = TRUE, option = "D")+
                      geom_line()+
                      geom_errorbar(width = 0, size = 0.5)+
                      facet_wrap(~clone)+
                      scale_x_continuous(breaks=c(0,0.25,0.5)) +
                      labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                      ylim(0,60) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm")) 

plot_ind3_nteeth <- ggplot(transform(indSum_3_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_3, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
                      geom_point(size = 2)+
                      geom_line(size = 1)+
                      scale_color_viridis(discrete = TRUE, option = "D")+
                      geom_line()+
                      geom_errorbar(width = 0, size = 0.5)+
                      facet_wrap(~clone)+
                      scale_x_continuous(breaks=c(0,0.25,0.5)) +
                      labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                      ylim(0,60) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm")) 

plot_ind4_nteeth <- ggplot(transform(indSum_4_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_4, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
                      geom_point(size = 2)+
                      geom_line(size = 1)+
                      scale_color_viridis(discrete = TRUE, option = "D")+
                      geom_line()+
                      geom_errorbar(width = 0, size = 0.5)+
                      facet_wrap(~clone)+
                      scale_x_continuous(breaks=c(0,0.25,0.5)) +
                      labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                      ylim(0,60) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm")) 


## INDUCTION plot ---------------------------------------------------

patchwork_plots_induction <- (plot_ind1|plot_ind2|plot_ind3|plot_ind4) / 
                                (plot_ind1_ped|plot_ind2_ped|plot_ind3_ped|plot_ind4_ped) / 
                                    (plot_ind1_nteeth|plot_ind2_nteeth|plot_ind3_nteeth|plot_ind4_nteeth)
patchwork_plots_induction

patchwork_plots_induction + 
  plot_annotation(
        title = 'maximal induction @ days 1-4',
        subtitle = '(A-D) overall, (E-H) pedestal only, (I-L) neckteeth only', 
        tag_levels = 'A') + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

ggsave("patchwork_plots_induction.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image


## >> strongest induction in high inducers in instars 1 and 2 (slightly depending on clonal background)
## >> development of neckteeth and pedestal show congruent pattern, i.e. defence morph is not driven by one feature alone


## INDUCTION lmer() ---------------------------------------------------
# maximal induction across all instars (including first instar, w/ induction score == pedestal + nteeth)
glimpse(IndDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(IndDat_use, aes(x = juju, y = maxInd_total))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
ind_mod_full <- lmer(maxInd_total ~ poly(copper,2) + poly(juju,2) + copper:juju + (1|cloneID), data = IndDat_use)

summary(ind_mod_full)
# no evidence for 2nd order copper in coefficient
# no evidence for copper:juju either

Confint(ind_mod_full)  
# confirmed by Confint from car package, there is no evidence to keep the interaction (see final model)

# Type II ANOVA using car with Kenward-Roger df via F test.
Anova(ind_mod_full, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(ind_mod_full, type = "II")
# Pretty much heading towards the conclusion that copper has no 2nd order effect or in interaction

# Reduce the copper polynomial
ind_mod_copper_linear <- lmer(maxInd_total ~ copper + poly(juju,2) + copper:juju + (1|cloneID), data = IndDat_use)

anova(ind_mod_full, ind_mod_copper_linear)  # no difference, i.e. loss of term is justified by likelihood ratio test

summary(ind_mod_copper_linear) 
Confint(ind_mod_copper_linear)

# Include copper just as a single main effect (i.e., exclude interaction term)
ind_mod_copper_linear_noInt <- lmer(maxInd_total ~ copper + poly(juju,2) + (1|cloneID), data = IndDat_use)

anova(ind_mod_copper_linear, ind_mod_copper_linear_noInt) # no difference, i.e. loss of term is justified by likelihood ratio test

summary(ind_mod_copper_linear_noInt)
Confint(ind_mod_copper_linear_noInt)

Anova(ind_mod_copper_linear_noInt)

# NO COPPER model 
ind_mod_no_copper <- lmer(maxInd_total ~ poly(juju,2) + (1|cloneID), data = IndDat_use)
anova(ind_mod_copper_linear_noInt, ind_mod_no_copper)  # linear copper term is better fit - keep! 

pbkrtest::KRmodcomp(ind_mod_copper_linear_noInt, ind_mod_no_copper)
pbkrtest::PBmodcomp(ind_mod_copper_linear_noInt, ind_mod_no_copper)
# tests confirm that linear copper term is better fit - keep! 


# Model to use: linear copper with no interaction effect
Anova(ind_mod_copper_linear_noInt, type = "II", test = "F") # we see that copper linear
Confint(ind_mod_copper_linear_noInt) # suggests even linear effect of copper might not be sig

## for REPORTING, we go back to the full model and the model with linear copper
Anova(ind_mod_full, type = "II", test = "F")
Anova(ind_mod_copper_linear, type = "II", test = "F") # we see that copper linear
summary(ind_mod_copper_linear) # no evidence for interaction
Confint(ind_mod_copper_linear) # suggests even linear effect of copper might not be sig

# pbkrtests as alternative
# pbkrtest::KRmodcomp(ind_mod_full, ind_mod_copper_linear) # is the polynomial copper significant 
# pbkrtest::PBmodcomp(ind_mod_full, ind_mod_copper_linear) # parametric bootstrap

# diagnostics
diagnose <- fortify.merMod(ind_mod_copper_linear_noInt)
d1 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()
# the residual diagnostics are fine

d2 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()+
  facet_wrap(~clone)

d3 <- ggplot(diagnose, aes(sample = .resid))+
  stat_qq_line()+stat_qq()

(d1+d3)/d2


# fixed & random effects
plot_model(ind_mod_copper_linear_noInt, vline.color = "black", type="est", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, sort.est=FALSE)
plot_model(ind_mod_copper_linear_noInt, vline.color = "black", type="re", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, grid = TRUE) #axis.lim=c(-50,50), sort.est=TRUE


# expand data
ind_newX <- expand.grid(
              juju = seq(0,0.5,length = 10),
              copper = seq(0,25, length = 10),  
              cloneID = unique(IndDat$cloneID)
            )

# fixed effect is average among clones & clone effects are the clone-specific 'deviations' from average, their own effects
ind_fixed_effect <- predict(ind_mod_copper_linear_noInt, newdata = ind_newX, re.form = NA)
ind_clone_effect <- predict(ind_mod_copper_linear_noInt, newdata = ind_newX, re.form = ~(1|cloneID))

# housekeeping
ind_plotData <- data.table(ind_newX, fixed_effect=ind_fixed_effect, clone_effect=ind_clone_effect)
ind_plotData

# plot the average and clone specifics
ind_Mod_average <- ggplot(ind_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

ind_Mod_byClone <- ggplot(transform(ind_plotData,
                                    cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                    aes(x = juju, y = copper))+
                          geom_raster(aes(fill = clone_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          facet_wrap(~cloneID)+
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

ind_average <- ggplot(transform(IndDat.ag, 
                                cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                      aes(x = juju, y = meanInd, group = cloneID, colour = cloneID))+
                  geom_pointrange(aes(x = juju, y = meanInd, ymin=meanInd-seInd, ymax=meanInd+seInd, colour = cloneID, group = cloneID), size=0.5, shape = 19)+
                  geom_line(aes(linetype=cloneID),size = 0.7) +
                  scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                  scale_color_viridis(discrete = TRUE, option = "C")+
                  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                  ylim(0,105) +
                  facet_wrap(~copper, nrow = 1) +
                  theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=12, family='Arial'),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm"))


# fix copper to 0 
ind_newX2 <- expand.grid(
                juju = seq(0,0.5,length = 10),
                copper = 0,
                cloneID = unique(IndDat$cloneID)
              )

# get predictions for clone effects at 0 Cu
ind_clone_effect_0Cu <- predict(ind_mod_copper_linear_noInt, newdata = ind_newX2, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
ind_plotData_0Cu <- as.data.table(ind_newX2, ind_clone_effect_0Cu) %>% 
                        mutate(ind_clone_effect_0Cu = ifelse(ind_clone_effect_0Cu < 0, 0, ind_clone_effect_0Cu))

# fix juju to 0 
ind_newX3 <- expand.grid(
                juju = 0,
                copper = seq(0,25,length = 10),
                cloneID = unique(IndDat$cloneID)
              )

# get predictions for clone effects at 0 JuJu
ind_clone_effect_0JuJu <- predict(ind_mod_copper_linear_noInt, newdata = ind_newX3, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
ind_plotData_0JuJu <- as.data.table(ind_newX3, ind_clone_effect_0JuJu) %>% 
                          mutate(ind_clone_effect_0JuJu = ifelse(ind_clone_effect_0JuJu < 0, 0, ind_clone_effect_0JuJu))

# plot
ind_juju_clone <- ggplot(transform(ind_plotData_0Cu, 
                                   cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                        aes(x = juju, y = ind_clone_effect_0Cu, group = cloneID, colour = cloneID))+
                      geom_line(aes(linetype=cloneID),size = 1) + 
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      ylim(0,105) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))

ind_copper_clone <- ggplot(transform(ind_plotData_0JuJu, 
                                     cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = copper, y = ind_clone_effect_0JuJu, group = cloneID, colour = cloneID))+
                      geom_line(aes(linetype=cloneID),size = 1) + 
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(mean~induction[max]), x=expression(copper~(mg~L^{-1}))) +
                      scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                      ylim(0,105) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))


## INDUCTION lmer() plot ---------------------------------------------------

patchwork_plots_induction_model <- ind_average / (ind_Mod_average | ind_Mod_byClone) / (ind_juju_clone | ind_copper_clone | guide_area() )

patchwork_plots_induction_model + 
  plot_annotation(
    title = 'maximal induction @ days 1-4',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_induction_model.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image


## SIZE @ MATURITY lmer() ---------------------------------------------------
glimpse(SizeDat_use)

# plot w/ polynomial order 2 fit in each panel
ggplot(SizeDat_use, aes(x = juju, y = size_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
size_mod_full <- lmer(size_mat ~ poly(copper,2) + poly(juju,2) + copper:juju + (1|cloneID), data = SizeDat_use)

summary(size_mod_full)
Confint(size_mod_full) 
# no evidence for 2nd order copper in coefficient

# Type II ANOVA using car with KR df via F test.
Anova(size_mod_full, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(size_mod_full, type = "II")

# Pretty much heading towards the conclusion that copper has no 2nd order effect

# Reduce the copper polynomial
size_mod_copper_linear <- lmer(size_mat ~ copper + poly(juju,2) + copper:juju + (1|cloneID), data = SizeDat_use)

anova(size_mod_full, size_mod_copper_linear)  # no difference, loss of term is justified by likelihood ratio test

summary(size_mod_copper_linear) 
Confint(size_mod_copper_linear)

# two alternative routes to the same answers
pbkrtest::KRmodcomp(size_mod_full, size_mod_copper_linear) # is the polynomial copper significant 
pbkrtest::PBmodcomp(size_mod_full, size_mod_copper_linear) # parametric bootstrap


# Model to use: linear copper with interaction effect
Anova(size_mod_copper_linear, type = "II", test = "F") # copper, juju and interaction are strongly affecting size @ maturity
# summary(size_mod_copper_linear) 
Confint(size_mod_copper_linear) 

## for REPORTING, we go back to the full model 
Anova(size_mod_full, type = "II", test = "F")
Confint(size_mod_full) 


# diagnostics
diagnose <- fortify.merMod(size_mod_copper_linear)
d1 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()

d2 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()+
  facet_wrap(~clone)

d3 <- ggplot(diagnose, aes(sample = .resid))+
  stat_qq_line()+stat_qq()

(d1+d3)/d2


# fixed & random effects
plot_model(size_mod_copper_linear, vline.color = "black", type="est", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, sort.est=FALSE)
plot_model(size_mod_copper_linear, vline.color = "black", type="re", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, grid = TRUE) #axis.lim=c(-50,50), sort.est=TRUE


# expand data
size_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(SizeDat$cloneID)
)

# fixed effect is average among clones & clone effects are the clone-specific 'deviations' from average, their own effects
size_fixed_effect <- predict(size_mod_copper_linear, newdata = size_newX, re.form = NA)
size_clone_effect <- predict(size_mod_copper_linear, newdata = size_newX, re.form = ~(1|cloneID))

# housekeeping
size_plotData <- data.table(size_newX, fixed_effect=size_fixed_effect, clone_effect=size_clone_effect)
size_plotData


# plot the average and clone specifics
size_Mod_average <- ggplot(size_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

size_Mod_byClone <- ggplot(transform(size_plotData,
                                    cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = juju, y = copper))+
                          geom_raster(aes(fill = clone_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          facet_wrap(~cloneID)+
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

size_average <- ggplot(transform(SizeDat.ag, 
                                cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                    aes(x = juju, y = meanSize, group = cloneID, colour = cloneID))+
                      geom_pointrange(aes(x = juju, y = meanSize, ymin=meanSize-seSize, ymax=meanSize+seSize, colour = cloneID, group = cloneID), size=0.5, shape = 19)+
                      geom_line(aes(linetype=cloneID),size = 0.7) +
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(size[maturity]~(cm)), x=expression(juju~(µl~ml^{-1}))) +
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      ylim(1.3,2.3) +
                      facet_wrap(~copper, nrow = 1) +
                      theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


# fix copper to 0 
size_newX2 <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = 0,
  cloneID = unique(SizeDat$cloneID)
)

# get predictions for clone effects at 0 Cu
size_clone_effect_0Cu <- predict(size_mod_copper_linear, newdata = size_newX2, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
size_plotData_0Cu <- as.data.table(size_newX2, size_clone_effect_0Cu) %>% 
  mutate(size_clone_effect_0Cu = ifelse(size_clone_effect_0Cu < 0, 0, size_clone_effect_0Cu))

# fix juju to 0 
size_newX3 <- expand.grid(
  juju = 0,
  copper = seq(0,25,length = 10),
  cloneID = unique(SizeDat$cloneID)
)

# get predictions for clone effects at 0 JuJu
size_clone_effect_0JuJu <- predict(size_mod_copper_linear, newdata = size_newX3, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
size_plotData_0JuJu <- as.data.table(size_newX3, size_clone_effect_0JuJu) %>% 
  mutate(size_clone_effect_0JuJu = ifelse(size_clone_effect_0JuJu < 0, 0, size_clone_effect_0JuJu))

# plot
size_juju_clone <- ggplot(transform(size_plotData_0Cu, 
                                   cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                        aes(x = juju, y = size_clone_effect_0Cu, group = cloneID, colour = cloneID))+
                        geom_line(aes(linetype=cloneID),size = 1) + 
                        scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                        scale_color_viridis(discrete = TRUE, option = "C")+
                        labs(y=expression(size[maturity]~(cm)), x=expression(juju~(µl~ml^{-1}))) +
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        ylim(1.3,2.3) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))

size_copper_clone <- ggplot(transform(size_plotData_0JuJu, 
                                     cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = copper, y = size_clone_effect_0JuJu, group = cloneID, colour = cloneID))+
                          geom_line(aes(linetype=cloneID),size = 1) + 
                          scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                          scale_color_viridis(discrete = TRUE, option = "C")+
                          labs(y=expression(size[maturity]~(cm)), x=expression(copper~(mg~L^{-1}))) +
                          scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                          ylim(1.3,2.3) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))



## SIZE @ MATURITY lmer() plot ---------------------------------------------------

patchwork_plots_size_model <- size_average / (size_Mod_average | size_Mod_byClone) / (size_juju_clone | size_copper_clone | guide_area() )

patchwork_plots_size_model + 
  plot_annotation(
    title = 'size @ maturity',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_sizeMat_model.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image


## GROWTHRATE MATURITY lmer() ---------------------------------------------------
glimpse(GrowthDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(GrowthDat_use, aes(x = juju, y = growth_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
growth_mod_full <- lmer(growth_mat ~ poly(copper,2) + poly(juju,2) + copper:juju + (1|cloneID), data = GrowthDat_use)

summary(growth_mod_full)
Confint(growth_mod_full) 
# no evidence for 2nd order copper in coefficient
# no evidence for 2nd order juju in coefficient (if any at all)

# Type II ANOVA using car with KR df via F test.
Anova(growth_mod_full, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(growth_mod_full, type = "II")
# Pretty much heading towards the conclusion that copper and juju have no 2nd order effects

# Reduce the copper and juju polynomials
growth_mod_linear <- lmer(growth_mat ~ copper + juju + copper:juju + (1|cloneID), data = GrowthDat_use)

anova(growth_mod_full, growth_mod_linear)  # no difference, i.e. loss of terms is justified by likelihood ratio test

Anova(growth_mod_linear, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in


# Model to use: linear copper and juju with interaction effect
Anova(growth_mod_linear, type = "II", test = "F") # copper is strongly affecting size @ maturity; there is also a strong interaction effect even though juju itself is not affecting the trait
# summary(growth_mod_linear) 
Confint(growth_mod_linear) 

## for REPORTING, we go back to the full model 
Anova(growth_mod_full, type = "II", test = "F")
Confint(growth_mod_full) 


# The residual diagnostics are OK
diagnose <- fortify.merMod(growth_mod_linear)
d1 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()

d2 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()+
  facet_wrap(~clone)

d3 <- ggplot(diagnose, aes(sample = .resid))+
  stat_qq_line()+stat_qq()

(d1+d3)/d2


# fixed & random effects
plot_model(growth_mod_linear, vline.color = "black", type="est", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, sort.est=FALSE)
plot_model(growth_mod_linear, vline.color = "black", type="re", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, grid = TRUE) #axis.lim=c(-50,50), sort.est=TRUE


# expand data
growth_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(GrowthDat$cloneID)
)

# fixed effect is average among clones & clone effects are the clone specific 'deviations' from average, their own effects
growth_fixed_effect <- predict(growth_mod_linear, newdata = growth_newX, re.form = NA)
growth_clone_effect <- predict(growth_mod_linear, newdata = growth_newX, re.form = ~(1|cloneID))

# housekeeping
growth_plotData <- data.frame(growth_newX, fixed_effect=growth_fixed_effect, clone_effect=growth_clone_effect)

# plot the average and clone specifics
growth_Mod_average <- ggplot(growth_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

growth_Mod_byClone <- ggplot(transform(growth_plotData,
                                     cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = juju, y = copper))+
                            geom_raster(aes(fill = clone_effect), interpolate = TRUE)+
                            scale_fill_continuous(type = "viridis")+
                            scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                            labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                            facet_wrap(~cloneID)+
                            theme(rect = element_rect(fill = "transparent"),
                                  panel.grid.major = element_line(colour = "grey70", size=0.25),
                                  panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  plot.background = element_rect(fill = "transparent",colour = NA), 
                                  axis.line = element_line(size = 1),
                                  axis.title.x = element_text(size=12, family='Arial'),
                                  axis.title.y = element_text(size=12, family='Arial'),
                                  axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                  axis.text = element_text(size=12, family='Arial'),
                                  strip.text.x = element_text(size =10, color = "black"),
                                  strip.text.y = element_text(size =10, color = "black"),
                                  panel.spacing.x = unit(4, "mm"),
                                  panel.spacing.y = unit(1, "mm"))

growth_average <- ggplot(transform(GrowthDat.ag, 
                                 cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                      aes(x = juju, y = meanGrowth, group = cloneID, colour = cloneID))+
                      geom_pointrange(aes(x = juju, y = meanGrowth, ymin=meanGrowth-seGrowth, ymax=meanGrowth+seGrowth, colour = cloneID, group = cloneID), size=0.5, shape = 19)+
                      geom_line(aes(linetype=cloneID),size = 0.7) +
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(Growth~rate[maturity]), x=expression(juju~(µl~ml^{-1}))) +
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      ylim(0.05,0.21) +
                      facet_wrap(~copper, nrow = 1) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))


# fix copper to 0 
growth_newX2 <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = 0,
  cloneID = unique(GrowthDat$cloneID)
)

# get predictions for clone effects at 0 Cu
growth_clone_effect_0Cu <- predict(growth_mod_linear, newdata = growth_newX2, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
growth_plotData_0Cu <- as.data.table(growth_newX2, growth_clone_effect_0Cu) %>% 
  mutate(growth_clone_effect_0Cu = ifelse(growth_clone_effect_0Cu < 0, 0, growth_clone_effect_0Cu))

# fix juju to 0 
growth_newX3 <- expand.grid(
  juju = 0,
  copper = seq(0,25,length = 10),
  cloneID = unique(GrowthDat$cloneID)
)

# get predictions for clone effects at 0 JuJu
growth_clone_effect_0JuJu <- predict(growth_mod_linear, newdata = growth_newX3, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
growth_plotData_0JuJu <- as.data.table(growth_newX3, growth_clone_effect_0JuJu) %>% 
  mutate(growth_clone_effect_0JuJu = ifelse(growth_clone_effect_0JuJu < 0, 0, growth_clone_effect_0JuJu))

# plot
growth_juju_clone <- ggplot(transform(growth_plotData_0Cu, 
                                    cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                        aes(x = juju, y = growth_clone_effect_0Cu, group = cloneID, colour = cloneID))+
                          geom_line(aes(linetype=cloneID),size = 1) + 
                          scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                          scale_color_viridis(discrete = TRUE, option = "C")+
                          labs(y=expression(Growth~rate[maturity]), x=expression(juju~(µl~ml^{-1}))) +
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          ylim(0.05,0.21) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))

growth_copper_clone <- ggplot(transform(growth_plotData_0JuJu, 
                                      cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = copper, y = growth_clone_effect_0JuJu, group = cloneID, colour = cloneID))+
                          geom_line(aes(linetype=cloneID),size = 1) + 
                          scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                          scale_color_viridis(discrete = TRUE, option = "C")+
                          labs(y=expression(Growth~rate[maturity]), x=expression(copper~(mg~L^{-1}))) +
                          scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                          ylim(0.05,0.21) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))



## GROWTHRATE MATURITY lmer() plot ---------------------------------------------------

patchwork_plots_growth_model <- growth_average / (growth_Mod_average | growth_Mod_byClone) / (growth_juju_clone | growth_copper_clone | guide_area() )

patchwork_plots_growth_model + 
  plot_annotation(
    title = 'growth rate',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_growth_model.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image



## AGE @ MATURITY lmer() ---------------------------------------------------
glimpse(AgeDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(AgeDat_use, aes(x = juju, y = age_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
age_mod_full <- lmer(age_mat ~ poly(copper,2) + poly(juju,2) + copper:juju + (1|cloneID), data = AgeDat_use)

summary(age_mod_full)
Confint(age_mod_full) 
# no evidence for 2nd order copper in coefficient
# no evidence for 2nd order juju in coefficient (if any at all)

# Type II ANOVA using car with KR df via F test.
Anova(age_mod_full, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(age_mod_full, type = "II")

# Reduce the copper and juju polynomials
age_mod_linear <- lmer(age_mat ~ copper + juju + copper:juju + (1|cloneID), data = AgeDat_use)

anova(age_mod_full, age_mod_linear) # no difference, i.e. loss of these terms is justified by likelihood ratio test

summary(age_mod_linear) 
Confint(age_mod_linear)

Anova(age_mod_linear, type = "II", test = "F") 
# Juju has no effect on age @ maturity but there is a strong copper and copper*juju interaction effect

# two alternative routes to the same answers
pbkrtest::KRmodcomp(age_mod_full, age_mod_linear) # is the polynomial copper significant 
pbkrtest::PBmodcomp(age_mod_full, age_mod_linear) # parametric bootstrap


# Model to use: linear copper and juju with interaction effect
Anova(age_mod_linear, type = "II", test = "F") # copper is strongly affecting size @ maturity; there is also a strong interaction effect even though juju itself is not affecting the trait
# summary(age_mod_linear) 
Confint(age_mod_linear) 

## for REPORTING, we go back to the full model 
Anova(age_mod_full, type = "II", test = "F")
Confint(age_mod_full) 


# diagnostics
diagnose <- fortify.merMod(age_mod_linear)
d1 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()

d2 <- ggplot(diagnose, aes(x = .fitted, y = .resid))+
  geom_point()+
  facet_wrap(~clone)

d3 <- ggplot(diagnose, aes(sample = .resid))+
  stat_qq_line()+stat_qq()

(d1+d3)/d2


# fixed & random effects
plot_model(age_mod_linear, vline.color = "black", type="est", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, sort.est=FALSE)
plot_model(age_mod_linear, vline.color = "black", type="re", show.values=TRUE, show.p=TRUE, value.offset = .3, ci.lvl = 0.9, grid = TRUE) #axis.lim=c(-50,50), sort.est=TRUE

# expand data
age_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(AgeDat$cloneID)
)

# fixed effect is average among clones & clone effects are the clone specific 'deviations' from average, their own effects
age_fixed_effect <- predict(age_mod_linear, newdata = age_newX, re.form = NA)
age_clone_effect <- predict(age_mod_linear, newdata = age_newX, re.form = ~(1|cloneID))

# housekeeping
age_plotData <- data.table(age_newX, fixed_effect=age_fixed_effect, clone_effect=age_clone_effect)

# plot the average and clone specifics
age_Mod_average <- ggplot(age_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                        scale_fill_continuous(type = "viridis")+
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))

age_Mod_byClone <- ggplot(transform(age_plotData,
                                       cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                            aes(x = juju, y = copper))+
                      geom_raster(aes(fill = clone_effect), interpolate = TRUE)+
                      scale_fill_continuous(type = "viridis")+
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                      facet_wrap(~cloneID)+
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))

age_average <- ggplot(transform(AgeDat.ag, 
                                   cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                        aes(x = juju, y = meanAge, group = cloneID, colour = cloneID))+
                    geom_pointrange(aes(x = juju, y = meanAge, ymin=meanAge-seAge, ymax=meanAge+seAge, colour = cloneID, group = cloneID), size=0.5, shape = 19)+
                    geom_line(aes(linetype=cloneID),size = 0.7) +
                    scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                    scale_color_viridis(discrete = TRUE, option = "C")+
                    labs(y=expression(Age[maturity]~(days)), x=expression(juju~(µl~ml^{-1}))) +
                    scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                    ylim(4,13) +
                    facet_wrap(~copper, nrow = 1) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


# fix copper to 0 
age_newX2 <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = 0,
  cloneID = unique(AgeDat$cloneID)
)

# get predictions for clone effects at 0 Cu
age_clone_effect_0Cu <- predict(age_mod_linear, newdata = age_newX2, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
age_plotData_0Cu <- as.data.table(age_newX2, age_clone_effect_0Cu) %>% 
  mutate(age_clone_effect_0Cu = ifelse(age_clone_effect_0Cu < 0, 0, age_clone_effect_0Cu))

# fix juju to 0 
age_newX3 <- expand.grid(
  juju = 0,
  copper = seq(0,25,length = 10),
  cloneID = unique(AgeDat$cloneID)
)

# get predictions for clone effects at 0 JuJu
age_clone_effect_0JuJu <- predict(age_mod_linear, newdata = age_newX3, re.form = ~(1|cloneID) )

# collect and make slightly negative predictions = 0
age_plotData_0JuJu <- as.data.table(age_newX3, age_clone_effect_0JuJu) %>% 
  mutate(age_clone_effect_0JuJu = ifelse(age_clone_effect_0JuJu < 0, 0, age_clone_effect_0JuJu))

# plot
age_juju_clone <- ggplot(transform(age_plotData_0Cu, 
                                      cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                          aes(x = juju, y = age_clone_effect_0Cu, group = cloneID, colour = cloneID))+
                      geom_line(aes(linetype=cloneID),size = 1) + 
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(Age[maturity]~(days)), x=expression(juju~(µl~ml^{-1}))) +
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      ylim(4,13) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))

age_copper_clone <- ggplot(transform(age_plotData_0JuJu, 
                                        cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                                            aes(x = copper, y = age_clone_effect_0JuJu, group = cloneID, colour = cloneID))+
                        geom_line(aes(linetype=cloneID),size = 1) + 
                        scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                        scale_color_viridis(discrete = TRUE, option = "C")+
                        labs(y=expression(Age[maturity]~(days)), x=expression(copper~(mg~L^{-1}))) +
                        scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                        ylim(4,13) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))



## AGE @ MATURITY lmer() plot ---------------------------------------------------

patchwork_plots_age_model <- age_average / (age_Mod_average | age_Mod_byClone) / (age_juju_clone | age_copper_clone | guide_area() )

patchwork_plots_age_model + 
  plot_annotation(
    title = 'age @ maturity',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_age_model.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image



## SURVIVAL ANALYSIS ---------------------------------------------------

survival <- AllData[, c("clone","clone_type","cloneID","juju","copper","rep","survival_time","alive_dead")]
## survival data based on 'maturity data', i.e. data was recorded until animals were mature; 0 = alive, 1 = dead; 

## Cox Hazard Model (use 'coxph' when no random effects, and 'coxme' when using random effects)

## clone as random effect
mod1 <- coxme(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + (1|clone), data = na.omit(survival))
summary(mod1)
confint(mod1)
Anova(mod1, type = "II") # strong copper effect, small juju*copper interaction effect

##clone as frailty term 
mod2 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + frailty(clone), data = na.omit(survival))
summary(mod2)
confint(mod2)
Anova(mod2, type = "II")

anova(mod1, mod2)  # frailty term model slightly better fit


# expand grid
surv_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(survival$cloneID)
)

# fixed effects
surv_fixed_effects_lp <- predict(mod2, type = "lp", newdata = surv_newX, re.form = NA)  # linear prediction
surv_fixed_effects_risk <- predict(mod2, type = "risk", newdata = surv_newX, re.form = NA)  # risk factor 

# housekeeping
surv_fixed_plotData <- data.table(surv_newX, fixed_effect_lp=surv_fixed_effects_lp, fixed_effect_risk=surv_fixed_effects_risk)


## ISSUE (!!)
# How can one extract predicitons from frailty term model (model 1) in order to get clone-specific effects???
# Based on https://stat.ethz.ch/R-manual/R-devel/library/survival/html/predict.coxph.html this might not be possible:
# "Models that contain a frailty term are a special case: due to the technical difficulty, 
# when there is a newdata argument the predictions will always be for a random effect of zero."

# Potential solution - use model w/ copper, juju and interaction term separately for each clone and extract predictions

## clone-specific effects 
# C14
mod3_C14 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "high_C14"]))
plot_model(mod3_C14, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_C14 <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "high_C14")

surv_C14_effects_lp <- predict(mod3_C14, type = "lp", newdata = surv_newX_C14, re.form = NA)  # linear prediction
surv_C14_effects_risk <- predict(mod3_C14, type = "risk", newdata = surv_newX_C14, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_C14 <- data.table(surv_newX_C14, clone_effect_lp=surv_C14_effects_lp, clone_effect_risk=surv_C14_effects_risk)


# Chard
mod3_Chard <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "high_Chard"]))
plot_model(mod3_Chard, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_Chard <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "high_Chard")

surv_Chard_effects_lp <- predict(mod3_Chard, type = "lp", newdata = surv_newX_Chard, re.form = NA)  # linear prediction
surv_Chard_effects_risk <- predict(mod3_Chard, type = "risk", newdata = surv_newX_Chard, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_Chard <- data.table(surv_newX_Chard, clone_effect_lp=surv_Chard_effects_lp, clone_effect_risk=surv_Chard_effects_risk)


# LD33
mod3_LD33 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "high_LD33"]))
plot_model(mod3_LD33, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_LD33 <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "high_LD33")

surv_LD33_effects_lp <- predict(mod3_LD33, type = "lp", newdata = surv_newX_LD33, re.form = NA)  # linear prediction
surv_LD33_effects_risk <- predict(mod3_LD33, type = "risk", newdata = surv_newX_LD33, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_LD33 <- data.table(surv_newX_LD33, clone_effect_lp=surv_LD33_effects_lp, clone_effect_risk=surv_LD33_effects_risk)


# D86A
mod3_D86A <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "low_D86A"]))
plot_model(mod3_D86A, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_D86A <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "low_D86A")

surv_D86A_effects_lp <- predict(mod3_D86A, type = "lp", newdata = surv_newX_D86A, re.form = NA)  # linear prediction
surv_D86A_effects_risk <- predict(mod3_D86A, type = "risk", newdata = surv_newX_D86A, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_D86A <- data.table(surv_newX_D86A, clone_effect_lp=surv_D86A_effects_lp, clone_effect_risk=surv_D86A_effects_risk)


# D87A
mod3_D87A <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "low_D87A"]))
plot_model(mod3_D87A, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_D87A <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "low_D87A")

surv_D87A_effects_lp <- predict(mod3_D87A, type = "lp", newdata = surv_newX_D87A, re.form = NA)  # linear prediction
surv_D87A_effects_risk <- predict(mod3_D87A, type = "risk", newdata = surv_newX_D87A, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_D87A <- data.table(surv_newX_D87A, clone_effect_lp=surv_D87A_effects_lp, clone_effect_risk=surv_D87A_effects_risk)


# Cyril
mod3_Cyril <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[cloneID == "low_Cyril"]))
plot_model(mod3_Cyril, type = "pred", terms = c("juju", "copper"), show.data = F, ci.lvl = .95)  # Predicted values (marginal effects) for specific model terms. 

surv_newX_Cyril <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = "low_Cyril")

surv_Cyril_effects_lp <- predict(mod3_Cyril, type = "lp", newdata = surv_newX_Cyril, re.form = NA)  # linear prediction
surv_Cyril_effects_risk <- predict(mod3_Cyril, type = "risk", newdata = surv_newX_Cyril, re.form = NA)  # risk factor 

# housekeeping
surv_plotData_Cyril <- data.table(surv_newX_Cyril, clone_effect_lp=surv_Cyril_effects_lp, clone_effect_risk=surv_Cyril_effects_risk)


# combine all clone_effects
surv_clone_plotData <- as.data.table(rbind(surv_plotData_C14, surv_plotData_Chard, surv_plotData_LD33,
                                           surv_plotData_D86A, surv_plotData_D87A, surv_plotData_Cyril))


# plot the average and clone specifics
surv_Mod_average <- ggplot(surv_fixed_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = fixed_effect_risk), interpolate = TRUE)+  # fill = surv_fixed_effects_risk
                        scale_fill_continuous(type = "viridis")+
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


surv_Mod_byClone <- ggplot(transform(surv_clone_plotData,
                                     cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                              aes(x = juju, y = copper, fill = clone_effect_risk))+
                        geom_raster(interpolate=T) +  
                        scale_fill_continuous(type = "viridis") +
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                        facet_wrap(~cloneID, scales = "free")+
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


# # alternative plot for clone-specific responses: each panel gets its own colour legend.... (might be handy if we need to plot risk factors - instead of linear predictions)
# transform(surv_clone_plotData,
#                             cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))) %>%
#                   
#                             group_split(cloneID) %>%
#                   
#                             map(
#                               ~ggplot(., aes(juju, copper)) +
#                                 geom_raster(aes(fill = clone_effect_risk), interpolate = TRUE) +
#                                 scale_fill_continuous(type = "viridis") +
#                                 scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
#                                 labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
#                                 facet_grid(~ cloneID) +
#                                 theme(rect = element_rect(fill = "transparent"),
#                                       panel.grid.major = element_line(colour = "grey70", size=0.25),
#                                       panel.grid.minor = element_line(colour = "grey90", size=0.1),
#                                       panel.background = element_rect(fill = "transparent",colour = NA),
#                                       plot.background = element_rect(fill = "transparent",colour = NA),
#                                       axis.line = element_line(size = 1),
#                                       axis.title.x = element_text(size=12, family='Arial'),
#                                       axis.title.y = element_text(size=12, family='Arial'),
#                                       axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
#                                       axis.text = element_text(size=12, family='Arial'),
#                                       strip.text.x = element_text(size =10, color = "black"),
#                                       strip.text.y = element_text(size =10, color = "black"),
#                                       panel.spacing.x = unit(4, "mm"),
#                                       panel.spacing.y = unit(1, "mm"))
#                   
#                             ) %>%
#                             plot_grid()  # can't get this to work: plot_grid(plotlist = ., align = 'hv', ncol = 3)


# fix copper and juju to 0 & plot
surv_juju_clone <- ggplot(transform(surv_clone_plotData[copper == 0], 
                                     cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                          aes(x = juju, y = clone_effect_risk, group = cloneID, colour = cloneID))+
                      geom_line(aes(linetype=cloneID),size = 1) + 
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(survival[~risk~factor]), x=expression(juju~(µl~ml^{-1}))) +
                      scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                      ylim(0,16) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))


surv_copper_clone <- ggplot(transform(surv_clone_plotData[juju == 0], 
                                    cloneID = factor(cloneID, levels=c('high_C14', 'high_Chard', 'high_LD33', 'low_D86A', 'low_D87A', 'low_Cyril'))), 
                          aes(x = copper, y = clone_effect_risk, group = cloneID, colour = cloneID))+
                      geom_line(aes(linetype=cloneID),size = 1) + 
                      scale_linetype_manual(values=c("twodash","twodash","twodash","solid","solid","solid")) +
                      scale_color_viridis(discrete = TRUE, option = "C")+
                      labs(y=expression(survival[~risk~factor]), x=expression(copper~(mg~L^{-1}))) +
                      scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                      ylim(0,16) +
                      theme(rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.line = element_line(size = 1),
                            axis.title.x = element_text(size=12, family='Arial'),
                            axis.title.y = element_text(size=12, family='Arial'),
                            axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                            axis.text = element_text(size=12, family='Arial'),
                            strip.text.x = element_text(size =10, color = "black"),
                            strip.text.y = element_text(size =10, color = "black"),
                            panel.spacing.x = unit(4, "mm"),
                            panel.spacing.y = unit(1, "mm"))



## SURVIVAL ANALYSIS plot ---------------------------------------------------

patchwork_plots_surv_model <- plot_spacer() / (surv_Mod_average | surv_Mod_byClone) / (surv_juju_clone | surv_copper_clone | guide_area() )

patchwork_plots_surv_model + 
  plot_annotation(
    title = 'survival up to maturity',
    subtitle = '(A) average effect of Cu and JuJu among clones, and (B) clone specific effects of Cu and JuJu, (C) effect of JuJu among clones, (D) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(5,6,3), guides = 'collect') & theme(legend.position = "right")  # note re plotting - increased row height due to plot_spacer()

ggsave("patchwork_plots_surv_model.tiff", dpi = 300, device = "tiff")  # Saving 14.3 x 11 in image





#################################
#### CONTRASTING LOW vs HIGH ####
#################################

## INDUCTION lmer() extended ---------------------------------------------------
# maximal induction across all instars (including first instar, w/ induction score == pedestal + nteeth)

# The full RSM model
ind_mod_full <- lm(maxInd_total ~ copper + I(copper^2) + juju + I(juju^2) + copper*juju*clone_type, data = IndDat_use)
Anova(ind_mod_full, type = "II", test = "F") # strong effects: juju, juju^2, clone_type, copper:clone_type, juju:clone_type
Confint(ind_mod_full)  # no evidence for 2nd order copper in coefficient

# Linear Copper (no polynomial)
ind_mod_copper_linear <- lm(maxInd_total ~ copper + juju + I(juju^2) + copper*juju*clone_type, data = IndDat_use)
Anova(ind_mod_copper_linear, type = "II", test = "F") # (strong) effects: copper, juju, juju^2, clone_type, copper:clone_type, juju:clone_type
Confint(ind_mod_copper_linear)  # no evidence for juju:copper in coefficient

anova(ind_mod_full, ind_mod_copper_linear)  # no difference, i.e. loss of term is justified by likelihood ratio test

# Linear Copper and no interaction of copper and juju (but interaction terms with clone_types)
ind_mod_copper_linear_noInt <- lm(maxInd_total ~ copper + juju + I(juju^2) + copper*clone_type + juju*clone_type, data = IndDat_use)
Anova(ind_mod_copper_linear_noInt, type = "II", test = "F") # (strong) effects: copper, juju, juju^2, clone_type, copper:clone_type, juju:clone_type
Confint(ind_mod_copper_linear_noInt) # 

anova(ind_mod_copper_linear, ind_mod_copper_linear_noInt)  # no difference, i.e. loss of term is justified by likelihood ratio test

# No Copper model 
ind_mod_no_copper <- lm(maxInd_total ~ juju + I(juju^2) + juju*clone_type, data = IndDat_use)
Anova(ind_mod_no_copper, type = "II", test = "F") # 
Confint(ind_mod_no_copper) # 

anova(ind_mod_copper_linear_noInt, ind_mod_no_copper)  # loss of term is NOT justified by likelihood ratio test; i.e. use linear copper & no interaction model (as before)

# MODEL to use: ind_mod_copper_linear_noInt
summary(ind_mod_copper_linear_noInt)
Anova(ind_mod_copper_linear_noInt, type = "II") # Type II ANOVA using car with Wald test
Confint(ind_mod_copper_linear_noInt) 


# expand data
ind_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  clone_type = unique(IndDat$clone_type)
)

# fixed effect 
ind_fixed_effect <- predict(ind_mod_copper_linear_noInt, newdata = ind_newX, re.form = "response")

# housekeeping
ind_plotData <- data.table(ind_newX, fixed_effect=ind_fixed_effect)
ind_plotData


# plots
ind_Mod_average <- ggplot(ind_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          facet_wrap(~clone_type) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))


ind_juju_clone <- ggplot(ind_plotData[copper == 0], aes(x = juju, y = fixed_effect, group = clone_type, colour = clone_type))+
                          geom_line(aes(linetype=clone_type),size = 1) + 
                          scale_linetype_manual(values=c("twodash","twodash")) +
                          scale_color_manual(values=c("#000000","#FF0000"))+
                          labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          ylim(0,105) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))


ind_copper_clone <- ggplot(ind_plotData[juju == 0], aes(x = copper, y = fixed_effect, group = clone_type, colour = clone_type))+
                            geom_line(aes(linetype=clone_type),size = 1) +  
                            scale_linetype_manual(values=c("twodash","twodash")) +
                            scale_color_manual(values=c("#000000","#FF0000"))+
                            labs(y=expression(mean~induction[max]), x=expression(copper~(mg~L^{-1}))) +
                            scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                            ylim(0,105) +
                            theme(rect = element_rect(fill = "transparent"),
                                  panel.grid.major = element_line(colour = "grey70", size=0.25),
                                  panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  plot.background = element_rect(fill = "transparent",colour = NA), 
                                  axis.line = element_line(size = 1),
                                  axis.title.x = element_text(size=12, family='Arial'),
                                  axis.title.y = element_text(size=12, family='Arial'),
                                  axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                  axis.text = element_text(size=12, family='Arial'),
                                  strip.text.x = element_text(size =10, color = "black"),
                                  strip.text.y = element_text(size =10, color = "black"),
                                  panel.spacing.x = unit(4, "mm"),
                                  panel.spacing.y = unit(1, "mm"))


## INDUCTION lmer() extended plot ---------------------------------------------------

patchwork_plots_induction_model_HvsL <- ind_Mod_average  / (ind_juju_clone | ind_copper_clone ) 

patchwork_plots_induction_model_HvsL + 
  plot_annotation(
    title = 'maximal induction @ days 1-4',
    subtitle = '(A) average effect of Cu and JuJu among clone types, 
(B) effect of JuJu among clone types, 
(C) effect of Cu among clone types', 
    tag_levels = 'A') + 
  plot_layout(nrow=2, height = c(6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_induction_model_HvsL.tiff", dpi = 300, device = "tiff")  # Saving 6.21 x 6.89 in image



## SIZE @ MATURITY lm() extended ---------------------------------------------------

# The full RSM model
size_mod_full <- lm(size_mat ~ copper + I(copper^2) + juju + I(juju^2) + copper*juju*clone_type, data = SizeDat_use)
Anova(size_mod_full, type = "II", test = "F") # (strong) effects: copper, juju, juju^2, clone_type, copper:juju, copper:clone_type
Confint(size_mod_full) # no evidence for 2nd order copper in coefficient, and limited evidence for copper:juju interaction...

# Linear Copper (no polynomial)
size_mod_copper_linear <- lm(size_mat ~ copper + juju + I(juju^2) + copper*juju*clone_type, data = SizeDat_use)
Anova(size_mod_copper_linear, type = "II", test = "F") # (strong) effects: copper, juju, juju^2, clone_type, copper:juju, copper:clone_type
Confint(size_mod_copper_linear) # no evidence for copper:juju interaction...

anova(size_mod_full, size_mod_copper_linear)  # no difference, i.e. loss of term is justified by likelihood ratio test

# Linear Copper and no interaction of copper and juju (but interaction terms with clone_types)
size_mod_copper_linear_noInt <- lm(size_mat ~ copper + juju + I(juju^2) + copper*clone_type + juju*clone_type, data = SizeDat_use)
Anova(size_mod_copper_linear_noInt, type = "II", test = "F") # (strong) effects: copper, juju, clone_type, copper:clone_type
Confint(size_mod_copper_linear_noInt) # # no evidence for 2nd order juju in coefficient anymore...

anova(size_mod_copper_linear, size_mod_copper_linear_noInt)  # no difference, i.e. loss of term is justified by likelihood ratio test (but marginal so...)

# Linear Copper AND juju, and no interaction of copper and juju (but interaction terms with clone_types)
size_mod_copperANDjuju_linear_noInt <- lm(size_mat ~ copper + juju + copper*clone_type + juju*clone_type, data = SizeDat_use)
Anova(size_mod_copperANDjuju_linear_noInt, type = "II", test = "F") # (strong) effects: copper, juju, clone_type, copper:clone_type
Confint(size_mod_copperANDjuju_linear_noInt) # 

anova(size_mod_copper_linear_noInt, size_mod_copperANDjuju_linear_noInt)  # no difference, i.e. loss of term is justified by likelihood ratio test (but marginal so...)

# MODEL to use: size_mod_copperANDjuju_linear_noInt
summary(size_mod_copperANDjuju_linear_noInt)
Anova(size_mod_copperANDjuju_linear_noInt, type = "II") # Type II ANOVA using car with Wald test
Confint(size_mod_copperANDjuju_linear_noInt) 


# expand data
size_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  clone_type = unique(SizeDat$clone_type))


# fixed effect 
size_fixed_effect <- predict(size_mod_copperANDjuju_linear_noInt, newdata = size_newX, type = "response")

# housekeeping
size_plotData <- data.table(size_newX, fixed_effect=size_fixed_effect)


# plots
size_Mod_average <- ggplot(size_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                          scale_fill_continuous(type = "viridis")+
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                          facet_wrap(~clone_type) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))
                        

size_juju_clone <- ggplot(size_plotData[copper == 0], aes(x = juju, y = fixed_effect, group = clone_type, colour = clone_type))+
                          geom_line(aes(linetype=clone_type),size = 1) + 
                          scale_linetype_manual(values=c("twodash","twodash")) +
                          scale_color_manual(values=c("#000000","#FF0000"))+
                          labs(y=expression(size[maturity]), x=expression(juju~(µl~ml^{-1}))) +
                          scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                          ylim(1.3,2.3) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))


size_copper_clone <- ggplot(size_plotData[juju == 0], aes(x = copper, y = fixed_effect, group = clone_type, colour = clone_type))+
                          geom_line(aes(linetype=clone_type),size = 1) +  
                          scale_linetype_manual(values=c("twodash","twodash")) +
                          scale_color_manual(values=c("#000000","#FF0000"))+
                          labs(y=expression(size[maturity]), x=expression(copper~(mg~L^{-1}))) +
                          scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                          ylim(1.3,2.3) +
                          theme(rect = element_rect(fill = "transparent"),
                                panel.grid.major = element_line(colour = "grey70", size=0.25),
                                panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA), 
                                axis.line = element_line(size = 1),
                                axis.title.x = element_text(size=12, family='Arial'),
                                axis.title.y = element_text(size=12, family='Arial'),
                                axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                axis.text = element_text(size=12, family='Arial'),
                                strip.text.x = element_text(size =10, color = "black"),
                                strip.text.y = element_text(size =10, color = "black"),
                                panel.spacing.x = unit(4, "mm"),
                                panel.spacing.y = unit(1, "mm"))


## SIZE @ MATURITY lmer() extended plot ---------------------------------------------------

patchwork_plots_size_model_HvsL <- size_Mod_average  / (size_juju_clone | size_copper_clone ) 

patchwork_plots_size_model_HvsL + 
  plot_annotation(
    title = 'size @ maturity',
    subtitle = '(A) average effect of Cu and JuJu among clone types, 
    (B) effect of JuJu among clone types, 
    (C) effect of Cu among clone types', 
    tag_levels = 'A') + 
  plot_layout(nrow=2, height = c(6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_sizeMat_model_HvsL.tiff", dpi = 300, device = "tiff")  # Saving 6.21 x 6.89 in image



## GROWTHRATE @ MATURITY lm() extended ---------------------------------------------------
glimpse(GrowthDat_use)

# The full RSM model
growth_mod_full <- lm(growth_mat ~ copper + I(copper^2) + juju + I(juju^2) + copper*juju*clone_type, data = GrowthDat_use)
Anova(growth_mod_full, type = "II", test = "F") # strong effects: copper, copper:juju, copper:clone_type
Confint(growth_mod_full) # no evidence for 2nd order copper and juju in coefficient, limited evidence of copper:juju and copper:clone_type

# Linear Copper and juju (no polynomials)
growth_mod_copperANDjuju_linear <- lm(growth_mat ~ copper + juju + copper*juju*clone_type, data = GrowthDat_use)
Anova(growth_mod_copperANDjuju_linear, type = "II", test = "F") # strong effects: copper, copper:juju, copper:clone_type
Confint(growth_mod_copperANDjuju_linear) # no evidence of copper:juju and copper:clone_type

# Linear Copper AND juju, and no interaction of copper and juju (but interaction terms with clone_types)
growth_mod_copperANDjuju_linear_noInt <- lm(growth_mat ~ copper + juju + copper*clone_type + juju*clone_type, data = GrowthDat_use)
Anova(growth_mod_copperANDjuju_linear_noInt, type = "II", test = "F") # strong effects: copper, copper:clone_type
Confint(growth_mod_copperANDjuju_linear_noInt) # 


# expand data
growth_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  clone_type = unique(GrowthDat_use$clone_type))


# fixed effect 
growth_fixed_effect <- predict(growth_mod_copperANDjuju_linear_noInt, newdata = growth_newX, re.form = "response")

# housekeeping
growth_plotData <- data.table(growth_newX, fixed_effect=growth_fixed_effect)


# plots
growth_Mod_average <- ggplot(growth_plotData, aes(x = juju, y = copper))+
                              geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                              scale_fill_continuous(type = "viridis")+
                              scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                              labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                              facet_wrap(~clone_type) +
                              theme(rect = element_rect(fill = "transparent"),
                                    panel.grid.major = element_line(colour = "grey70", size=0.25),
                                    panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                    panel.background = element_rect(fill = "transparent",colour = NA),
                                    plot.background = element_rect(fill = "transparent",colour = NA), 
                                    axis.line = element_line(size = 1),
                                    axis.title.x = element_text(size=12, family='Arial'),
                                    axis.title.y = element_text(size=12, family='Arial'),
                                    axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                    axis.text = element_text(size=12, family='Arial'),
                                    strip.text.x = element_text(size =10, color = "black"),
                                    strip.text.y = element_text(size =10, color = "black"),
                                    panel.spacing.x = unit(4, "mm"),
                                    panel.spacing.y = unit(1, "mm"))


growth_juju_clone <- ggplot(growth_plotData[copper == 0], aes(x = juju, y = fixed_effect, group = clone_type, colour = clone_type))+
                              geom_line(aes(linetype=clone_type),size = 1) + 
                              scale_linetype_manual(values=c("twodash","twodash")) +
                              scale_color_manual(values=c("#000000","#FF0000"))+
                              labs(y=expression(Growth~rate[maturity]), x=expression(juju~(µl~ml^{-1}))) +
                              scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                              ylim(0.05,0.21) +
                              theme(rect = element_rect(fill = "transparent"),
                                    panel.grid.major = element_line(colour = "grey70", size=0.25),
                                    panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                    panel.background = element_rect(fill = "transparent",colour = NA),
                                    plot.background = element_rect(fill = "transparent",colour = NA), 
                                    axis.line = element_line(size = 1),
                                    axis.title.x = element_text(size=12, family='Arial'),
                                    axis.title.y = element_text(size=12, family='Arial'),
                                    axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                    axis.text = element_text(size=12, family='Arial'),
                                    strip.text.x = element_text(size =10, color = "black"),
                                    strip.text.y = element_text(size =10, color = "black"),
                                    panel.spacing.x = unit(4, "mm"),
                                    panel.spacing.y = unit(1, "mm"))


growth_copper_clone <- ggplot(growth_plotData[juju == 0], aes(x = copper, y = fixed_effect, group = clone_type, colour = clone_type))+
                              geom_line(aes(linetype=clone_type),size = 1) +  
                              scale_linetype_manual(values=c("twodash","twodash")) +
                              scale_color_manual(values=c("#000000","#FF0000"))+
                              labs(y=expression(Growth~rate[maturity]), x=expression(copper~(mg~L^{-1}))) +
                              scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                              ylim(0.05,0.21) +
                              theme(rect = element_rect(fill = "transparent"),
                                    panel.grid.major = element_line(colour = "grey70", size=0.25),
                                    panel.grid.minor = element_line(colour = "grey90", size=0.1),
                                    panel.background = element_rect(fill = "transparent",colour = NA),
                                    plot.background = element_rect(fill = "transparent",colour = NA), 
                                    axis.line = element_line(size = 1),
                                    axis.title.x = element_text(size=12, family='Arial'),
                                    axis.title.y = element_text(size=12, family='Arial'),
                                    axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                                    axis.text = element_text(size=12, family='Arial'),
                                    strip.text.x = element_text(size =10, color = "black"),
                                    strip.text.y = element_text(size =10, color = "black"),
                                    panel.spacing.x = unit(4, "mm"),
                                    panel.spacing.y = unit(1, "mm"))


## GROWTHRATE MATURITY lmer() extended plot ---------------------------------------------------

patchwork_plots_growth_model_HvsL <- growth_Mod_average  / (growth_juju_clone | growth_copper_clone ) 

patchwork_plots_growth_model_HvsL + 
  plot_annotation(
    title = 'growth rate',
    subtitle = '(A) average effect of Cu and JuJu among clone types, 
    (B) effect of JuJu among clone types, 
    (C) effect of Cu among clone types', 
    tag_levels = 'A') + 
  plot_layout(nrow=2, height = c(6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_growth_model_HvsL.tiff", dpi = 300, device = "tiff")  # Saving 6.21 x 6.89 in image



## AGE @ MATURITY lm() extended ---------------------------------------------------

# The full RSM model
age_mod_full <- lm(age_mat ~ copper + I(copper^2) + juju + I(juju^2) + copper*juju*clone_type, data = AgeDat_use)
Anova(age_mod_full, type = "II", test = "F") # (strong) effects: copper, copper:juju, copper:clone_type, copper:juju:clone_type
Confint(age_mod_full) # no evidence for 2nd order copper and juju in coefficient, limited evidence of copper:juju and copper:clone_type

# Linear Copper and juju (no polynomials)
age_mod_copperANDjuju_linear <- lm(age_mat ~ copper + juju + copper*juju*clone_type, data = AgeDat_use)
Anova(age_mod_copperANDjuju_linear, type = "II", test = "F") # (strong) effects: copper, copper:juju, copper:clone_type, copper:juju:clone_type
Confint(age_mod_copperANDjuju_linear) # no evidence of copper:juju and copper:clone_type

anova(age_mod_full, age_mod_copperANDjuju_linear)  # no difference, i.e. loss of terms is justified by likelihood ratio test

# Linear Copper and juju, and no interaction of copper and juju (but interaction terms with clone_types)
age_mod_copperANDjuju_linear_noInt <- lm(age_mat ~ copper + juju + copper*clone_type + juju*clone_type, data = AgeDat_use)
Anova(age_mod_copperANDjuju_linear_noInt, type = "II", test = "F") # (strong) effects: copper, copper:clone_type
Confint(age_mod_copperANDjuju_linear_noInt) #

anova(age_mod_copperANDjuju_linear, age_mod_copperANDjuju_linear_noInt)  # loss of term is NOT justified by likelihood ratio test; i.e. use model with linear copper and juju terms

# MODEL to use: age_mod_copperANDjuju_linear
summary(age_mod_copperANDjuju_linear)
Anova(age_mod_copperANDjuju_linear, type = "II") # Type II ANOVA using car with Wald test
Confint(age_mod_copperANDjuju_linear) 


# expand data
age_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  clone_type = unique(AgeDat$clone_type))

# fixed effect 
age_fixed_effect <- predict(age_mod_copperANDjuju_linear, newdata = age_newX, re.form = "response")

# housekeeping
age_plotData <- data.table(age_newX, fixed_effect=age_fixed_effect)


# plots
age_Mod_average <- ggplot(age_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = fixed_effect), interpolate = TRUE)+
                        scale_fill_continuous(type = "viridis")+
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
                        facet_wrap(~clone_type) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


age_juju_clone <- ggplot(age_plotData[copper == 0], aes(x = juju, y = fixed_effect, group = clone_type, colour = clone_type))+
                        geom_line(aes(linetype=clone_type),size = 1) + 
                        scale_linetype_manual(values=c("twodash","twodash")) +
                        scale_color_manual(values=c("#000000","#FF0000"))+
                        labs(y=expression(age[maturity]), x=expression(juju~(µl~ml^{-1}))) +
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        ylim(4,13) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


age_copper_clone <- ggplot(age_plotData[juju == 0], aes(x = copper, y = fixed_effect, group = clone_type, colour = clone_type))+
                        geom_line(aes(linetype=clone_type),size = 1) +  
                        scale_linetype_manual(values=c("twodash","twodash")) +
                        scale_color_manual(values=c("#000000","#FF0000"))+
                        labs(y=expression(age[maturity]), x=expression(copper~(mg~L^{-1}))) +
                        scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                        ylim(4,13) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


## AGE @ MATURITY lmer() extended plot ---------------------------------------------------

patchwork_plots_age_model_HvsL <- age_Mod_average / (age_juju_clone | age_copper_clone )

patchwork_plots_age_model_HvsL + 
  plot_annotation(
    title = 'age @ maturity',
    subtitle = '(A) average effect of Cu and JuJu among clone types, 
    (B) effect of JuJu among clone types, 
    (C) effect of Cu among clone types', 
    tag_levels = 'A') + 
  plot_layout(nrow=2, height = c(6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_age_model_HvsL.tiff", dpi = 300, device = "tiff")  # Saving 6.21 x 6.9 in image



## SURVIVAL ANALYSIS extended  ---------------------------------------------------

survival <- AllData[, c("clone","clone_type","cloneID","juju","copper","rep","survival_time","alive_dead")]
## survival data based on 'maturity data', i.e. data was recorded until animals were mature; 0 = alive, 1 = dead; 

## Cox Hazard Model (use 'coxph' when no random effects, and 'coxme' when using random effects)

#
mod1_HvsL <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper*clone_type, data = na.omit(survival))
summary(mod1_HvsL)
confint(mod1_HvsL)
Anova(mod1_HvsL, type = "II") # strong copper and copper*clone effect, small juju*copper interaction 


# expand grid
surv_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  clone_type = unique(survival$clone_type)
)

# fixed effects
surv_fixed_effects_HvsL_lp <- predict(mod1_HvsL, type = "lp", newdata = surv_newX, re.form = NA)  # linear prediction
surv_fixed_effects_HvsL_risk <- predict(mod1_HvsL, type = "risk", newdata = surv_newX, re.form = NA)  # risk factor 

# housekeeping
surv_plotData <- data.table(surv_newX, fixed_effect_lp=surv_fixed_effects_HvsL_lp, fixed_effect_risk=surv_fixed_effects_HvsL_risk)


# plots
surv_Mod_average <- ggplot(surv_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = fixed_effect_risk), interpolate = TRUE)+  # fill = surv_fixed_effects_risk
                        scale_fill_continuous(type = "viridis")+
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) + 
                        facet_wrap(~clone_type) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))
                      

surv_juju_clone <- ggplot(surv_plotData[copper == 0], aes(x = juju, y = fixed_effect_risk, group = clone_type, colour = clone_type))+
                        geom_line(aes(linetype=clone_type),size = 1) + 
                        scale_linetype_manual(values=c("twodash","twodash")) +
                        scale_color_manual(values=c("#000000","#FF0000"))+
                        labs(y=expression(survival[~risk~factor]), x=expression(juju~(µl~ml^{-1}))) +
                        scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
                        ylim(0,16) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


surv_copper_clone <- ggplot(surv_plotData[juju == 0], aes(x = copper, y = fixed_effect_risk, group = clone_type, colour = clone_type))+
                        geom_line(aes(linetype=clone_type),size = 1) + 
                        scale_linetype_manual(values=c("twodash","twodash")) +
                        scale_color_manual(values=c("#000000","#FF0000"))+
                        labs(y=expression(survival[~risk~factor]), x=expression(copper~(mg~L^{-1}))) +
                        scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
                        ylim(0,16) +
                        theme(rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))



## SURVIVAL ANALYSIS plot ---------------------------------------------------

patchwork_plots_surv_model_HvsL <- surv_Mod_average / (surv_juju_clone | surv_copper_clone )

patchwork_plots_surv_model_HvsL + 
  plot_annotation(
    title = 'survival up to maturity',
    subtitle = '(A) average effect of Cu and JuJu among clone types, 
    (B) effect of JuJu among clone types, 
    (C) effect of Cu among clone types', 
    tag_levels = 'A') + 
  plot_layout(nrow=2, height = c(6,3), guides = 'collect') & theme(legend.position = "right")

ggsave("patchwork_plots_surv_model_HvsL.tiff", dpi = 300, device = "tiff")  # Saving 6.21 x 6.9 in image





## Induction vs other life history traits (i.e. age and size at maturity and somatic growth rate) --------------------------

# aggregate data on LOW and HIGH and calculate SE among clones for morpho induction and LH traits (use for 'error graphs')

IndDat_use_LOWvsHIGH <- as.data.table(na.omit(IndDat_use) %>% 
  group_by(copper, juju, clone_type) %>% 
  summarise(
    meanInd = mean(maxInd_total), 
    seInd = sd(maxInd_total)/sqrt(sum(!is.na(maxInd_total)))))

IndDat_use_LOWvsHIGH[, Ind_min := meanInd-seInd]
IndDat_use_LOWvsHIGH[, Ind_max := meanInd+seInd]

setkey(IndDat_use_LOWvsHIGH, copper, juju, clone_type)



SizeDat_use_LOWvsHIGH <- as.data.table(na.omit(SizeDat_use) %>% 
  group_by(copper, juju, clone_type) %>% 
  summarise(
    meanSize = mean(size_mat), 
    seSize = sd(size_mat)/sqrt(sum(!is.na(size_mat)))))

SizeDat_use_LOWvsHIGH[, Size_min := meanSize-seSize]
SizeDat_use_LOWvsHIGH[, Size_max := meanSize+seSize]

setkey(SizeDat_use_LOWvsHIGH, copper, juju, clone_type)



GrowthDat_use_LOWvsHIGH <- as.data.table(na.omit(GrowthDat_use) %>% 
  group_by(copper, juju, clone_type) %>% 
  summarise(
    meanGrowth = mean(growth_mat), 
    seGrowth = sd(growth_mat)/sqrt(sum(!is.na(growth_mat)))))

GrowthDat_use_LOWvsHIGH[, Growth_min := meanGrowth-seGrowth]
GrowthDat_use_LOWvsHIGH[, Growth_max := meanGrowth+seGrowth]

setkey(GrowthDat_use_LOWvsHIGH, copper, juju, clone_type)



AgeDat_use_LOWvsHIGH <- as.data.table(na.omit(AgeDat_use) %>% 
  group_by(copper, juju, clone_type) %>% 
  summarise(
    meanAge = mean(age_mat), 
    seAge = sd(age_mat)/sqrt(sum(!is.na(age_mat)))))

AgeDat_use_LOWvsHIGH[, Age_min := meanAge-seAge]
AgeDat_use_LOWvsHIGH[, Age_max := meanAge+seAge]

setkey(AgeDat_use_LOWvsHIGH, copper, juju, clone_type)



Ind_Size <- merge(IndDat_use_LOWvsHIGH, SizeDat_use_LOWvsHIGH)
Ind_Size[, copper_juju := c('Cu0 - Juju0', 
                            'Cu0 - Juju0', 
                            'Cu0 - Juju0.1',
                            'Cu0 - Juju0.1',
                            'Cu0 - Juju0.25',
                            'Cu0 - Juju0.25',
                            'Cu0 - Juju0.5',
                            'Cu0 - Juju0.5',
                           
                            'Cu5 - Juju0', 
                            'Cu5 - Juju0', 
                            'Cu5 - Juju0.1',
                            'Cu5 - Juju0.1',
                            'Cu5 - Juju0.25',
                            'Cu5 - Juju0.25',
                            'Cu5 - Juju0.5',
                            'Cu5 - Juju0.5',
                            
                            'Cu10 - Juju0', 
                            'Cu10 - Juju0', 
                            'Cu10 - Juju0.1',
                            'Cu10 - Juju0.1',
                            'Cu10 - Juju0.25',
                            'Cu10 - Juju0.25',
                            'Cu10 - Juju0.5',
                            'Cu10 - Juju0.5',
                            
                            'Cu25 - Juju0', 
                            'Cu25 - Juju0', 
                            'Cu25 - Juju0.1',
                            'Cu25 - Juju0.1',
                            'Cu25 - Juju0.25',
                            'Cu25 - Juju0.25',
                            'Cu25 - Juju0.5',
                            'Cu25 - Juju0.5')]

Ind_Growth <- merge(IndDat_use_LOWvsHIGH, GrowthDat_use_LOWvsHIGH)
Ind_Growth[, copper_juju := c('Cu0 - Juju0', 
                            'Cu0 - Juju0', 
                            'Cu0 - Juju0.1',
                            'Cu0 - Juju0.1',
                            'Cu0 - Juju0.25',
                            'Cu0 - Juju0.25',
                            'Cu0 - Juju0.5',
                            'Cu0 - Juju0.5',
                            
                            'Cu5 - Juju0', 
                            'Cu5 - Juju0', 
                            'Cu5 - Juju0.1',
                            'Cu5 - Juju0.1',
                            'Cu5 - Juju0.25',
                            'Cu5 - Juju0.25',
                            'Cu5 - Juju0.5',
                            'Cu5 - Juju0.5',
                            
                            'Cu10 - Juju0', 
                            'Cu10 - Juju0', 
                            'Cu10 - Juju0.1',
                            'Cu10 - Juju0.1',
                            'Cu10 - Juju0.25',
                            'Cu10 - Juju0.25',
                            'Cu10 - Juju0.5',
                            'Cu10 - Juju0.5',
                            
                            'Cu25 - Juju0', 
                            'Cu25 - Juju0', 
                            'Cu25 - Juju0.1',
                            'Cu25 - Juju0.1',
                            'Cu25 - Juju0.25',
                            'Cu25 - Juju0.25',
                            'Cu25 - Juju0.5',
                            'Cu25 - Juju0.5')]


Ind_Age <- merge(IndDat_use_LOWvsHIGH, AgeDat_use_LOWvsHIGH)
Ind_Age[, copper_juju := c('Cu0 - Juju0', 
                              'Cu0 - Juju0', 
                              'Cu0 - Juju0.1',
                              'Cu0 - Juju0.1',
                              'Cu0 - Juju0.25',
                              'Cu0 - Juju0.25',
                              'Cu0 - Juju0.5',
                              'Cu0 - Juju0.5',
                              
                              'Cu5 - Juju0', 
                              'Cu5 - Juju0', 
                              'Cu5 - Juju0.1',
                              'Cu5 - Juju0.1',
                              'Cu5 - Juju0.25',
                              'Cu5 - Juju0.25',
                              'Cu5 - Juju0.5',
                              'Cu5 - Juju0.5',
                              
                              'Cu10 - Juju0', 
                              'Cu10 - Juju0', 
                              'Cu10 - Juju0.1',
                              'Cu10 - Juju0.1',
                              'Cu10 - Juju0.25',
                              'Cu10 - Juju0.25',
                              'Cu10 - Juju0.5',
                              'Cu10 - Juju0.5',
                              
                              'Cu25 - Juju0', 
                              'Cu25 - Juju0', 
                              'Cu25 - Juju0.1',
                              'Cu25 - Juju0.1',
                              'Cu25 - Juju0.25',
                              'Cu25 - Juju0.25',
                              'Cu25 - Juju0.5',
                              'Cu25 - Juju0.5')]


## Induction vs other life history traits - plots --------------------------

# plot extremes (i.e., [copper == 0][juju == 0], [copper == 0][juju == 0.5], [copper == 25][juju == 0], [copper == 25][juju == 0])

# Induction vs Size
Ind_size_p1 <- ggplot(Ind_Size[copper == 0][juju == 0][, -"clone_type"], aes(x = meanInd, y = meanSize)) + 
                    geom_line(linetype="twodash", colour="darkgrey") + 
                    geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                    geom_errorbar(aes(ymin = Size_min, ymax = Size_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    ylim(1.5, 2.05) +
                    xlim(0, 100) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


Ind_size_p2 <- ggplot(Ind_Size[copper == 0][juju == 0.5][, -"clone_type"], aes(x = meanInd, y = meanSize)) + 
                    geom_line(linetype="twodash", colour="darkgrey") + 
                    geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                    geom_errorbar(aes(ymin = Size_min, ymax = Size_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    ylim(1.5, 2.05) +
                    xlim(0, 100) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


Ind_size_p3 <- ggplot(Ind_Size[copper == 25][juju == 0][, -"clone_type"], aes(x = meanInd, y = meanSize)) + 
                    geom_line(linetype="twodash", colour="darkgrey") + 
                    geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                    geom_errorbar(aes(ymin = Size_min, ymax = Size_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    ylim(1.5, 2.05) +
                    xlim(0, 100) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


Ind_size_p4 <- ggplot(Ind_Size[copper == 25][juju == 0.5][, -"clone_type"], aes(x = meanInd, y = meanSize)) + 
                    geom_line(linetype="twodash", colour="darkgrey") + 
                    geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                    geom_errorbar(aes(ymin = Size_min, ymax = Size_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                    ylim(1.5, 2.05) +
                    xlim(0, 100) +
                    theme(rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=12, family='Arial'),
                          axis.title.y = element_text(size=12, family='Arial'),
                          axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                          axis.text = element_text(size=12, family='Arial'),
                          strip.text.x = element_text(size =10, color = "black"),
                          strip.text.y = element_text(size =10, color = "black"),
                          panel.spacing.x = unit(4, "mm"),
                          panel.spacing.y = unit(1, "mm"))


(Ind_size_p3 / Ind_size_p1) | (Ind_size_p4 / Ind_size_p2)



# Induction vs Age
Ind_age_p1 <- ggplot(Ind_Age[copper == 0][juju == 0][, -"clone_type"], aes(x = meanInd, y = meanAge)) + 
                  geom_line(linetype="twodash", colour="darkgrey") + 
                  geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                  geom_errorbar(aes(ymin = Age_min, ymax = Age_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                  geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                  ylim(1.5, 2.05) +
                  xlim(0, 100) +
                  theme(rect = element_rect(fill = "transparent"),
                        panel.grid.major = element_line(colour = "grey70", size=0.25),
                        panel.grid.minor = element_line(colour = "grey90", size=0.1),
                        panel.background = element_rect(fill = "transparent",colour = NA),
                        plot.background = element_rect(fill = "transparent",colour = NA), 
                        axis.line = element_line(size = 1),
                        axis.title.x = element_text(size=12, family='Arial'),
                        axis.title.y = element_text(size=12, family='Arial'),
                        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                        axis.text = element_text(size=12, family='Arial'),
                        strip.text.x = element_text(size =10, color = "black"),
                        strip.text.y = element_text(size =10, color = "black"),
                        panel.spacing.x = unit(4, "mm"),
                        panel.spacing.y = unit(1, "mm"))


Ind_age_p2 <- ggplot(Ind_Age[copper == 0][juju == 0.5][, -"clone_type"], aes(x = meanInd, y = meanAge)) + 
                  geom_line(linetype="twodash", colour="darkgrey") + 
                  geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                  geom_errorbar(aes(ymin = Age_min, ymax = Age_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                  geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                  ylim(1.5, 2.05) +
                  xlim(0, 100) +
                  theme(rect = element_rect(fill = "transparent"),
                        panel.grid.major = element_line(colour = "grey70", size=0.25),
                        panel.grid.minor = element_line(colour = "grey90", size=0.1),
                        panel.background = element_rect(fill = "transparent",colour = NA),
                        plot.background = element_rect(fill = "transparent",colour = NA), 
                        axis.line = element_line(size = 1),
                        axis.title.x = element_text(size=12, family='Arial'),
                        axis.title.y = element_text(size=12, family='Arial'),
                        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                        axis.text = element_text(size=12, family='Arial'),
                        strip.text.x = element_text(size =10, color = "black"),
                        strip.text.y = element_text(size =10, color = "black"),
                        panel.spacing.x = unit(4, "mm"),
                        panel.spacing.y = unit(1, "mm"))


Ind_age_p3 <- ggplot(Ind_Age[copper == 25][juju == 0][, -"clone_type"], aes(x = meanInd, y = meanAge)) + 
                geom_line(linetype="twodash", colour="darkgrey") + 
                geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                geom_errorbar(aes(ymin = Age_min, ymax = Age_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                ylim(1.5, 2.05) +
                xlim(0, 100) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=12, family='Arial'),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm"))


Ind_age_p4 <- ggplot(Ind_Age[copper == 25][juju == 0.5][, -"clone_type"], aes(x = meanInd, y = meanAge)) + 
                geom_line(linetype="twodash", colour="darkgrey") + 
                geom_point(colour = c("#FF0000", "#000000"), size = 3) + 
                geom_errorbar(aes(ymin = Age_min, ymax = Age_max), width = 1, colour = c("#FF0000", "#000000"), size = 0.2) + 
                geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max), height = 0.007, colour = c("#FF0000", "#000000"), size = 0.2) + 
                ylim(1.5, 2.05) +
                xlim(0, 100) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=12, family='Arial'),
                      axis.title.y = element_text(size=12, family='Arial'),
                      axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                      axis.text = element_text(size=12, family='Arial'),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm"))


(Ind_age_p3 / Ind_age_p1) | (Ind_age_p4 / Ind_age_p2)

########








Ind_Size_plot <- ggplot(transform(Ind_Size,
                        copper_juju = factor(copper_juju, levels=c('Cu0 - Juju0', 'Cu5 - Juju0', 'Cu10 - Juju0', 'Cu25 - Juju0',
                                                                   'Cu0 - Juju0.1', 'Cu5 - Juju0.1', 'Cu10 - Juju0.1', 'Cu25 - Juju0.1',
                                                                   'Cu0 - Juju0.25', 'Cu5 - Juju0.25', 'Cu10 - Juju0.25', 'Cu25 - Juju0.25',
                                                                   'Cu0 - Juju0.5', 'Cu5 - Juju0.5', 'Cu10 - Juju0.5', 'Cu25 - Juju0.5'))),
                        aes(x = meanInd,y = meanSize, colour = as.factor(clone_type))) + 
                        geom_point(size = 3, pch = 21) + 
                        geom_errorbar(aes(ymin = Size_min, ymax = Size_max)) + 
                        geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max)) + 
                        facet_wrap(~copper_juju) + 
                        theme(legend.position="top", , 
                              rect = element_rect(fill = "transparent"),
                              panel.grid.major = element_line(colour = "grey70", size=0.25),
                              panel.grid.minor = element_line(colour = "grey90", size=0.1),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              #strip.text.x = element_blank(),
                              #axis.text.x = element_blank(), 
                              #axis.title.x = element_blank(), 
                              #axis.title.y = element_blank(),
                              axis.line = element_line(size = 1),
                              axis.title.x = element_text(size=12, family='Arial'),
                              axis.title.y = element_text(size=12, family='Arial'),
                              axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                              axis.text = element_text(size=12, family='Arial'),
                              strip.text.x = element_text(size =10, color = "black"),
                              strip.text.y = element_text(size =10, color = "black"),
                              panel.spacing.x = unit(4, "mm"),
                              panel.spacing.y = unit(1, "mm"))


tiff(file = "Ind_Size_plot.tiff", width = 1000, height = 800)
Ind_Size_plot  + labs(x = "mean induction score", y = "mean size @ maturity") + scale_x_continuous(limits=c(0,100), breaks=c(0,25,50,75,100)) + scale_color_manual(values=c("#FF0000","#000000"))
dev.off()



Ind_Age_plot <- ggplot(transform(Ind_Age,
                                  copper_juju = factor(copper_juju, levels=c('Cu0 - Juju0', 'Cu5 - Juju0', 'Cu10 - Juju0', 'Cu25 - Juju0',
                                                                             'Cu0 - Juju0.1', 'Cu5 - Juju0.1', 'Cu10 - Juju0.1', 'Cu25 - Juju0.1',
                                                                             'Cu0 - Juju0.25', 'Cu5 - Juju0.25', 'Cu10 - Juju0.25', 'Cu25 - Juju0.25',
                                                                             'Cu0 - Juju0.5', 'Cu5 - Juju0.5', 'Cu10 - Juju0.5', 'Cu25 - Juju0.5'))),
                        aes(x = meanInd,y = meanAge, colour = as.factor(clone_type))) + 
  geom_point(size = 3, pch = 21) + 
  geom_errorbar(aes(ymin = Age_min, ymax = Age_max)) + 
  geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max)) + 
  facet_wrap(~copper_juju) + 
  theme(legend.position="top", , 
        rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        #strip.text.x = element_blank(),
        #axis.text.x = element_blank(), 
        #axis.title.x = element_blank(), 
        #axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))

tiff(file = "Ind_Age_plot.tiff", width = 1000, height = 800)
Ind_Age_plot  + labs(x = "mean induction score", y = "mean age @ maturity") + scale_x_continuous(limits=c(0,100), breaks=c(0,25,50,75,100)) + scale_color_manual(values=c("#FF0000","#000000"))
dev.off()


Ind_Growth_plot <- ggplot(transform(Ind_Growth,
                                 copper_juju = factor(copper_juju, levels=c('Cu0 - Juju0', 'Cu5 - Juju0', 'Cu10 - Juju0', 'Cu25 - Juju0',
                                                                            'Cu0 - Juju0.1', 'Cu5 - Juju0.1', 'Cu10 - Juju0.1', 'Cu25 - Juju0.1',
                                                                            'Cu0 - Juju0.25', 'Cu5 - Juju0.25', 'Cu10 - Juju0.25', 'Cu25 - Juju0.25',
                                                                            'Cu0 - Juju0.5', 'Cu5 - Juju0.5', 'Cu10 - Juju0.5', 'Cu25 - Juju0.5'))),
                       aes(x = meanInd,y = meanGrowth, colour = as.factor(clone_type))) + 
  geom_point(size = 3, pch = 21) + 
  geom_errorbar(aes(ymin = Growth_min, ymax = Growth_max)) + 
  geom_errorbarh(aes(xmin = Ind_min, xmax = Ind_max)) + 
  facet_wrap(~copper_juju) + 
  theme(legend.position="top", , 
        rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        #strip.text.x = element_blank(),
        #axis.text.x = element_blank(), 
        #axis.title.x = element_blank(), 
        #axis.title.y = element_blank(),
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))

tiff(file = "Ind_Growth_plot.tiff", width = 1000, height = 800)
Ind_Growth_plot  + labs(x = "mean induction score", y = "mean growth @ maturity") + scale_x_continuous(limits=c(0,100), breaks=c(0,25,50,75,100)) + scale_color_manual(values=c("#FF0000","#000000"))
dev.off()













################
### OLD CODE ###
################


# run separate models for LOW and HIGH
### CHECK: in EXPAND DATA: clone_type = unique(SizeDat$clone_type)

## DATA IMPORT ---------------------------------------------------
AllData <- read.csv("data/response_surface_final_April2021.csv")
AllData <- as.data.table(AllData)
head(AllData)
names(AllData)
AllData[, cloneID := paste0(clone_type,"_",clone)]

# induction data
Ind <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
              ind_1_total = induction_1, 
              ind_2_total = induction_2, 
              ind_3_total = induction_3,
              ind_4_total = induction_4,
              
              ind_1_nteeth = induction_1_nteeth, 
              ind_2_nteeth = induction_2_nteeth, 
              ind_3_nteeth = induction_3_nteeth,
              ind_4_nteeth = induction_4_nteeth,
              
              ind_1_pedestal = induction_1_pedestal, 
              ind_2_pedestal = induction_2_pedestal, 
              ind_3_pedestal = induction_3_pedestal,
              ind_4_pedestal = induction_4_pedestal,
              
              # max induction across all days, including instar 1 
              maxInd_total = max_induction_total_I1,
              maxInd_nteeth = max_induction_nteeth_I1,
              maxInd_pedestal = max_induction_pedestal_I1,
              
              # max induction instars 2 and later 
              maxInd_total_I2 = max_induction_total_I2,
              maxInd_nteeth_I2 = max_induction_nteeth_I2,
              maxInd_pedestal_I2 = max_induction_pedestal_I2)

IndDat <- Ind[!maxInd_total=="NA"]

# induction across all instars
IndDat_use <- IndDat %>% 
  select(clone, clone_type, cloneID, juju, copper, maxInd_total) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

IndDat.ag <- na.omit(IndDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanInd = mean(maxInd_total), 
    seInd = sd(maxInd_total)/sqrt(sum(!is.na(maxInd_total))))

# induction in instars 2 and older
IndDat2_use <- IndDat %>% 
  select(clone, clone_type, cloneID, juju, copper, maxInd_total_I2) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

IndDat2.ag <- na.omit(IndDat2_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanInd = mean(maxInd_total_I2), 
    seInd = sd(maxInd_total_I2)/sqrt(sum(!is.na(maxInd_total_I2))))


# size data
Size <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
               size_mat = size_maturity)

SizeDat <- Size[!size_mat=="NA"]

SizeDat_use <- SizeDat %>% 
  select(clone, clone_type, cloneID, juju, copper, size_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

SizeDat.ag <- na.omit(SizeDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanSize = mean(size_mat), 
    seSize = sd(size_mat)/sqrt(sum(!is.na(size_mat))))


# growth data
Growth <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
                 growth_mat = somaticGrowth_maturity, 
                 growth_3d = somaticGrowth_3d, 
                 growth_4d = somaticGrowth_4d)

GrowthDat <- Growth[!growth_mat=="NA"]

GrowthDat_use <- GrowthDat %>% 
  select(clone, clone_type, cloneID, juju, copper, growth_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

GrowthDat.ag <- na.omit(GrowthDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanGrowth = mean(growth_mat), 
    seGrowth = sd(growth_mat)/sqrt(sum(!is.na(growth_mat))))


# age data
Age <- select(AllData, clone, clone_type, cloneID, juju, copper, rep, 
              age_mat = age_maturity)

AgeDat <- Age[!age_mat=="NA"]

AgeDat_use <- AgeDat %>% 
  select(clone, clone_type, cloneID, juju, copper, age_mat) %>% 
  mutate(copper2 = copper - mean(copper),
         juju2 = juju - mean(juju))

AgeDat.ag <- na.omit(AgeDat_use) %>% 
  group_by(copper, juju, cloneID) %>% 
  summarise(
    meanAge = mean(age_mat), 
    seAge = sd(age_mat)/sqrt(sum(!is.na(age_mat))))


## INDUCTION ---------------------------------------------------

# total induction (i.e., pedestal and neckteeth)

indSum_1 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_total),
    seInd_1 = sd(ind_1_total)/sqrt(sum(!is.na(ind_1_total))))

indSum_2 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_total),
    seInd_2 = sd(ind_2_total)/sqrt(sum(!is.na(ind_2_total))))

indSum_3 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_total),
    seInd_3 = sd(ind_3_total)/sqrt(sum(!is.na(ind_3_total))))

indSum_4 <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_total),
    seInd_4 = sd(ind_4_total)/sqrt(sum(!is.na(ind_4_total))))


plot_ind1 <- ggplot(transform(indSum_1, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_1, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,110) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind2 <- ggplot(transform(indSum_2, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_2, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,110) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind3 <- ggplot(transform(indSum_3, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_3, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,110) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind4 <- ggplot(transform(indSum_4, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                    aes(x = juju, y = meanInd_4, 
                        colour = factor(copper), group = factor(copper), 
                        ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,110) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 


# pedestal only

indSum_1_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_pedestal),
    seInd_1 = sd(ind_1_pedestal)/sqrt(sum(!is.na(ind_1_pedestal))))

indSum_2_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_pedestal),
    seInd_2 = sd(ind_2_pedestal)/sqrt(sum(!is.na(ind_2_pedestal))))

indSum_3_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_pedestal),
    seInd_3 = sd(ind_3_pedestal)/sqrt(sum(!is.na(ind_3_pedestal))))

indSum_4_ped <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_pedestal),
    seInd_4 = sd(ind_4_pedestal)/sqrt(sum(!is.na(ind_4_pedestal))))


plot_ind1_ped <- ggplot(transform(indSum_1_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_1, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind2_ped <- ggplot(transform(indSum_2_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_2, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind3_ped <- ggplot(transform(indSum_3_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_3, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind4_ped <- ggplot(transform(indSum_4_ped, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                        aes(x = juju, y = meanInd_4, 
                            colour = factor(copper), group = factor(copper), 
                            ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) + 
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_blank(),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 


# neckteeth only

indSum_1_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_1 = mean(ind_1_nteeth),
    seInd_1 = sd(ind_1_nteeth)/sqrt(sum(!is.na(ind_1_nteeth))))

indSum_2_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_2 = mean(ind_2_nteeth),
    seInd_2 = sd(ind_2_nteeth)/sqrt(sum(!is.na(ind_2_nteeth))))

indSum_3_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_3 = mean(ind_3_nteeth),
    seInd_3 = sd(ind_3_nteeth)/sqrt(sum(!is.na(ind_3_nteeth))))

indSum_4_nteeth <- na.omit(IndDat) %>% 
  group_by(clone, juju, copper) %>%
  summarise(
    meanInd_4 = mean(ind_4_nteeth),
    seInd_4 = sd(ind_4_nteeth)/sqrt(sum(!is.na(ind_4_nteeth))))

plot_ind1_nteeth <- ggplot(transform(indSum_1_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_1, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_1 - seInd_1, ymax = meanInd_1 + seInd_1))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind2_nteeth <- ggplot(transform(indSum_2_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_2, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_2 - seInd_2, ymax = meanInd_2 + seInd_2))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind3_nteeth <- ggplot(transform(indSum_3_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_3, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_3 - seInd_3, ymax = meanInd_3 + seInd_3))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 

plot_ind4_nteeth <- ggplot(transform(indSum_4_nteeth, clone = factor(clone, levels=c('C14', 'Chard', 'LD33', 'D86A', 'D87A', 'Cyril'))),
                           aes(x = juju, y = meanInd_4, 
                               colour = factor(copper), group = factor(copper), 
                               ymin = meanInd_4 - seInd_4, ymax = meanInd_4 + seInd_4))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  geom_line()+
  geom_errorbar(width = 0, size = 0.5)+
  facet_wrap(~clone)+
  scale_x_continuous(breaks=c(0,0.25,0.5)) +
  labs(y=expression(mean~induction[max]), x=expression(juju~(µl~ml^{-1}))) +
  ylim(0,60) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm")) 


## INDUCTION plot ---------------------------------------------------

patchwork_plots_induction <- (plot_ind1|plot_ind2|plot_ind3|plot_ind4) /
  (plot_ind1_ped|plot_ind2_ped|plot_ind3_ped|plot_ind4_ped) / 
  (plot_ind1_nteeth|plot_ind2_nteeth|plot_ind3_nteeth|plot_ind4_nteeth)

# tiff(file = "patchwork_plots_induction.tiff", width = 1000, height = 800)

patchwork_plots_induction + 
  plot_annotation(
    title = 'maximal induction @ days 1-4',
    subtitle = '(A-D) overall, (E-H) pedestal only, (I-L) neckteeth only', 
    tag_levels = 'A') + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

# dev.off()

## >> strongest induction in high inducers in instars 1 and 2 (slightly depending on clonal background)
## >> development of neckteeth and pedestal show congruent pattern, i.e. defence morph is not driven by one feature alone



## INDUCTION lm() ---------------------------------------------------
# maximal induction across all instars (including first instar, w/ induction score == pedestal + nteeth)
glimpse(IndDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(IndDat_use, aes(x = juju, y = maxInd_total))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
ind_mod_full_LOW <- lm(maxInd_total ~ poly(copper,2) + poly(juju,2) + copper:juju, data = IndDat_use[clone_type == "low"])
ind_mod_full_HIGH <- lm(maxInd_total ~ poly(copper,2) + poly(juju,2) + copper:juju, data = IndDat_use[clone_type == "high"])

summary(ind_mod_full_LOW)
# no evidence for 2nd order copper in coefficient
# no evidence for copper:juju either

summary(ind_mod_full_HIGH)
# no evidence for 2nd order copper in coefficient
# no evidence for copper:juju either


Confint(ind_mod_full_LOW)  
# confirmed by Confint from car package, there is no evidence to keep the interaction (see final model)

Confint(ind_mod_full_HIGH)  
# confirmed by Confint from car package, there is no evidence to keep the interaction (see final model)


# Type II ANOVA using car with Kenward-Roger df via F test.
Anova(ind_mod_full_LOW, type = "II", test = "F")
Anova(ind_mod_full_HIGH, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(ind_mod_full_LOW, type = "II")
Anova(ind_mod_full_HIGH, type = "II")
# Pretty much heading towards the conclusion that copper has no 2nd order effect or in interaction

# Reduce the copper polynomial
ind_mod_copper_linear_LOW <- lm(maxInd_total ~ copper + poly(juju,2) + copper:juju, data = IndDat_use[clone_type == "low"])
ind_mod_copper_linear_HIGH <- lm(maxInd_total ~ copper + poly(juju,2) + copper:juju, data = IndDat_use[clone_type == "high"])

anova(ind_mod_full_LOW, ind_mod_copper_linear_LOW)  # no difference, i.e. loss of term is justified by likelihood ratio test
anova(ind_mod_full_HIGH, ind_mod_copper_linear_HIGH)  # no difference, i.e. loss of term is justified by likelihood ratio test

summary(ind_mod_copper_linear_LOW) 
summary(ind_mod_copper_linear_HIGH) 

Confint(ind_mod_copper_linear_LOW)
Confint(ind_mod_copper_linear_HIGH)


# Include copper just as a single main effect (i.e., exclude interaction term)
ind_mod_copper_linear_noInt_LOW <- lm(maxInd_total ~ copper + poly(juju,2), data = IndDat_use[clone_type == "low"])
ind_mod_copper_linear_noInt_HIGH <- lm(maxInd_total ~ copper + poly(juju,2), data = IndDat_use[clone_type == "high"])

anova(ind_mod_copper_linear_LOW, ind_mod_copper_linear_noInt_LOW) # no difference, i.e. loss of term is justified by likelihood ratio test
anova(ind_mod_copper_linear_HIGH, ind_mod_copper_linear_noInt_HIGH) # no difference, i.e. loss of term is justified by likelihood ratio test

summary(ind_mod_copper_linear_noInt_LOW)
summary(ind_mod_copper_linear_noInt_HIGH)

Confint(ind_mod_copper_linear_noInt_LOW)
Confint(ind_mod_copper_linear_noInt_HIGH)

Anova(ind_mod_copper_linear_noInt_LOW)
Anova(ind_mod_copper_linear_noInt_HIGH)

# NO COPPER model 
ind_mod_no_copper_LOW <- lm(maxInd_total ~ poly(juju,2), data = IndDat_use[clone_type == "low"])
ind_mod_no_copper_HIGH <- lm(maxInd_total ~ poly(juju,2), data = IndDat_use[clone_type == "high"])

anova(ind_mod_copper_linear_noInt_LOW, ind_mod_no_copper_LOW)  # no difference, i.e. loss of term is justified by likelihood ratio test
anova(ind_mod_copper_linear_noInt_HIGH, ind_mod_no_copper_HIGH)  # linear copper term is better fit - keep! 

Confint(ind_mod_copper_linear_noInt_LOW)
Confint(ind_mod_copper_linear_noInt_HIGH)


# expand data
ind_newX_LOW <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(IndDat_use[clone_type == "low"])$clone_type)


ind_newX_HIGH <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(IndDat_use[clone_type == "high"])$clone_type)



# fixed effect is average among clones & clone effects are the clone-specific 'deviations' from average, their own effects
ind_fixed_effect_LOW_noInt <- predict(ind_mod_copper_linear_noInt_LOW, newdata = ind_newX_LOW, re.form = NA)
ind_fixed_effect_LOW_noCu <- predict(ind_mod_no_copper_LOW, newdata = ind_newX_LOW, re.form = NA)
ind_fixed_effect_HIGH_noInt <- predict(ind_mod_copper_linear_noInt_HIGH, newdata = ind_newX_HIGH, re.form = NA)

# housekeeping
ind_plotData_LOW_noInt <- data.frame(ind_newX_LOW, ind_fixed_effect = ind_fixed_effect_LOW_noInt)
ind_plotData_LOW_noCu <- data.frame(ind_newX_LOW, ind_fixed_effect = ind_fixed_effect_LOW_noCu)
ind_plotData_HIGH_noInt <- data.frame(ind_newX_HIGH, ind_fixed_effect = ind_fixed_effect_HIGH_noInt)

ind_plotData_noInt <- rbind(ind_plotData_LOW_noInt,ind_plotData_HIGH_noInt)
ind_plotData_noInt_noCu <- rbind(ind_plotData_LOW_noCu,ind_plotData_HIGH_noInt)


# plot the average - use noInt models (no huge difference between noInt and noCu model for LOW clones)
ind_Mod_average_noInt <- ggplot(ind_plotData_noInt, aes(x = juju, y = copper))+
  geom_raster(aes(fill = ind_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## INDUCTION lm() plot ---------------------------------------------------

ind_Mod_average_noInt <- ggplot(ind_plotData_noInt, aes(x = juju, y = copper))+
  geom_raster(aes(fill = ind_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## SIZE @ MATURITY lm() ---------------------------------------------------
glimpse(SizeDat_use)

# plot w/ polynomial order 2 fit in each panel
ggplot(SizeDat_use, aes(x = juju, y = size_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
size_mod_full_LOW <- lm(size_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = SizeDat_use[clone_type == "low"])
size_mod_full_HIGH <- lm(size_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = SizeDat_use[clone_type == "high"])

summary(size_mod_full_LOW)
summary(size_mod_full_HIGH)

Confint(size_mod_full_LOW) 
Confint(size_mod_full_HIGH) 
# no evidence for 2nd order copper in coefficient

# Type II ANOVA using car with KR df via F test.
Anova(size_mod_full_LOW, type = "II", test = "F")
Anova(size_mod_full_HIGH, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(size_mod_full_LOW, type = "II")
Anova(size_mod_full_HIGH, type = "II")
# Pretty much heading towards the conclusion that copper has no 2nd order effect

# Reduce the copper polynomial
size_mod_copper_linear_LOW <- lm(size_mat ~ copper + poly(juju,2) + copper:juju, data = SizeDat_use[clone_type == "low"])
size_mod_copper_linear_HIGH <- lm(size_mat ~ copper + poly(juju,2) + copper:juju, data = SizeDat_use[clone_type == "high"])

anova(size_mod_full_LOW, size_mod_copper_linear_LOW)  # no difference, loss of term is justified by likelihood ratio test
anova(size_mod_full_HIGH, size_mod_copper_linear_HIGH)  # no difference, loss of term is justified by likelihood ratio test

summary(size_mod_copper_linear_LOW) 
summary(size_mod_copper_linear_HIGH) 

Confint(size_mod_copper_linear_LOW)
Confint(size_mod_copper_linear_HIGH)


# Reduce the juju polynomial
size_mod_copper_juju_linear_LOW <- lm(size_mat ~ copper + juju + copper:juju, data = SizeDat_use[clone_type == "low"])
size_mod_copper_juju_linear_HIGH <- lm(size_mat ~ copper + juju + copper:juju, data = SizeDat_use[clone_type == "high"])

anova(size_mod_copper_linear_LOW, size_mod_copper_juju_linear_LOW)  # no difference, loss of term is justified by likelihood ratio test
anova(size_mod_copper_linear_HIGH, size_mod_copper_juju_linear_HIGH)  # no difference, loss of term is justified by likelihood ratio test

summary(size_mod_copper_juju_linear_LOW) 
summary(size_mod_copper_juju_linear_HIGH) 

Confint(size_mod_copper_juju_linear_LOW)
Confint(size_mod_copper_juju_linear_HIGH)


# NO JUJU model
size_mod_noJuju_LOW <- lm(size_mat ~ copper, data = SizeDat_use[clone_type == "low"])
size_mod_noJuju_HIGH <- lm(size_mat ~ copper, data = SizeDat_use[clone_type == "high"])

anova(size_mod_copper_juju_linear_LOW, size_mod_noJuju_LOW)  # linear juju term is better fit - keep! 
anova(size_mod_copper_juju_linear_HIGH, size_mod_noJuju_HIGH)  # linear juju term is better fit - keep! 

summary(size_mod_copper_juju_linear_LOW) 
summary(size_mod_copper_juju_linear_HIGH) 

Confint(size_mod_copper_juju_linear_LOW)
Confint(size_mod_copper_juju_linear_HIGH)


# NO interaction
size_mod_noInt_LOW <- lm(size_mat ~ copper + juju, data = SizeDat_use[clone_type == "low"])
size_mod_noInt_HIGH <- lm(size_mat ~ copper + juju, data = SizeDat_use[clone_type == "high"])

anova(size_mod_copper_juju_linear_LOW, size_mod_noInt_LOW)  # no difference, loss of term is justified by likelihood ratio test
anova(size_mod_copper_juju_linear_HIGH, size_mod_noInt_HIGH)  # no difference (but borderline 0.057), loss of term is justified by likelihood ratio test

summary(size_mod_noInt_LOW) 
summary(size_mod_noInt_HIGH) 

Confint(size_mod_noInt_LOW)
Confint(size_mod_noInt_HIGH)


## REPORTING
Anova(size_mod_noInt_LOW, type = "II", test = "F") # copper, juju and interaction are strongly affecting size @ maturity
Anova(size_mod_noInt_HIGH, type = "II", test = "F") # copper, juju and interaction are strongly affecting size @ maturity

summary(size_mod_noInt_LOW) 
summary(size_mod_noInt_HIGH) 

Confint(size_mod_noInt_LOW) 
Confint(size_mod_noInt_HIGH) 


# expand data
size_newX_LOW <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(SizeDat[clone_type == "low"])$clone_type)


size_newX_HIGH <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(SizeDat[clone_type == "high"])$clone_type)



# fixed effect is average among clones & clone effects are the clone-specific 'deviations' from average, their own effects
size_fixed_effect_LOW <- predict(size_mod_noInt_LOW, newdata = size_newX_LOW, re.form = NA)
size_fixed_effect_HIGH <- predict(size_mod_noInt_HIGH, newdata = size_newX_HIGH, re.form = NA)

# housekeeping
size_plotData_LOW <- data.frame(size_newX_LOW, size_fixed_effect=size_fixed_effect_LOW)
size_plotData_HIGH <- data.frame(size_newX_HIGH, size_fixed_effect=size_fixed_effect_HIGH)

size_plotData <- rbind(size_plotData_LOW,size_plotData_HIGH)


# plot the average and clone specifics
size_Mod_average <- ggplot(size_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = size_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## SIZE @ MATURITY lm() plot ---------------------------------------------------

size_Mod_average <- ggplot(size_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = size_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## GROWTHRATE MATURITY lm() ---------------------------------------------------
glimpse(GrowthDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(GrowthDat_use, aes(x = juju, y = growth_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
growth_mod_full_LOW <- lm(growth_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = GrowthDat_use[clone_type == "low"])
growth_mod_full_HIGH <- lm(growth_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = GrowthDat_use[clone_type == "high"])

summary(growth_mod_full_LOW)
Confint(growth_mod_full_LOW) 
# no evidence for 2nd order copper in coefficient
# no evidence for 2nd order juju in coefficient (if any at all)

summary(growth_mod_full_HIGH)
Confint(growth_mod_full_HIGH) 
# no evidence for 2nd order copper in coefficient
# no evidence for interaction term

# Type II ANOVA using car with KR df via F test.
Anova(growth_mod_full_LOW, type = "II", test = "F")
Anova(growth_mod_full_HIGH, type = "II", test = "F")

# Type II ANOVA using car with Wald test
Anova(growth_mod_full_LOW, type = "II")
Anova(growth_mod_full_HIGH, type = "II")

# Reduce the copper polynomial
growth_mod_linear_copper_LOW <- lm(growth_mat ~ copper + poly(juju,2) + copper:juju, data = GrowthDat_use[clone_type == "low"])
growth_mod_linear_copper_HIGH <- lm(growth_mat ~ copper + poly(juju,2) + copper:juju, data = GrowthDat_use[clone_type == "high"])

anova(growth_mod_full_LOW, growth_mod_linear_copper_LOW)  # no difference, i.e. loss of terms is justified by likelihood ratio test
anova(growth_mod_full_HIGH, growth_mod_linear_copper_HIGH)  # no difference, i.e. loss of terms is justified by likelihood ratio test

summary(growth_mod_linear_copper_LOW) 
summary(growth_mod_linear_copper_HIGH) 

Confint(growth_mod_linear_copper_LOW)
Confint(growth_mod_linear_copper_HIGH)

Anova(growth_mod_linear_copper_LOW, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in
Anova(growth_mod_linear_copper_HIGH, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in


# Reduce the juju polynomial (or entirely) for LOW and INTERACTION for HIGH
growth_mod_linear_copper_juju_LOW <- lm(growth_mat ~ copper + juju + copper:juju, data = GrowthDat_use[clone_type == "low"])
growth_mod_linear_copper_NoJuju_LOW <- lm(growth_mat ~ copper, data = GrowthDat_use[clone_type == "low"])
growth_mod_linear_copper_noInt_HIGH <- lm(growth_mat ~ copper + poly(juju,2), data = GrowthDat_use[clone_type == "high"])

anova(growth_mod_linear_copper_LOW, growth_mod_linear_copper_juju_LOW)  # no difference, i.e. loss of terms is justified by likelihood ratio test
anova(growth_mod_linear_copper_juju_LOW, growth_mod_linear_copper_NoJuju_LOW)  # keep juju in, better fit !
anova(growth_mod_linear_copper_HIGH, growth_mod_linear_copper_noInt_HIGH)  # no difference, i.e. loss of terms is justified by likelihood ratio test

summary(growth_mod_linear_copper_juju_LOW) 
summary(growth_mod_linear_copper_noInt_HIGH) 

Confint(growth_mod_linear_copper_juju_LOW)
Confint(growth_mod_linear_copper_noInt_HIGH)

Anova(growth_mod_linear_copper_juju_LOW, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in
Anova(growth_mod_linear_copper_noInt_HIGH, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in


# expand data
growth_newX_LOW <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(GrowthDat[clone_type == "low"])$clone_type)

growth_newX_HIGH <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(GrowthDat[clone_type == "high"])$clone_type)


# fixed effect is average among clones & clone effects are the clone specific 'deviations' from average, their own effects
growth_fixed_effect_LOW <- predict(growth_mod_linear_copper_juju_LOW, newdata = growth_newX_LOW, re.form = NA)
growth_fixed_effect_HIGH <- predict(growth_mod_linear_copper_noInt_HIGH, newdata = growth_newX_HIGH, re.form = NA)

# housekeeping
growth_plotData_LOW <- data.frame(growth_newX_LOW, growth_fixed_effect = growth_fixed_effect_LOW)
growth_plotData_HIGH <- data.frame(growth_newX_HIGH, growth_fixed_effect = growth_fixed_effect_HIGH)

growth_plotData <- rbind(growth_plotData_LOW, growth_plotData_HIGH)

# plot the average and clone specifics
growth_Mod_average <- ggplot(growth_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = growth_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))



## GROWTHRATE MATURITY lm() plot ---------------------------------------------------

growth_Mod_average <- ggplot(growth_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = growth_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))



## AGE @ MATURITY lm() ---------------------------------------------------
glimpse(AgeDat_use)

# The picture we are modelling, with polynomial order 2 fit in each panel
ggplot(AgeDat_use, aes(x = juju, y = age_mat))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se=FALSE)+
  facet_grid(clone ~ copper)

# The full RSM model
age_mod_full_LOW <- lm(age_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = AgeDat_use[clone_type == "low"])
age_mod_full_HIGH <- lm(age_mat ~ poly(copper,2) + poly(juju,2) + copper:juju, data = AgeDat_use[clone_type == "high"])

summary(age_mod_full_LOW)
Confint(age_mod_full_LOW) 
# no evidence for 2nd order copper and juju in coefficient

summary(age_mod_full_HIGH)
Confint(age_mod_full_HIGH) 
# no evidence for 2nd order copper in coefficient


# Reduce the copper polynomial (LOW) and interaction term (HIGH)
age_mod_linear_copper_LOW <- lm(age_mat ~ copper + poly(juju,2) + copper:juju, data = AgeDat_use[clone_type == "low"])
age_mod_linear_copper_juju_LOW <- lm(age_mat ~ copper + juju + copper:juju, data = AgeDat_use[clone_type == "low"])
age_mod_noInt_HIGH <- lm(age_mat ~ poly(copper,2) + poly(juju,2), data = AgeDat_use[clone_type == "high"])

anova(age_mod_full_LOW, age_mod_linear_copper_LOW) # no difference, i.e. loss of these terms is justified by likelihood ratio test
anova(age_mod_linear_copper_LOW, age_mod_linear_copper_juju_LOW) # no difference, i.e. loss of these terms is justified by likelihood ratio test
anova(age_mod_full_HIGH, age_mod_noInt_HIGH) # no difference, i.e. loss of these terms is justified by likelihood ratio test

summary(age_mod_linear_copper_juju_LOW) 
Confint(age_mod_linear_copper_juju_LOW)

summary(age_mod_noInt_HIGH) 
Confint(age_mod_noInt_HIGH)

Anova(age_mod_linear_copper_juju_LOW, type = "II", test = "F") 
Anova(age_mod_noInt_HIGH, type = "II", test = "F") 


# expand data
age_newX_LOW <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(AgeDat[clone_type == "low"])$clone_type)

age_newX_HIGH <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(AgeDat[clone_type == "high"])$clone_type)


# fixed effect is average among clones & clone effects are the clone specific 'deviations' from average, their own effects
age_fixed_effect_LOW <- predict(age_mod_linear_copper_juju_LOW, newdata = age_newX_LOW, re.form = NA)
age_fixed_effect_HIGH <- predict(age_mod_noInt_HIGH, newdata = age_newX_HIGH, re.form = NA)

# housekeeping
age_plotData_LOW <- data.frame(age_newX_LOW, age_fixed_effect = age_fixed_effect_LOW)
age_plotData_HIGH <- data.frame(age_newX_HIGH, age_fixed_effect = age_fixed_effect_HIGH)

age_plotData <- rbind(age_plotData_LOW,age_plotData_HIGH)


# plot the average and clone specifics
age_Mod_average <- ggplot(age_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = age_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## AGE @ MATURITY lm() plot ---------------------------------------------------

age_Mod_average <- ggplot(age_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = age_fixed_effect), interpolate = TRUE)+
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## SURVIVAL ANALYSIS ---------------------------------------------------

survival <- AllData[, c("clone","cloneID","clone_type","juju","copper","rep","survival_time","alive_dead")]
## survival data based on 'maturity data', i.e. data was recorded until animals were mature; 0 = alive, 1 = dead; 

# Cox Hazard Model (use 'coxph' when no random effects, 'coxme' when using random effects)
library(coxme)
mod1_LOW <- coxph(Surv(survival_time,alive_dead) ~ poly(juju,2) + poly(copper,2) + juju*copper, data = na.omit(survival[clone_type == "low"]))
mod1_HIGH <- coxph(Surv(survival_time,alive_dead) ~ poly(juju,2) + poly(copper,2) + juju*copper, data = na.omit(survival[clone_type == "high"]))

summary(mod1_LOW)
confint(mod1_LOW)

summary(mod1_HIGH)
confint(mod1_HIGH)

# Reduce polynomial for juju
mod2_LOW_linear_juju <- coxph(Surv(survival_time,alive_dead) ~ juju + poly(copper,2) + juju*copper, data = na.omit(survival[clone_type == "low"]))
mod2_HIGH_linear_juju <- coxph(Surv(survival_time,alive_dead) ~ juju + poly(copper,2) + juju*copper, data = na.omit(survival[clone_type == "high"]))

anova(mod1_LOW,mod2_LOW_linear_juju)   # no difference, drop of term justified
anova(mod1_HIGH,mod2_HIGH_linear_juju)   # no difference, drop of term justified

summary(mod2_LOW_noJuju)
confint(mod2_LOW_noJuju)

summary(mod2_HIGH_noJuju)
confint(mod2_HIGH_noJuju)


# Reduce polynomial for copper and interaction (LOW only) or juju entirely (LOW only)
mod3_LOW_linear_juju_copper <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[clone_type == "low"]))
mod4_LOW_linear_noInt <- coxph(Surv(survival_time,alive_dead) ~ juju + copper, data = na.omit(survival[clone_type == "low"]))
mod5_LOW_linear_copper_noJuju <- coxph(Surv(survival_time,alive_dead) ~ copper, data = na.omit(survival[clone_type == "low"]))
mod3_HIGH_linear_juju_copper <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival[clone_type == "high"]))

anova(mod2_LOW_linear_juju,mod3_LOW_linear_juju_copper)   # no difference, drop of term justified
anova(mod3_LOW_linear_juju_copper,mod4_LOW_linear_noInt)   # no difference, drop of term justified
anova(mod4_LOW_linear_noInt,mod5_LOW_linear_copper_noJuju)   # no difference, drop of term justified
anova(mod2_HIGH_linear_juju,mod3_HIGH_linear_juju_copper)   # no difference, drop of term justified

summary(mod5_LOW_linear_copper_noJuju)
confint(mod5_LOW_linear_copper_noJuju)

summary(mod3_HIGH_linear_juju_copper)
confint(mod3_HIGH_linear_juju_copper)

Anova(mod5_LOW_linear_copper_noJuju, type = "II")
Anova(mod3_HIGH_linear_juju_copper, type = "II")



# expand grid
surv_newX_LOW <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(survival[clone_type == "low"])$clone_type)

surv_newX_HIGH <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(survival[clone_type == "high"])$clone_type)


# fixed effects 
# surv_fixed_effects_lp_LOW <- predict(mod5_LOW_linear_copper_noJuju, type = "lp", newdata = surv_newX_LOW, re.form = NA)  # linear prediction
# surv_fixed_effects_lp_HIGH <- predict(mod3_HIGH_linear_juju_copper, type = "lp", newdata = surv_newX_HIGH, re.form = NA)  # linear prediction
surv_fixed_effects_risk_LOW <- predict(mod5_LOW_linear_copper_noJuju, type = "risk", newdata = surv_newX_LOW, re.form = NA)  # risk factor 
surv_fixed_effects_risk_HIGH <- predict(mod3_HIGH_linear_juju_copper, type = "risk", newdata = surv_newX_HIGH, re.form = NA)  # risk factor 

# housekeeping
surv_fixed_plotData_LOW <- data.table(surv_newX_LOW, surv_fixed_effects_risk = surv_fixed_effects_risk_LOW)
surv_fixed_plotData_HIGH <- data.table(surv_newX_HIGH, surv_fixed_effects_risk = surv_fixed_effects_risk_HIGH)

surv_fixed_plotData <- rbind(surv_fixed_plotData_LOW, surv_fixed_plotData_HIGH)


# plot the average and clone specifics
surv_Mod_average <- ggplot(surv_fixed_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = surv_fixed_effects_risk), interpolate = TRUE)+  # fill = surv_fixed_effects_risk
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))


## SURVIVAL ANALYSIS plot ---------------------------------------------------

surv_Mod_average <- ggplot(surv_fixed_plotData, aes(x = juju, y = copper))+
  geom_raster(aes(fill = surv_fixed_effects_risk), interpolate = TRUE)+  # fill = surv_fixed_effects_risk
  scale_fill_continuous(type = "viridis")+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  labs(y=expression(copper~(mg~L^{-1})), x=expression(juju~(µl~ml^{-1}))) +
  facet_wrap(~cloneID) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12, family='Arial'),
        axis.title.y = element_text(size=12, family='Arial'),
        axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12, family='Arial'),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))








