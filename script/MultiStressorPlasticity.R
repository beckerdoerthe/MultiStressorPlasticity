# surface design predation & copper / SHF-DFG project
# update DB 21 May 2021 

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

tiff(file = "patchwork_plots_induction.tiff", width = 1000, height = 800)

patchwork_plots_induction + 
  plot_annotation(
        title = 'maximal induction @ days 1-4',
        subtitle = '(A-D) overall, (E-H) pedestal only, (I-L) neckteeth only', 
        tag_levels = 'A') + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

dev.off()

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

Confint(ind_mod_copper_linear_noInt)

## for REPORTING, we go back to the full model and the model with linear copper
Anova(ind_mod_full, type = "II", test = "F")
Anova(ind_mod_copper_linear, type = "II", test = "F") # we see that copper linear
summary(ind_mod_copper_linear) # no evidence for interaction
Confint(ind_mod_copper_linear) # suggests even linear effect of copper might not be sig

# pbkrtests as alternative
pbkrtest::KRmodcomp(ind_mod_full, ind_mod_copper_linear) # is the polynomial copper significant 
pbkrtest::PBmodcomp(ind_mod_full, ind_mod_copper_linear) # parametric bootstrap

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
ind_plotData <- data.frame(ind_newX, ind_fixed_effect, ind_clone_effect)

# plot the average and clone specifics
ind_Mod_average <- ggplot(ind_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = ind_fixed_effect), interpolate = TRUE)+
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
                          geom_raster(aes(fill = ind_clone_effect), interpolate = TRUE)+
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
                juju = seq(0,0.5,length = 100),
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
                copper = seq(0,25,length = 100),
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

tiff(file = "patchwork_plots_induction_model.tiff", width = 800, height = 800)

patchwork_plots_induction_model <- ind_average / (ind_Mod_average | ind_Mod_byClone) / (ind_juju_clone | ind_copper_clone | guide_area() )

patchwork_plots_induction_model + 
  plot_annotation(
    title = 'maximal induction @ days 1-4',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

dev.off()


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

## REPORTING
Anova(size_mod_copper_linear, type = "II", test = "F") # copper, juju and interaction are strongly affecting size @ maturity
summary(size_mod_copper_linear) 
Confint(size_mod_copper_linear) 

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
size_plotData <- data.frame(size_newX, size_fixed_effect, size_clone_effect)

# plot the average and clone specifics
size_Mod_average <- ggplot(size_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = size_fixed_effect), interpolate = TRUE)+
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
                          geom_raster(aes(fill = size_clone_effect), interpolate = TRUE)+
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
  juju = seq(0,0.5,length = 100),
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
  copper = seq(0,25,length = 100),
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

tiff(file = "patchwork_plots_size_model.tiff", width = 800, height = 800)

patchwork_plots_size_model <- size_average / (size_Mod_average | size_Mod_byClone) / (size_juju_clone | size_copper_clone | guide_area() )

patchwork_plots_size_model + 
  plot_annotation(
    title = 'size @ maturity',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

dev.off()


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

summary(growth_mod_linear) 
Confint(growth_mod_linear)

Anova(growth_mod_linear, type = "II", test = "F") # Juju has no effect on growthrate, but there is a strong interaction effect > leave juju in

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
growth_plotData <- data.frame(growth_newX, growth_fixed_effect, growth_clone_effect)

# plot the average and clone specifics
growth_Mod_average <- ggplot(growth_plotData, aes(x = juju, y = copper))+
                          geom_raster(aes(fill = growth_fixed_effect), interpolate = TRUE)+
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
                            geom_raster(aes(fill = growth_clone_effect), interpolate = TRUE)+
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
  juju = seq(0,0.5,length = 100),
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
  copper = seq(0,25,length = 100),
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

tiff(file = "patchwork_plots_growth_model.tiff", width = 800, height = 800)

patchwork_plots_growth_model <- growth_average / (growth_Mod_average | growth_Mod_byClone) / (growth_juju_clone | growth_copper_clone | guide_area() )

patchwork_plots_growth_model + 
  plot_annotation(
    title = 'growth rate',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

dev.off()


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
age_plotData <- data.frame(age_newX, age_fixed_effect, age_clone_effect)

# plot the average and clone specifics
age_Mod_average <- ggplot(age_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = age_fixed_effect), interpolate = TRUE)+
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
                      geom_raster(aes(fill = age_clone_effect), interpolate = TRUE)+
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
  juju = seq(0,0.5,length = 100),
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
  copper = seq(0,25,length = 100),
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

tiff(file = "patchwork_plots_age_model.tiff", width = 800, height = 800)

patchwork_plots_age_model <- age_average / (age_Mod_average | age_Mod_byClone) / (age_juju_clone | age_copper_clone | guide_area() )

patchwork_plots_age_model + 
  plot_annotation(
    title = 'age @ maturity',
    subtitle = '(A) effect of Cu and Juju on genetically distinct clones (B) average effect of Cu and JuJu among clones, (C) clone specific effects of Cu and JuJu, 
    (D) effect of JuJu among clones, (E) effect of Cu among clones', 
    tag_levels = 'A') + 
  plot_layout(nrow=3, height = c(3,6,3), guides = 'collect') & theme(legend.position = "right")

dev.off()


## SURVIVAL ANALYSIS ---------------------------------------------------

survival <- AllData[, c("clone","cloneID","juju","copper","rep","survival_time","alive_dead")]
## survival data based on 'maturity data', i.e. data was recorded until animals were mature; 0 = alive, 1 = dead; 

# get distance 
dd <- datadist(survival)
options(datadist="dd")

## PSM model: parametric survival model; assumes weibull distribution; takes a two column response variable Surv(time, code)

# juju, copper, and interaction
mod1 <- psm(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival))
summary(mod1)
anova(mod1)

# clone as a fixed effect 
mod2 <- psm(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + clone, data = na.omit(survival))  
summary(mod2)
anova(mod2)

# mod2 <- update(mod2, x=TRUE,y=TRUE)       
# validate(mod2, B=300)  # increase to 300 or so
# plot(Predict(mod2, juju, copper, clone)) 

# response surface
# low responders
par(mfrow = c(3,4))
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "D86A")
title("Cu = 0, D86A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "D86A")
title("Cu = 5, D86A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "D86A")
title("Cu = 10, D86A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "D86A")
title("Cu = 25, D86A")

survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "D87A")
title("Cu = 0, D87A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "D87A")
title("Cu = 5, D87A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "D87A")
title("Cu = 10, D87A")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "D87A")
title("Cu = 25, D87A")

survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "Cyril")
title("Cu = 0, Cyril")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "Cyril")
title("Cu = 5, Cyril")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "Cyril")
title("Cu = 10, Cyril")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "Cyril")
title("Cu = 25, Cyril")


# high responders
par(mfrow = c(3,4))
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "C14")
title("Cu = 0, C14")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "C14")
title("Cu = 5, C14")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "C14")
title("Cu = 10, C14")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "C14")
title("Cu = 25, C14")

survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "Chard")
title("Cu = 0, Chard")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "Chard")
title("Cu = 5, Chard")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "Chard")
title("Cu = 10, Chard")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "Chard")
title("Cu = 25, Chard")

survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 0, clone = "LD33")
title("Cu = 0, LD33")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 5, clone = "LD33")
title("Cu = 5, LD33")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 10, clone = "LD33")
title("Cu = 10, LD33")
survplot(mod2, juju=c(0,0.1,0.25,0.5), copper = 25, clone = "LD33")
title("Cu = 25, LD33")

## for copper 0, 5, and 10, INCREASE in juju DECREASES survival probability.
## for copper 25, INCREASE in juju INCREASES survival probability.
## this is true for all clones, even though we do see clone-specific survival rates across treatments.


# Cox Hazard Model ('coxph': no random effects, 'coxme' when using random effects)
library(coxme)
# clone as fixed effect
mod3 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + clone, data = na.omit(survival))
summary(mod3)

# clone as random effect
mod4 <- coxme(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + (1|clone), data = na.omit(survival))
summary(mod4)

anova(mod3,mod4)  # random effect model (mod4) is better


# exclude clone entirely ?
mod5 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper, data = na.omit(survival))
summary(mod5)

anova(mod4,mod5)  # excluding clone (mod5) fits better

# model with robust SE via clustering
mod6 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + cluster(clone), data = na.omit(survival))
summary(mod6)

anova(mod5,mod6)  # same, use mod 5

## model with a frailty term 
mod7 <- coxph(Surv(survival_time,alive_dead) ~ juju + copper + juju*copper + frailty(clone), data = na.omit(survival))
summary(mod7)
anova(mod6,mod7)  # frailty term model fits better

Anova(mod7, type = "II")

## juju effect is n.s., while copper has a strong effect on survival; there is also a small juju:copper interaction effect. 


## @Andrew - what is the frailty(clone) term ACTUALLY telling us. Not sure I do get the concept of fraily entirely.... 
## I thought the frailty(clone) is pretty much a random effect on our clones (i.e., + (1|clone) ), but why is mod7 then any 
## better than mod4 ???
anova(mod7,mod4)


# ftest <- cox.zph(mod7)  # test Cox fit; NOT working with frailty term model...
ggforest(mod7, data=survival) ## copper has a strong effect on survival, juju effect is n.s.

# these look a bit messy.... 
ggcoxdiagnostics(mod7, type = "deviance", ox.scale = "linear.predictions")  
ggcoxdiagnostics(mod7, type = "deviance", ox.scale = "observation.id")  


# expand grid
surv_newX <- expand.grid(
  juju = seq(0,0.5,length = 10),
  copper = seq(0,25, length = 10),  
  cloneID = unique(survival$cloneID)
)

surv_effects_lp_avg <- predict(mod7, type = "lp", newdata = surv_newX)  # linear prediction
surv_effects_risk <- predict(mod7, type = "risk", newdata = surv_newX)  # risk factor 

## @Andrew - I can't figure out how to get the clone-specific predictions from the model. I.e., something similar to 
## what we used before: predict(model, type=XX, newdata=YY, re.form = ~(1|clone)) when using an lmer() w/ random effect.
## And, given that the frailty(clone) term is significant, I guess we should plot the resonse surface for each clone separately. 

surv_plotData <- data.frame(surv_newX, surv_effects_lp, risk=surv_effects_risk)

## plot --------------------------------------------------------
surv_Mod_average <- ggplot(surv_plotData, aes(x = juju, y = copper))+
                        geom_raster(aes(fill = risk), interpolate = TRUE)+
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


