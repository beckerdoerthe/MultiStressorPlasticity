### FIG 7

# packages
library(data.table)
library(tidyverse)
library(patchwork)
library(viridis)
library(cowplot)

load(file = "Fig7_data.RData")

## Induction vs other life history traits - plots --------------------------

# plot extremes (i.e., [copper == 0][juju == 0], [copper == 0][juju == 0.5], [copper == 25][juju == 0], [copper == 25][juju == 0.5])

# Induction vs Age
Ind_age <- ggplot(transform(Ind_Age_data, copper_juju = factor(copper_juju, levels=c("Cu0 - Juju0", "Cu0 - Juju0.5", "Cu25 - Juju0.5", "Cu25 - Juju0"))),
                  aes(y = meanInd, x = meanAge, colour = clone_type)) + 
                geom_line(aes(group = clone_type), linetype="twodash", colour="darkgrey") + 
                geom_point(colour = rep(c("#FF0000", "#000000"),4), size = 3) + 
                geom_errorbarh(aes(xmin = Age_min, xmax = Age_max), height = 2, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
                geom_errorbar(aes(ymin = Ind_min, ymax = Ind_max), width = 0.01, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
                labs(x=expression(age[maturity]~(days)), y=expression(mean~induction[max])) +
                xlim(4, 13) +
                ylim(0, 100) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=12),
                      axis.title.y = element_text(size=12),
                      # axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                      axis.text = element_text(size=12),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm"))


# Induction vs Size
Ind_size <- ggplot(transform(Ind_Size_data, copper_juju = factor(copper_juju, levels=c("Cu0 - Juju0", "Cu0 - Juju0.5", "Cu25 - Juju0.5", "Cu25 - Juju0"))),
                   aes(y = meanInd, x = meanSize, colour = clone_type)) + 
                geom_line(aes(group = clone_type), linetype="twodash", colour="darkgrey") + 
                geom_point(colour = rep(c("#FF0000", "#000000"),4), size = 3) + 
                geom_errorbarh(aes(xmin = Size_min, xmax = Size_max), height = 2, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
                geom_errorbar(aes(ymin = Ind_min, ymax = Ind_max), width = 0.01, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
                labs(x=expression(size[maturity]~(mm)), y=expression(mean~induction[max])) +
                xlim(1.4, 2.1) +
                ylim(0, 100) +
                theme(rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=12),
                      axis.title.y = element_text(size=12),
                      # axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
                      axis.text = element_text(size=12),
                      strip.text.x = element_text(size =10, color = "black"),
                      strip.text.y = element_text(size =10, color = "black"),
                      panel.spacing.x = unit(4, "mm"),
                      panel.spacing.y = unit(1, "mm"))


# Induction vs GrowthRate
Ind_Growth_data <- Ind_Growth[copper_juju %in% c("Cu0 - Juju0", "Cu25 - Juju0", "Cu0 - Juju0.5", "Cu25 - Juju0.5")]

Ind_growth <- ggplot(transform(Ind_Growth_data, copper_juju = factor(copper_juju, levels=c("Cu0 - Juju0", "Cu0 - Juju0.5", "Cu25 - Juju0.5", "Cu25 - Juju0"))),
                     aes(y = meanInd, x = meanGrowth, colour = clone_type)) + 
  geom_line(aes(group = clone_type), linetype="twodash", colour="darkgrey") + 
  geom_point(colour = rep(c("#FF0000", "#000000"),4), size = 3) + 
  geom_errorbarh(aes(xmin = Growth_min, xmax = Growth_max), height = 2, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
  geom_errorbar(aes(ymin = Ind_min, ymax = Ind_max), width = 0.01, colour = rep(c("#FF0000", "#000000"),4), size = 0.5) + 
  labs(x=expression(growth~rate~(mm~d^{-1})), y=expression(mean~induction[max])) +
  xlim(0.05, 0.21) +
  ylim(0, 100) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(colour = "grey70", size=0.25),
        panel.grid.minor = element_line(colour = "grey90", size=0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(size = 1),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        # axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size =10, color = "black"),
        strip.text.y = element_text(size =10, color = "black"),
        panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(1, "mm"))




## all plots together
patchwork_plots_HvsL <- Ind_size | Ind_growth | Ind_age

setEPS()
postscript("Fig7.eps", colormodel = "RGB", width = 8, height = 5)
patchwork_plots_HvsL +  
  plot_annotation(
    tag_levels = 'A') +
  plot_layout(ncol=3, width = c(1,1,1), guides = 'collect') & theme(legend.position = "bottom")
dev.off()


## stats for figure 7 ---------------------

## Do LOW and HIGH responders differ?

# data to use:
Size_data_fig7 <- na.omit(SizeDat_use[copper == 0 & juju == 0 | copper == 25 & juju == 0 | copper == 0 & juju == 0.5 | copper == 25 & juju == 0.5])
Growth_data_fig7 <- na.omit(GrowthDat_NoLn_use[copper == 0 & juju == 0 | copper == 25 & juju == 0 | copper == 0 & juju == 0.5 | copper == 25 & juju == 0.5])
Age_data_fig7 <- na.omit(AgeDat_use[copper == 0 & juju == 0 | copper == 25 & juju == 0 | copper == 0 & juju == 0.5 | copper == 25 & juju == 0.5])

## SIZE @ MATURITY
# data normally distributed? 
shapiro.test(Size_data_fig7$size_mat) # NO(!) > use wilcox.test()

# t.test(size_mat ~ clone_type, data = Size_data_fig7)
wilcox_SizeMat <- wilcox.test(size_mat ~ clone_type, data = Size_data_fig7[copper == 0 & juju == 0.5]) 
wilcox_SizeMat

Zscore = qnorm(wilcox_SizeMat$p.value/2)
Zscore

effectSize <- wilcox_effsize(size_mat ~ clone_type, data = Size_data_fig7[copper == 0 & juju == 0.5])
effectSize

# [copper == 0 & juju == 0]
# (Z = -2.402412, effect size = 0.467, p-value = 0.01629) **

# [copper == 25 & juju == 0]
# (Z = -0.08365173, effect size = 0.0439, p-value = 0.9333)

# [copper == 25 & juju == 0.5]
# (Z = -2.486406, effect size = 0.312, p-value = 0.2325)

# [copper == 0 & juju == 0.5]
# (Z = -0.3717383, effect size = 0.0948, p-value = 0.7101)


## SOMATIC GROWTH RATE
# data normally distributed? 
shapiro.test(Growth_data_fig7$growth_mat_NoLn) # NO(!) > use wilcox.test()

# t.test(growth_mat_NoLn ~ clone_type, data = Growth_data_fig7)
wilcox_GrowthMat <- wilcox.test(growth_mat_NoLn ~ clone_type, data = Growth_data_fig7[copper == 0 & juju == 0.5]) 
wilcox_GrowthMat

Zscore = qnorm(wilcox_GrowthMat$p.value/2)
Zscore

effectSize <- wilcox_effsize(growth_mat_NoLn ~ clone_type, data = Growth_data_fig7[copper == 0 & juju == 0.5])
effectSize

# [copper == 0 & juju == 0]
# (Z = -2.770555, effect size = 0.523, p-value = 0.005596) ***

# [copper == 25 & juju == 0]
# (Z = -1.57922, effect size = 0.439, p-value = 0.1143)

# [copper == 25 & juju == 0.5]
# (Z = -3.481024, effect size = 0.786, p-value = 0.0004995) ***

# [copper == 0 & juju == 0.5]
# (Z = -0.2014007, effect size = 0.0568, p-value = 0.8404)


## AGE @ MATURITY
# data normally distributed? 
shapiro.test(Age_data_fig7$age_mat) # NO(!) > use wilcox.test()

# t.test(age_mat ~ clone_type, data = Age_data_fig7)
wilcox_AgeMat <- wilcox.test(age_mat ~ clone_type, data = Age_data_fig7[copper == 0 & juju == 0.5]) 
wilcox_AgeMat

Zscore = qnorm(wilcox_AgeMat$p.value/2)
Zscore

effectSize <- wilcox_effsize(age_mat ~ clone_type, data = Age_data_fig7[copper == 0 & juju == 0.5])
effectSize

# [copper == 0 & juju == 0]
# (Z = -0.1743115, effect size = 0.0391, p-value = 0.8616) 

# [copper == 25 & juju == 0]
# (Z = -1.653857, effect size = 0.449, p-value = 0.09816)

# [copper == 25 & juju == 0.5]
# (Z = -3.069801, effect size = 0.781, p-value = 0.002142) ***

# [copper == 0 & juju == 0.5]
# (Z = -0.6389078, effect size = 0.160, p-value = 0.5229)


## Do different treatments induce a shift in low and high responders ?

Size_data_fig7[, trt_group := ifelse(juju == 0 & copper == 0, "control", 
                                     ifelse(juju == 0.5 & copper == 0, "highPred_noCu", 
                                            ifelse(juju == 0 & copper == 25, "highCu_noPred", 
                                                   ifelse(juju == 0.5 & copper == 25, "highCu_highPred", NA))))]

Growth_data_fig7[, trt_group := ifelse(juju == 0 & copper == 0, "control", 
                                       ifelse(juju == 0.5 & copper == 0, "highPred_noCu", 
                                              ifelse(juju == 0 & copper == 25, "highCu_noPred", 
                                                     ifelse(juju == 0.5 & copper == 25, "highCu_highPred", NA))))]


Age_data_fig7[, trt_group := ifelse(juju == 0 & copper == 0, "control", 
                                    ifelse(juju == 0.5 & copper == 0, "highPred_noCu", 
                                           ifelse(juju == 0 & copper == 25, "highCu_noPred", 
                                                  ifelse(juju == 0.5 & copper == 25, "highCu_highPred", NA))))]


library(car)
## SIZE @ MATURITY
# LOWs
SizeLOW_levene <- leveneTest(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "low"]) 
SizeLOW_levene # significant difference between variances
SizeLOW_oneway <- oneway.test(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "low"])
SizeLOW_oneway # *** differences between treatment groups

# ANOVA assumptions not ok, use Kruskal Wallis test
kruskal.test(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "low"])
pairwise.wilcox.test(Size_data_fig7[clone_type == "low"]$size_mat,
                     Size_data_fig7[clone_type == "low"]$trt_group,
                     p.adjust.method = "BH")  # Copper exposure (with no added predator cue) is different from other three treatments


# HIGHs
SizeHIGH_levene <- leveneTest(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "high"]) 
SizeHIGH_levene # not significant difference between variances (but use same test for all groups, i.e. Kruskal Wallis)
SizeHIGH_oneway <- oneway.test(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "high"])
SizeHIGH_oneway # *** differences between treatment groups

SizeHIGH_aov <- aov(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "high"])
summary(SizeHIGH_aov)
TukeyHSD(SizeHIGH_aov)


# ANOVA assumptions not ok, use Kruskal Wallis test
kruskal.test(size_mat ~ trt_group, data = Size_data_fig7[clone_type == "high"])
pairwise.wilcox.test(Size_data_fig7[clone_type == "high"]$size_mat,
                     Size_data_fig7[clone_type == "high"]$trt_group,
                     p.adjust.method = "BH")  # Copper exposure (with no added predator cue) is different from other three treatments


