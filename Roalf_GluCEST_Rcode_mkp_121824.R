#load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)
library(cowplot)
library(corrplot)
library(RColorBrewer)
library(ggsignif)
library(kableExtra)
library(gtsummary)
library(table1)
library(lmerTest)
library(visreg)
library(plyr)
library(sjPlot)
source("/Users/pecsok/Library/CloudStorage/Box-Box/GluCEST\ PhD/Manuscripts/Aging\ GluCEST/NbioAging/summary_se_function.R")
#source("https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/summary")
#####################################################################################
#read in data CSV
root="/Users/pecsok/Library/CloudStorage/Box-Box/GluCEST\ PhD/Manuscripts/Aging\ GluCEST/NbioAging"
#data cleaned 6/20/2023 see ~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/final_analysis/data/Aging_GluCEST_UNI_noMCI_cleaned_raw.xlsx for details
#the few MCI subjects were removed
glucest_age<-read.csv(file.path(root,"Aging_GluCEST_wide.csv"))
#read in cest long format:
glucest_age_long<-read.csv(file.path(root,"age_glucest_long_cestonly.csv"))
#read in motion data:
glucest_age_mvt<-read.csv(file.path(root,"aging_motion_fs.csv"))
#read in 3d hippo volume data:
#HOAhippo_3dvol <-read.csv(file.path(root,"Resubmission/all_aseg_stats_102524.csv")) 
#HYAhippo_3dvol <-read.csv(file.path(root,"Resubmission/young_bblids_aseg_stats_csvformat.csv"))
hippo_3dvol <-read.csv(file.path(root,"Resubmission/aseg_stats_compiled_110424.csv"))

cci_longdf <-read.csv(file.path(root,"age_glucest_long_cestonly_cci.csv"))

#####################################################################################
#add motion data
for (i in 1:nrow(glucest_age)) {
  match_row <- glucest_age_mvt[
    glucest_age_mvt$bblid == glucest_age$bblid[i] & 
      glucest_age_mvt$scanid == glucest_age$sessionid[i], 
  ]
    if (nrow(match_row) > 0) {
      print(match_row)
      glucest_age$motion_max[i] <- match_row$motion_max
      glucest_age$motion_mean[i] <- match_row$motion_mean
  }
}

#add 3d hippo data
#first, add bblid and scan id as aseg df only has inddids
hippo_3dvol$bblid = NaN
hippo_3dvol$scanid = NaN
for (i in 1:nrow(glucest_age_mvt)) {
  #add bblid and scanid based on inddid
  for (j in 1:nrow(hippo_3dvol)) {
    mvt_id = glucest_age_mvt$inddid[i]
    bbl_id = glucest_age_mvt$bblid[i]
    vol_id = hippo_3dvol$Measure.volume[j]
    if (!is.na(mvt_id) && !is.na(vol_id) && vol_id == mvt_id) {
      hippo_3dvol$bblid[j] = glucest_age_mvt$bblid[i]
      hippo_3dvol$scanid[j] = glucest_age_mvt$scanid[i]
    }
  }
}

#add 3d hippo data, aligning by bblid and scanid for HOA
glucest_age$bblid <- as.character(glucest_age$bblid)
hippo_3dvol$bblid <- as.character(hippo_3dvol$bblid)
glucest_age$sessionid <- as.character(glucest_age$sessionid)
hippo_3dvol$scanid <- as.character(hippo_3dvol$scanid)
glucest_age$LeftHippVol_3D<-NaN
glucest_age$RightHippVol_3D<-NaN
#
for (i in 1:nrow(glucest_age)) {
  #HOA
  match_row <- hippo_3dvol[hippo_3dvol$bblid == glucest_age$bblid[i] &
      hippo_3dvol$scanid == glucest_age$sessionid[i],]
  if (nrow(match_row) == 1) {
    glucest_age$LeftHippVol_3D[i] <- match_row$Left.Hippocampus
    glucest_age$RightHippVol_3D[i] <- match_row$Right.Hippocampus
    glucest_age$eTIV[i] <- match_row$EstimatedTotalIntraCranialVol
  }
  #HYA
  match_row <- hippo_3dvol[hippo_3dvol$Measure.volume == glucest_age$bblid[i],]
  if (nrow(match_row) == 1) {
    glucest_age$LeftHippVol_3D[i] <- match_row$Left.Hippocampus
    glucest_age$RightHippVol_3D[i] <- match_row$Right.Hippocampus
    glucest_age$eTIV[i] <- match_row$EstimatedTotalIntraCranialVol
  }
}
glucest_age$bblid<-as.character(glucest_age$bblid)
glucest_age_long$bblid<-as.character(glucest_age_long$bblid)


#Add 3D volume and Age to the long dataframe
merged_df <- glucest_age_long %>%
  left_join(glucest_age %>% select(bblid, LeftHippVol_3D, RightHippVol_3D, age, cci_total), by = "bblid") %>%
  mutate(Vol_3D = ifelse(Hemi == "Left", LeftHippVol_3D,
                         ifelse(Hemi == "Right", RightHippVol_3D, NA)))
glucest_age_long <- cbind(glucest_age_long, merged_df[c("Vol_3D", "age","cci_total")])


#####################################################################################
#make group and gender factors
glucest_age$group<-as.factor(glucest_age$group)
glucest_age$gender<-as.factor(glucest_age$gender)

#Demographics
glucest.agemean <-summarySE(data=glucest_age,measurevar="age",groupvars=c("group"),na.rm=T)
glucest.agemean.sex <-summarySE(data=glucest_age,measurevar="age",groupvars=c("group", "gender"),na.rm=T)
glucest.moca <-summarySE(data=glucest_age,measurevar="mocatotal",groupvars=c("group"),na.rm=T)
glucest.moca.sex <-summarySE(data=glucest_age,measurevar="mocatotal",groupvars=c("group", "gender"),na.rm=T)
glucest.cci <-summarySE(data=glucest_age,measurevar="cci_total",groupvars=c("group"),na.rm=T)

#simple t.tests for data exploration and comparison
t.test(mocatotal~group, data=glucest_age)
t.test(cci_total~group, data=glucest_age)
t.test(RightHippoGluCEST~group, data=glucest_age)
t.test(LeftHippoGluCEST~group, data=glucest_age)
t.test(motion_max~group, data=glucest_age)
t.test(motion_mean~group, data=glucest_age)
t.test(LeftHippVol_3D~group, data=glucest_age)
t.test(RightHippVol_3D~group, data=glucest_age)
t.test(eTIV~group, data=glucest_age)


t.test(Vol_3D~Hemi, data=glucest_age_long[glucest_age_long$group=="HOA",])
t.test(Vol_3D~Hemi, data=glucest_age_long[glucest_age_long$group=="HYA",])
t.test(CEST~Hemi, data=glucest_age_long[glucest_age_long$group=="HOA",])
t.test(CEST~Hemi, data=glucest_age_long[glucest_age_long$group=="HYA",])

# Gender Sensitivity analysis for Revision
t.test(CEST~Hemi, data=glucest_age_long[glucest_age_long$gender=="Male",])
t.test(CEST~Hemi, data=glucest_age_long[glucest_age_long$gender=="Female",])

#Limit race variable to White, Boc, and Other given small/no sample size in some races
glucest_age$race2<-glucest_age$race
glucest_age$race2[which(glucest_age$race2=="White")]<- "White"
glucest_age$race2[which(glucest_age$race2=="Black or African American")]<- "Black or African American"
glucest_age$race2[which(glucest_age$race2=="Hawaiian/Pacific Islander")]<- "Other"
glucest_age$race2[which(glucest_age$race2=="More than one race")]<- "Other"
glucest_age$race2[which(glucest_age$race2=="Asian")]<- "Other"

#Simple correlations b/w MoCA vs cest
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$mocatotal)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$mocatotal)
t.test(mocatotal~gender,glucest_age)
t.test(cci_total~gender,glucest_age[glucest_age$group == "HOA",])
t.test(LeftHippoGluCEST~gender,glucest_age)
t.test(RightHippoGluCEST~gender,glucest_age)

#Simple correlations b/w CCI vs cest
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$cci_total)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$cci_total)

#Simple correlations b/w 3D volume vs cest
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$LeftHippVol_3D)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$RightHippVol_3D)
glucest_old<-glucest_age[glucest_age$group=="HOA",] 
cor.test(glucest_old$LeftHippoGluCEST,glucest_old$LeftHippVol_3D)
cor.test(glucest_old$RightHippoGluCEST,glucest_old$RightHippVol_3D)

#GluCEST and Volume descriptives by age group and sex
glucest.age.lhipp <-summarySE(data=glucest_age,measurevar="LeftHippoGluCEST",groupvars=c("group"),na.rm=T)
glucest.age.sex.lhipp <-summarySE(data=glucest_age,measurevar="LeftHippoGluCEST",groupvars=c("group", "gender"),na.rm=T)
glucest.age.rhipp <-summarySE(data=glucest_age,measurevar="RightHippoGluCEST",groupvars=c("group"),na.rm=T)
glucest.age.sex.rhipp <-summarySE(data=glucest_age,measurevar="RightHippoGluCEST",groupvars=c("group", "gender"),na.rm=T)
glucest.asymm <-summarySE(data=glucest_age,measurevar="LminusRasymmGluCEST",groupvars=c("group"),na.rm=T)
glucest.sex.asymm <-summarySE(data=glucest_age,measurevar="LminusRasymmGluCEST",groupvars=c("group", "gender"),na.rm=T)
glucest.vol.lhipp<-summarySE(data=glucest_age,measurevar="LeftHippGluCESTVol",groupvars=c("group"),na.rm=T)
glucest.vol.sex.lhipp <-summarySE(data=glucest_age,measurevar="LeftHippGluCESTVol",groupvars=c("group", "gender"),na.rm=T)
glucest.vol.rhipp<-summarySE(data=glucest_age,measurevar="RightHippGluCESTVol",groupvars=c("group"),na.rm=T)
glucest.vol.sex.rhipp <-summarySE(data=glucest_age,measurevar="RightHippGluCESTVol",groupvars=c("group", "gender"),na.rm=T)
vol.age.asymm <-summarySE(data=glucest_age,measurevar="LminusRasymmVol",groupvars=c("group"),na.rm=T)
vol.age.sex.asymm <-summarySE(data=glucest_age,measurevar="LminusRasymmVol",groupvars=c("group", "gender"),na.rm=T)


glucest.maxmvt<-summarySE(data=glucest_age,measurevar="motion_max",groupvars=c("group"),na.rm=T)
glucest.sex.maxmvt <-summarySE(data=glucest_age,measurevar="motion_max",groupvars=c("group", "gender"),na.rm=T)
glucest.meanmvt<-summarySE(data=glucest_age,measurevar="motion_mean",groupvars=c("group"),na.rm=T)
glucest.sex.meanmvt <-summarySE(data=glucest_age,measurevar="motion_mean",groupvars=c("group", "gender"),na.rm=T)


#Table 1 demographic comparisons
## function for pvalues (categorical and continuous)
pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

#reformat cols for table1
glucest_age$Sex<-factor(glucest_age$gender)
glucest_age$Age<-glucest_age$age
glucest_age$Race<-factor(glucest_age$race2)
glucest_age$MoCA<-glucest_age$mocatotal
glucest_age$CCI<-glucest_age$cci_total
#make table1
table1(~ Age + Sex+ Race +MoCA+CCI | group, data=glucest_age, overall=F, extra.col=list(`P-value`=pvalue))


#####################################################################################
# Table 2: Imaging Variables

glucest_age$eTIV_sci <- formatC(glucest_age$eTIV, format = "e", digits = 3)

# Fix
tableS1<-table1(~ motion_mean + motion_max + LeftHippVol_3D + RightHippVol_3D + eTIV | group, data=glucest_age, overall=F, extra.col=list(`P-value`=pvalue))
tableS1_df<-as.data.frame(tableS1)
write.csv(tableS1_df,file.path(root,"Resubmission/Supplemental_TableS1.csv"))


#####################################################################################
# Table 3: Biological sex effects 

glucest_age$eTIV_sci <- formatC(glucest_age$eTIV, format = "e", digits = 3)

# Fix
tableS3<-table1(~ cci_total + mocatotal + LeftHippoGluCEST + RightHippoGluCEST | gender, data=glucest_age, overall=F, extra.col=list(`P-value`=pvalue))
tableS3_df<-as.data.frame(tableS3)

tableS3<-table1(~ cci_total + mocatotal + LeftHippoGluCEST + RightHippoGluCEST | gender, data=glucest_age[glucest_age$group=="HOA",], overall=F, extra.col=list(`P-value`=pvalue))


#write.csv(tableS1_df,file.path(root,"Resubmission/Supplemental_TableSx.csv"))





#####################################################################################
# FIGURE 2
#make gender and group factors in long data set
glucest_age_long$gender<-factor(glucest_age_long$gender)
glucest_age_long$group<-factor(glucest_age_long$group)

#summary stats for CEST by group and hemisphere
glucest.by.age <-summarySE(data=glucest_age_long,measurevar="CEST",groupvars=c("group"),na.rm=T)
glucest.by.hemi <-summarySE(data=glucest_age_long,measurevar="CEST",groupvars=c("Hemi"),na.rm=T)


##lme model
#library(lme4)
lmer2.cest<- lmer(CEST~ group * Hemi + gender+Vol+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)

visreg(lmer2.cest,"group",by="Hemi")
#post-hoc comparisons; use in plot function below. 
lmean.age.hemi<-lsmeans(lmer2.cest,pairwise~group * Hemi)
d<-summary(lsmeans(lmer2.cest, ~group * Hemi))

left.ls.data<-d[ which(d$Hemi=='Left'),]
right.ls.data<-d[ which(d$Hemi=='Right'),]

#Figure2A
left.model.lsmeans.plot <- ggplot(data = left.ls.data, aes(x = group, y = lsmean)) +
  geom_bar(aes(fill = group), stat = "identity", show.legend = FALSE, width = 0.75) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.8)) +
  geom_point(
    data = glucest_age,
    aes(x = group, y = LeftHippoGluCEST, color = group),  # Use group to map colors for dots
    position = position_jitter(width = 0.2),
    size = 3.5,
    alpha = 0.75,
    stroke = 0
  ) +
  scale_color_manual(values = c("#0d4f3b", "#823901")) +  # Darker shades of the bar colors
  ggtitle("Left Hippocampus") +
  ylim(0, 10.5) +
  xlab("Group") +
  ylab("GluCEST % Contrast") +
  scale_x_discrete(labels = c("HOA", "HYA")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank()
  )
ggsave(file.path(root,"Resubmission/Figure2A.pdf"), width=4, height=10)
#Figure2B
right.model.lsmeans.plot <- ggplot(data = right.ls.data, aes(x = group, y = lsmean)) +
  geom_bar(aes(fill = group), stat = "identity", show.legend = FALSE, width = 0.75) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.8)) +
  geom_point(
    data = glucest_age,
    aes(x = group, y = RightHippoGluCEST, color = group),  # Use group to map colors for dots
    position = position_jitter(width = 0.2),
    size = 3.5,
    alpha = 0.75,
    stroke = 0
  ) +
  scale_color_manual(values = c("#0d4f3b", "#823901")) +  # Darker shades of the bar colors
  ggtitle("Right Hippocampus") +
  ylim(0, 10.5) +
  xlab("Group") +
  ylab("GluCEST % Contrast") +
  scale_x_discrete(labels = c("HOA", "HYA")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank()
  )
ggsave(file.path(root,"Resubmission/Figure2B.pdf"), width=4, height=10)

#####################################################################################
#Figure2C Asymmetry Hippocampus GluCEST
#quick t test for asymmetry
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age)
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age[glucest_age$gender=="Female",])
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age[glucest_age$gender=="Male",])
                                                        
                                                        
#Convert to percent difference
glucest_age$AsymmCESTPercent <- ((glucest_age$LeftHippoGluCEST - glucest_age$RightHippoGluCEST) / ((glucest_age$LeftHippoGluCEST + glucest_age$RightHippoGluCEST) / 2))*100
glucest_age$AsymmVolPercent <- ((glucest_age$LeftHippGluCESTVol - glucest_age$RightHippGluCESTVol) / ((glucest_age$LeftHippGluCESTVol + glucest_age$RightHippGluCESTVol) / 2))*100
glucest_age$AsymmVol3DPercent <- ((glucest_age$LeftHippVol_3D - glucest_age$RightHippVol_3D) / ((glucest_age$LeftHippVol_3D + glucest_age$RightHippVol_3D) / 2))*100

#lm on asymmetry
#lm.asym<- lm(LminusRasymmGluCEST~ group + gender+LminusRasymmVol, data=glucest_age)
lm.asym<- lm(AsymmCESTPercent~ group + gender+AsymmVolPercent, data=glucest_age)
anova(lm.asym)
pval_group <- anova(lm.asym)["group", "Pr(>F)"]

lsmean.asym<-lsmeans(lm.asym,pairwise~group)
e<-summary(lsmeans(lm.asym, ~group))

Figure3<-ggplot(e, aes(x = group, y = lsmean, fill = group)) +
  geom_bar(aes(fill = group), stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), color = "black", width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(-10,30) +
  labs(
    title = "Left-Right Asymmetry",
    x = "Group",
    y = "% Asymmetry"
  ) +
  geom_point(
    data = glucest_age,
    aes(x = group, y = AsymmCESTPercent, color = group),  # Use group to map colors for dots
    position = position_jitter(width = 0.2),
    size = 3.5,
    alpha = 0.6,
    stroke = 0
  ) +
  scale_color_manual(values = c("#0d4f3b", "#823901")) +  # Darker shades of the bar colors
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())+
  geom_signif(
    annotation=paste("p=", formatC(pval_group, format = "f", digits = )),
    #annotation=formatC(anno.22q.hemi,digits =1), #Bug
    y_position=28.5, xmin =1, xmax = 2,
    tip_length = c(0.03, 0.15))
Figure3
ggsave(file.path(root,"Resubmission/Figure2C.pdf"), Figure3,
       width=4, height=10)


#####################################################################################
# FIGURE 3
library(sjPlot)
theme_set(theme_sjplot())

#CCI by GluCEST and group
cci_longdf$group<-as.factor(cci_longdf$group)
cci_longdf$Hemi<-as.factor(cci_longdf$Hemi)

cci_longdf <- cci_longdf %>% filter(!is.na(cci_total))

graph_df<-cci_longdf[cci_longdf$Hemi=="Left",]
lmer2.cci <- lm(cci_total ~ group * CEST, data=graph_df)




# FINAL MODEL 
lmer2.cci <- lm(cci_total ~ group * CEST*Hemi, data=cci_longdf)

anova(lmer2.cci)
fig = plot_model(lmer2.cci, type = "int")
fig[[1]]

# FIGURE 3 plot
model1b <- lm(cci_total ~ LeftHippoGluCEST *group + RightHippoGluCEST * group , data = glucest_age)
summary(model1b)
fig = plot_model(model1b, type = "int")
plotL = fig[[1]]
plotL <- plotL +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    title = "",
    x = "Left Hippocampus GluCEST",
    y = "Predicted CCI Total") +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/Figure3A_cci_interaction_plotL.png"), plot = plotL, height=6,width=6)

#CCI by GluCEST and group
model1a <- lm(cci_total ~ RightHippoGluCEST *group  + LeftHippoGluCEST * group , data = glucest_age)
summary(model1a)
fig = plot_model(model1a, type = "int")
plotR = fig[[1]]
plotR <- plotR +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    title = "",
    x = "Right Hippocampus GluCEST",
    y = "Predicted CCI Total") +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
#ggsave(file.path(root,"Resubmission/Figure3B_cci_interaction_plotR.png"), plot = plotR, height=6,width=6)


#model1b <- lm(cci_total ~ CEST *group, data = glucest_age_long[])
#summary(model1b)
#fig = plot_model(model1b, type = "int")
#plot1 = fig[[1]]
#ggsave(file.path(root,"Resubmission/Analysis/cci_interaction_plotL.png"), plot = plot1)





###############################################################################
#Figure 3: Behavioral Interaction
# Fit the model
model1a <- lm(cci_total ~ RightHippoGluCEST * group + LeftHippoGluCEST * group, data = glucest_age)

model1a <- lm(cci_total ~ RightHippoGluCEST * group + LeftHippoGluCEST * group, data = glucest_age)



# Plot the interaction with 95% confidence intervals using ggplot2
left_int <- ggplot(pred_data, aes(x = LeftHippoGluCEST, y = predicted_cci_total, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    #title = "Interaction Effect of Right and Left Hippocampus GluCEST on CCI Total",
    x = "Left Hippocampus GluCEST",
    y = "Predicted CCI Total"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
  
# Plot the interaction with 95% confidence intervals using ggplot2
right_int <- ggplot(pred_data, aes(x = RightHippoGluCEST, y = predicted_cci_total, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    #title = "Interaction Effect of Right and Left Hippocampus GluCEST on CCI Total",
    x = "Right Hippocampus GluCEST",
    y = "Predicted CCI Total"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())+
  
  theme(plot.title = element_text(hjust = 0.5))










#MoCA by GluCEST and group
model2a <- lm(mocatotal ~ RightHippoGluCEST *group  + LeftHippoGluCEST * group , data = glucest_age)
summary(model2a)
fig = plot_model(model2a, type = "int")
ggsave(file.path(root,"Resubmission/Analysis/moca_interaction_plotL.png"))


###############################################################################
#Figure 3 Sensivity: Just women


lmer2.cci_women <- lm(cci_total ~ group * CEST*Hemi, data=cci_longdf[cci_longdf$gender=="Female",])
anova(lmer2.cci_women)
lmer2.cci <- lm(cci_total ~ group * CEST*Hemi, data=cci_longdf)
anova(lmer2.cci)

# OLD model
model1b <- lm(cci_total ~ LeftHippoGluCEST *group + RightHippoGluCEST * group , data = glucest_age)
summary(model1b)
fig = plot_model(model1b, type = "int")
plotL = fig[[1]]
plotL <- plotL +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    title = "",
    x = "Left Hippocampus GluCEST",
    y = "Predicted CCI Total") +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/Figure3A_cci_interaction_plotL_WOMEN.png"), plot = plotL, height=6,width=6)

#CCI by GluCEST and group
model1a <- lm(cci_total ~ RightHippoGluCEST *group  + LeftHippoGluCEST * group , data = glucest_age)
summary(model1a)
fig = plot_model(model1a, type = "int")
plotR = fig[[1]]
plotR <- plotR +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  labs(
    title = "",
    x = "Right Hippocampus GluCEST",
    y = "Predicted CCI Total") +
  theme_minimal(base_size = 15) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "gray90"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
#ggsave(file.path(root,"Resubmission/Figure3B_cci_interaction_plotR.png"), plot = plotR, height=6,width=6)







#####################################################################################
#Sensitivity: Motion data
ggplot(data = glucest_age, aes(x = group, y = motion_mean, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5) +  # Box plot without outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Add jittered points for each subject
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +  # Custom colors
  labs(x = "Age Group", y = "Motion Mean", title = "Comparison of Motion Mean by Age Group") +
  theme_minimal() +  # Minimal theme for cleaner look
  theme(legend.position = "none")  # Hide legend since it's redundant
ggsave(file.path(root,"Resubmission/Analysis/age_by_motion_mean_boxplot.pdf"))


ggplot(data = glucest_age, aes(x = group, y = motion_max, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5) +  # Box plot without outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Add jittered points for each subject
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +  # Custom colors
  labs(x = "Age Group", y = "Motion Max", title = "Comparison of Motion Max by Age Group") +
  theme_minimal() +  # Minimal theme for cleaner look
  theme(legend.position = "none")  # Hide legend since it's redundant
ggsave(file.path(root,"Resubmission/Analysis/age_by_motion_max_boxplot.pdf"))

#####################################################################################
#Sensitivity: Volume

#Simple correlations b/w 2D volume vs cest
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$LeftHippGluCESTVol)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$RightHippGluCESTVol)

#Simple correlation b/w 3D volume and age
cor.test(glucest_age$LeftHippVol_3D,glucest_age$age)
cor.test(glucest_age$RightHippVol_3D,glucest_age$age)

#Simple correlations b/w 3D volume vs cest
young<-glucest_age[glucest_age$group=="HYA",]
old<-glucest_age[glucest_age$group=="HOA",]

cor.test(glucest_age$LeftHippoGluCEST,glucest_age$LeftHippVol_3D)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$RightHippVol_3D)
cor.test(young$LeftHippoGluCEST,young$LeftHippVol_3D)
cor.test(young$RightHippoGluCEST,young$RightHippVol_3D)
cor.test(old$LeftHippoGluCEST,old$LeftHippVol_3D)
cor.test(old$RightHippoGluCEST,old$RightHippVol_3D)

summary(lm(CEST ~ group + Vol_3D*Hemi, data = merged_df))

cor.test(glucest_age$RightHippoGluCEST,glucest_age$motion_max)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$motion_mean)
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$motion_max)
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$motion_mean)

glucest_old = glucest_age[glucest_age$group=="HOA",]
cor.test(glucest_old$RightHippoGluCEST,glucest_old$motion_max)
cor.test(glucest_old$RightHippoGluCEST,glucest_old$motion_mean)


glucest_young = glucest_age[glucest_age$group=="HYA",]
cor.test(glucest_young$RightHippoGluCEST,glucest_young$motion_max)
cor.test(glucest_young$RightHippoGluCEST,glucest_young$motion_mean)


############################################
# Figure 2 Analysis with 3D Volume: 

# Sensitivity analysis: Sensitivity_3DVol_CEST_Age
interaction_model <- lm(CEST ~ group + Hemi * Vol_3D, data = glucest_age_long)
summary(interaction_model)

# Since there's so much colinearity between age and 3d volume:
# 1a. Take residuals of 3d_vol vs group
vol_model <- lm(Vol_3D ~ Age, data = glucest_age_long, na.action = na.exclude)
summary(vol_model)
glucest_age_long$Vol_3D_resid <- residuals(vol_model)
# 1b. Center the Vol_3D variable.
glucest_age_long$centered_Vol_3D = glucest_age_long$Vol_3D - mean(glucest_age_long$Vol_3D, na.rm = TRUE)

# Now let's compare different variations of the model. Which one makes the most sense? 
#Original model:
lmer2.cest<- lmer(CEST~ group * Hemi + gender+ Vol+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)
#Main analysis with 3D volume added: Sensitivity_Figure2_anova_3dvol
lmer2.cest2 <- lmer(CEST~ group * Hemi + gender+ Vol_3D+(1|bblid), data=glucest_age_long)
anova(lmer2.cest2)
#Sensitivity_3DVol_Age_Hemi
interaction_model <- lm(Vol_3D ~ group * Hemi , data = glucest_age_long)
summary(interaction_model)
############################################
# Figure 2C Analysis with 3D Volume: 

# Create LminusR vol variable.
glucest_age['LminusRasymm_3DVol']<- glucest_age$LeftHippVol_3D - glucest_age$RightHippVol_3D 

# Now let's compare different variations of the model. Which one makes the most sense? 
#Original model:
lm.asym<- lm(LminusRasymmGluCEST~ group + gender+LminusRasymmVol, data=glucest_age)
anova(lm.asym)
#Asymmetry model with 3D volume asymmetry
lm.asym<- lm(LminusRasymmGluCEST~ group + gender + LminusRasymm_3DVol, data=glucest_age)
anova(lm.asym)
#Asymmetry model of 3D volume
lm.asym<- lm(LminusRasymm_3DVol~ group + gender, data=glucest_age)
anova(lm.asym)

#CCI by 3dVOLUME and group
model1b <- lm(cci_total ~ LeftHippVol_3D* group + RightHippVol_3D*group , data = glucest_age)
summary(model1b)
fig = plot_model(model1b, type = "int")
plot1 = fig[[1]]
ggsave(file.path(root,"Resubmission/Analysis/cci_interaction_plotL.png"), plot = plot1)



#####################################################################################
# FIGURE 2 SENSIVITY 

#lme model
#library(lme4)
lmer2.cest<- lmer(CEST~ group * Hemi + gender+Vol_3D+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)

visreg(lmer2.cest,"group",by="Hemi")
#post-hoc comparisons; use in plot function below. 
lmean.age.hemi<-lsmeans(lmer2.cest,pairwise~group * Hemi)
d<-summary(lsmeans(lmer2.cest, ~group * Hemi))

left.ls.data<-d[ which(d$Hemi=='Left'),]
right.ls.data<-d[ which(d$Hemi=='Right'),]

#Figure2A
left.model.lsmeans.plot<-ggplot(data=left.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ggtitle("Left Hippocampus")+
  ylim(0,10)+
  xlab("Group")+
  ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA")) + 
  theme_minimal(base_size = 15) +  # Use a minimal theme as a base
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/FigureS2A.pdf"), width=4, height=10)
#Figure2B
right.model.lsmeans.plot<-ggplot(data=right.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ggtitle("Right Hippocampus")+
  ylim(0,10)+
  xlab("Group")+
  ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA"))+
  theme_minimal(base_size = 15) +  # Use a minimal theme as a base
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/FigureS2B.pdf"), width=4, height=10)


#####################################################################################
#Figure2C Asymmetry Hippocampus GluCEST
#quick t test for asymmetry
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age)

#Convert to percent difference
glucest_age$AsymmCESTPercent <- ((glucest_age$LeftHippoGluCEST - glucest_age$RightHippoGluCEST) / ((glucest_age$LeftHippoGluCEST + glucest_age$RightHippoGluCEST) / 2))*100
glucest_age$AsymmVolPercent <- ((glucest_age$LeftHippGluCESTVol - glucest_age$RightHippGluCESTVol) / ((glucest_age$LeftHippGluCESTVol + glucest_age$RightHippGluCESTVol) / 2))*100
glucest_age$AsymmVol3DPercent <- ((glucest_age$LeftHippVol_3D - glucest_age$RightHippVol_3D) / ((glucest_age$LeftHippVol_3D + glucest_age$RightHippVol_3D) / 2))*100

#lm on asymmetry
#lm.asym<- lm(LminusRasymmGluCEST~ group + gender+LminusRasymmVol, data=glucest_age)
lm.asym<- lm(AsymmCESTPercent~ group + gender+AsymmVol3DPercent, data=glucest_age)
anova(lm.asym)
pval_group <- anova(lm.asym)["group", "Pr(>F)"]

lsmean.asym<-lsmeans(lm.asym,pairwise~group)
e<-summary(lsmeans(lm.asym, ~group))

Figure3<-ggplot(e, aes(x = group, y = lsmean, fill = group)) +
  geom_bar(aes(color = group), stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), color = "black", width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(-10,30) +
  labs(
    title = "Left-Right Asymmetry",
    x = "Group",
    y = "% Asymmetry"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())+
  geom_signif(
    annotation=paste("p=", formatC(pval_group, format = "f", digits = )),
    #annotation=formatC(anno.22q.hemi,digits =1), #Bug
    y_position=28.5, xmin =1, xmax = 2,
    tip_length = c(0.1, 0.20))

ggsave(file.path(root,"Resubmission/FigureS2C.pdf"), Figure3,
       width=4, height=10)



#####################################################################################
# FIGURE 2 SENSIVITY MEN VS WOMEN

#lme model
#library(lme4)
glucest_age_long_men <- glucest_age_long[glucest_age_long$gender=="Male",]
glucest_age_long_women <- glucest_age_long[glucest_age_long$gender=="Female",]


lmer2.cest<- lmer(CEST~ group * Hemi + gender+Vol_3D+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)

visreg(lmer2.cest,"group",by="Hemi")
#post-hoc comparisons; use in plot function below. 
lmean.age.hemi<-lsmeans(lmer2.cest,pairwise~group * Hemi)
d<-summary(lsmeans(lmer2.cest, ~group * Hemi))

left.ls.data<-d[ which(d$Hemi=='Left'),]
right.ls.data<-d[ which(d$Hemi=='Right'),]

#Figure2A
left.model.lsmeans.plot<-ggplot(data=left.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ggtitle("Left Hippocampus")+
  ylim(0,10)+
  xlab("Group")+
  ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA")) + 
  theme_minimal(base_size = 15) +  # Use a minimal theme as a base
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/FigureS2A.pdf"), width=4, height=10)
#Figure2B
right.model.lsmeans.plot<-ggplot(data=right.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ggtitle("Right Hippocampus")+
  ylim(0,10)+
  xlab("Group")+
  ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA"))+
  theme_minimal(base_size = 15) +  # Use a minimal theme as a base
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())
ggsave(file.path(root,"Resubmission/FigureS2B.pdf"), width=4, height=10)


#####################################################################################
#Figure2C Asymmetry Hippocampus GluCEST
#quick t test for asymmetry
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age)

#Convert to percent difference
glucest_age$AsymmCESTPercent <- ((glucest_age$LeftHippoGluCEST - glucest_age$RightHippoGluCEST) / ((glucest_age$LeftHippoGluCEST + glucest_age$RightHippoGluCEST) / 2))*100
glucest_age$AsymmVolPercent <- ((glucest_age$LeftHippGluCESTVol - glucest_age$RightHippGluCESTVol) / ((glucest_age$LeftHippGluCESTVol + glucest_age$RightHippGluCESTVol) / 2))*100
glucest_age$AsymmVol3DPercent <- ((glucest_age$LeftHippVol_3D - glucest_age$RightHippVol_3D) / ((glucest_age$LeftHippVol_3D + glucest_age$RightHippVol_3D) / 2))*100

#lm on asymmetry
#lm.asym<- lm(LminusRasymmGluCEST~ group + gender+LminusRasymmVol, data=glucest_age)
lm.asym<- lm(AsymmCESTPercent~ group + gender+AsymmVol3DPercent, data=glucest_age)
anova(lm.asym)
pval_group <- anova(lm.asym)["group", "Pr(>F)"]

lsmean.asym<-lsmeans(lm.asym,pairwise~group)
e<-summary(lsmeans(lm.asym, ~group))

Figure3<-ggplot(e, aes(x = group, y = lsmean, fill = group)) +
  geom_bar(aes(color = group), stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), color = "black", width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(-10,30) +
  labs(
    title = "Left-Right Asymmetry",
    x = "Group",
    y = "% Asymmetry"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "gray90", color = NA),  # Set gray background
    panel.grid.major = element_line(color = "white"),             # Light gray for major grid lines
    panel.grid.minor = element_blank())+
  geom_signif(
    annotation=paste("p=", formatC(pval_group, format = "f", digits = )),
    #annotation=formatC(anno.22q.hemi,digits =1), #Bug
    y_position=28.5, xmin =1, xmax = 2,
    tip_length = c(0.1, 0.20))

ggsave(file.path(root,"Resubmission/FigureS2C.pdf"), Figure3,
       width=4, height=10)







# Sensitivity analysis: Sensitivity_Figure3_asymm_3dvol
lm.asym<- lm(AsymmCESTPercent~ group + gender + AsymmVol3DPercent, data=glucest_age)
anova(lm.asym)
pval_group <- anova(lm.asym)["group", "Pr(>F)"]

# Sensitivity analysis: Sensitivity_Figure3_asymm_3dvol
lm.asym<- lm(AsymmVol3DPercent~ group + gender, data=glucest_age)
anova(lm.asym)
pval_group <- anova(lm.asym)["group", "Pr(>F)"]

# Sensitivity

lmer2.cest<- lmer(CEST~ group * Hemi + group*gender+Vol+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)
lmer2.cest<- lmer(CEST~ group * Hemi + gender+Vol+(1|bblid), data=glucest_age_long)


############ 
# Power analysis
library(pwr)
# Parameters
n1 <- 27  # Sample size for group 1 (old)
n2 <- 22  # Sample size for group 2 (young)
alpha <- 0.05  # Significance level
power <- 0.8  # Desired power
effect_size = 0.8 


# Calculate Cohen's f
sum_sq_effect <- 10.5877  # Sum of squares for group:Hemi interaction
sum_sq_error <- 46.757    # Residual sum of squares (error)
cohens_f <- sqrt(sum_sq_effect / (sum_sq_effect + sum_sq_error))

# Calculate the required sample size for each group

pow <- pwr.anova.test(k = 4, f = cohens_f, n = 24.5, sig.level = 0.006)


# Perform the power analysis
result <- pwr.t2n.test(n1 = n1, n2 = n2, sig.level = alpha, power = power)
result$d

# Output the required sample size for each group
ceiling(result$n)










