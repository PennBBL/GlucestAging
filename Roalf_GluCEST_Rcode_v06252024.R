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
source("~/Dropbox/R/summary_se_function.R")
#source("https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/summary")

#read in data CSV
#data cleaned 6/20/2023 see ~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/final_analysis/data/Aging_GluCEST_UNI_noMCI_cleaned_raw.xlsx for details
#the few MCI subjects were removed
glucest_age<-read.csv("~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/Manuscript/NbioAging/Aging_GluCEST_wide.csv")
#read in cest long format;
glucest_age_long<-read.csv("~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/Manuscript/NbioAging/age_glucest_long_cestonly.csv")

#make group and gender factors
glucest_age$group<-as.factor(glucest_age$group)
glucest_age$gender<-as.factor(glucest_age$gender)


#Demographics
glucest.agemean <-summarySE(data=glucest_age,measurevar="age",groupvars=c("group"),na.rm=T)
glucest.agemean.sex <-summarySE(data=glucest_age,measurevar="age",groupvars=c("group", "gender"),na.rm=T)
glucest.moca <-summarySE(data=glucest_age,measurevar="mocatotal",groupvars=c("group"),na.rm=T)
glucest.moca.sex <-summarySE(data=glucest_age,measurevar="mocatotal",groupvars=c("group", "gender"),na.rm=T)
glucest.cci <-summarySE(data=glucest_age,measurevar="cci_total",groupvars=c("group"),na.rm=T)

#Other scale not currently used
#glucest.tics <-summarySE(data=glucest_age,measurevar="ticstotal",groupvars=c("group"),na.rm=T)
#glucest.bdi <-summarySE(data=glucest_age,measurevar="bdi_total",groupvars=c("group"),na.rm=T)
#glucest.faq <-summarySE(data=glucest_age,measurevar="faq_total",groupvars=c("group"),na.rm=T)
#glucest.gds <-summarySE(data=glucest_age,measurevar="gds_total",groupvars=c("group"),na.rm=T)

#simple t.tests for data exploration and comparison
t.test(mocatotal~group, data=glucest_age)
t.test(cci_total~group, data=glucest_age)
t.test(RightHippoGluCEST~group, data=glucest_age)
t.test(LeftHippoGluCEST~group, data=glucest_age)

#not used
#t.test(faq_total~group, data=glucest_age)
#t.test(gds_total~group, data=glucest_age)

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

#Simple correlations b/w CCI vs cest
cor.test(glucest_age$LeftHippoGluCEST,glucest_age$cci_total)
cor.test(glucest_age$RightHippoGluCEST,glucest_age$cci_total)

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

#Table 1 demographic comparisons
library(table1)
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

glucest_age$Sex<-factor(glucest_age$gender)

glucest_age$Age<-glucest_age$age
glucest_age$Race<-factor(glucest_age$race2)

glucest_age$MoCA<-glucest_age$mocatotal
glucest_age$CCI<-glucest_age$cci_total
table1(~ Age + Sex+ Race +MoCA+CCI | group, data=glucest_age, overall=F, extra.col=list(`P-value`=pvalue))


#make gender and group factors in long data set
glucest_age_long$gender<-factor(glucest_age_long$gender)
glucest_age_long$group<-factor(glucest_age_long$group)

#summary stats for CEST by group and hemisphere
glucest.by.age <-summarySE(data=glucest_age_long,measurevar="CEST",groupvars=c("group"),na.rm=T)
glucest.by.hemi <-summarySE(data=glucest_age_long,measurevar="CEST",groupvars=c("Hemi"),na.rm=T)



##lme model
#library(lme4)
library(lmerTest)
lmer2.cest<- lmer(CEST~ group * Hemi + gender+Vol+(1|bblid), data=glucest_age_long)
anova(lmer2.cest)
library(visreg)
visreg(lmer2.cest,"group",by="Hemi")
#post-hoc comparisons; use in plot function below. 
lmean.age.hemi<-lsmeans(lmer2.cest,pairwise~group * Hemi)
d<-summary(lsmeans(lmer2.cest, ~group * Hemi))

left.ls.data<-d[ which(d$Hemi=='Left'),]
right.ls.data<-d[ which(d$Hemi=='Right'),]

left.model.lsmeans.plot<-ggplot(data=left.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ylim(0,10)+
  xlab("Groups")+
  ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA"))
  

right.model.lsmeans.plot<-ggplot(data=right.ls.data, aes(x=group, y=lsmean))+
  geom_bar( aes(fill = group),stat="identity", show.legend = FALSE, width = 0.75)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,position=position_dodge(.9))+
  ylim(0,10)+
  xlab("Groups")+
  #ylab("GluCEST % Constrast")+
  scale_x_discrete(labels=c("HOA","HYA"))+
  theme(axis.title.y = element_blank())

age.hemi.plot<-plot_grid(left.model.lsmeans.plot, right.model.lsmeans.plot)
ggsave("~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/final_analysis/data/age_by_hemi_GluCEST_boxplot.pdf")


#Asymmetry Hippocampus GluCEST
#quick t test for asymmetry
glucest.asymm.test <-t.test(LminusRasymmGluCEST~group, data=glucest_age)

#lm on asymmetry
lm.asym<- lm(LminusRasymmGluCEST~ group + gender+LminusRasymmVol, data=glucest_age)
anova(lm.asym)

lsmean.asym<-lsmeans(lm.asym,pairwise~group)
e<-summary(lsmeans(lm.asym, ~group))

Figure3<-ggplot(e, aes(x = group, y = lsmean, fill = group)) +
  geom_bar(aes(color = group), stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), color = "black", width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(-2.5,2.5) +
  labs(
    title = "GluCEST Asymmetry by Group",
    x = "Group",
    y = "GluCEST Asymmetry (% contrast)"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )+
  theme(legend.position = "none")

ggsave("~/Dropbox/Manuscripts_in_progress/Aging_GluCEST/final_analysis/data/Figure3.pdf",  Figure3,
       width=3.5, height=5)



#Glucest correlations with behavior
#correlations with MoCA and CCI

#MOCA vs left cest in HOA or HYA
cor.test(old_glu$LeftHippoGluCEST,old_glu$mocatotal)
cor.test(young_glu$LeftHippoGluCEST,young_glu$mocatotal)
#moca vs. right cest
cor.test(old_glu$RightHippoGluCEST,old_glu$mocatotal)
cor.test(young_glu$RightHippoGluCEST,young_glu$mocatotal)

#CCI vs left cest
cor.test(old_glu$LeftHippoGluCEST,old_glu$cci_total)
cor.test(young_glu$LeftHippoGluCEST,young_glu$cci_total)
#cci vs. right cest
cor.test(old_glu$RightHippoGluCEST,old_glu$cci_total)
cor.test(young_glu$RightHippoGluCEST,young_glu$cci_total)

#aysymm vs MoCA
cor.test(old_glu$LminusRasymmGluCEST,old_glu$mocatotal)
cor.test(young_glu$LminusRasymmGluCEST,young_glu$mocatotal)

