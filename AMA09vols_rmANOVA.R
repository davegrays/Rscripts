#clear all
rm(list=ls())
#libraries that need to be loaded: car, nlme
#name of working folder
setwd('/Users/dgrayson1/Dropbox/LABwork_amaral/labstuff_results/AMA09/R_analyses/')

#which data to load and how to save?
infile='FA_CCmsfg.txt'
outfile='statsum2_FA_CCmsfg.txt'
#'raw_volume_ratios.txt','corrected_volume_ratios.txt'
#'FA_uncinate.txt','FA_cingulum.txt','FA_CCmidcing.txt','FA_CCmsfg.txt'

#how are subjects arranged? Careful!!
#groups=c(rep("Control",8),rep("Lesion",6)) #for volumes files
groups=c(rep("Lesion",6),rep("Control",8)) #for FA files

#######################
#######################
#sink(outfile) #this pipes the script's output to the outfile

#load data
volco=read.table(infile)

#convert to format useable with rm-ANOVA
volco2=stack(volco)
#put subject names in
subs=rownames(volco)
volco2$subjects=subs
#put groups in
volco2$groups=groups

#reset column names
colnames(volco2)=c('volume','region','subject','group')

#run rm-ANOVA (Type I SS)
summary(aov(volume ~ group * region + Error(subject/region), data=volco2))

#attempting to use Type III SS
#Anova(lme(volume ~ group*region, random=~1 | subject, method="ML", data=volco2),type=3)
#Anova(lmer(volume ~ group*region + (1|subject), data=volco2),type=3,test.statistic="F")

#and type II SS
m1=lme(volume ~ group*region, random=~1 | subject, method="ML", data=volco2)
m2=lmer(volume ~ group*region + (1|subject), data=volco2)
Anova(m1,type=2) #chisq
Anova(m2,type=2,test.statistic="F") #F

#still not the same as SPSS; are the above mixed effects models specified correctly? something weird