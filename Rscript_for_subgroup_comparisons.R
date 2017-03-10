#-------------------
library(gtools)
library(lavaan)
library(Hmisc)
library(psych)
library(mice)
library(CCA)
library(corrplot)
library(GenABEL)
#-------------------
#if in lab
#setwd("/Users/grayson/Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")
#if from home
setwd("Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")

rm(list=ls())
d <- read.csv("7_measures.csv")
cms_subj_list <- read.csv("../../connectome_csvs/123connectomes_scale33parcellation/SUBS_123maleids.txt",header=FALSE)
cms_subgroups <- read.csv("../../connectome_csvs/ASD_2DWI_RF_subtypes_123subs.csv",header=FALSE)

#---------------------------------------------------------------------------------
#truncate the dataset:
#choose only males time 1's
trunc <- d
trunc <- trunc[trunc$subject_gender=="Male",]
trunc <- trunc[trunc$visit=="1",]

#------------------------------------------------------------------------------------------------
#match up the truncated dataset with the connectome data
trunc_with_connectomes<-match(as.character(cms_subj_list[,]),as.character(trunc$Case))
trunc<-trunc[trunc_with_connectomes,]
#---------------------------------------------------------------------------------
#choose only ASD #don't fuck up the ordering here
cms_subgroups <- as.data.frame(cms_subgroups[trunc$app_diag=="ASD",])
cms_subj_list <- as.data.frame(cms_subj_list[trunc$app_diag=="ASD",])
trunc <- trunc[trunc$app_diag=="ASD",]

#---------------------------------------------------------------------------------
#extract behavioral measures

#extract the total SRS, the 4 social-specific SRS measures, the average of the 4, and the mannerism
SRS_total <- trunc$TOTAL_T
SRS_4socials <- cbind(trunc$AWR_T,trunc$COG_T,trunc$COM_T,trunc$MOT_T)
SRS_AVGsocial <- rowMeans(SRS_4socials)
SRS_manner <- trunc$MANNER_T

#more social
ADOS_sa <- trunc$sa_total #this is the one that combines the imagination thing
ADOS_cst <- trunc$ados_cst #this is just social interaction + communication deficits
ADI_social <- trunc$SOCT_CS

#extract the section B scores (communication) of the ADI
ADI_2comms <- cbind(trunc$COMNVTCS,trunc$COMVT_CS)
ADI_2comms[is.na(ADI_2comms)] <- 0
ADI_comm <- rowSums(ADI_2comms)

#more communication
PPVT<-trunc$SCORE_STANDARD * -1
EOWPVT<-trunc$EOWPVT_SS * -1 #this is a censored measure

#extract RBS measures
ADOS_rrb <- trunc$rrb_total
RBS <- trunc$RBS_T_O
ADI_rbeh <- trunc$BEHT_CS

#and IQ and ADOS severity
#VDQ <- trunc$VDQ * -1
#NVDQ <- trunc$NVDQ * -1
#DQ <- cbind(VDQ,NVDQ)
DQ<-trunc$DQ * -1
ADOS_severity <- trunc$ados_severity

#short sensory profile
sensory<-trunc$TOT_RS * -1

#CBC
CBC_internalizing<-trunc$internalizing_t
CBC_externalizing<-trunc$externalizing_t
CBC_somatic<-trunc$somatic_complaints_t
CBC_attention<-trunc$attention_problems_t
CBC_aggressive<-trunc$aggressive_bx_t
CBC_reactive<-trunc$emotionally_reactive_t
CBC_sleep<-trunc$sleep_problems_t
CBC_anxiety<-trunc$anxious_depressed_t

#concatenate observed variables
#CBC<-cbind(CBC_externalizing,CBC_internalizing,CBC_aggressive,CBC_anxiety,CBC_attention,CBC_reactive,CBC_sleep,CBC_somatic)
CBC<-cbind(CBC_externalizing,CBC_internalizing)
social_obs <- cbind(SRS_AVGsocial)
rbs_obs <- cbind(SRS_manner,RBS)

#use just the ADOS and ADI here: 
#obs <- cbind(ADOS_cst,ADOS_rrb,ADI_social,ADI_comm,ADI_rbeh)
#ordering here makes sure ADOS and ADI are first variables for vizualization
#obs <- cbind(ADOS_cst,ADI_social,ADI_comm,ADOS_rrb,ADI_rbeh,social_obs,rbs_obs,PPVT,DQ,CBC,sensory)

#obs <- cbind(ADI_social,ADI_comm,ADI_rbeh)
#obs <- cbind(social_obs,rbs_obs,CBC,sensory)
obs <- cbind(PPVT,DQ)

#----------------------------------------------------------------------
#set here to test factor scores instead of using individual measures
#obs <- read.csv("ASDmales_89subs_imputed_6factors_try2.csv",header=FALSE)

#----------------------------------------------------------------------
#compare subgroups
subgroups <- list()
measure_t_stat <- list()
for (i in 1:max(cms_subgroups)){
  subgroups[[i]] <- obs[cms_subgroups==i,]
}
for (j in 1:length(obs[1,])){
  sg1_measure <- subgroups[[1]][,j]
  sg2_measure <- subgroups[[2]][,j]
  measure_t_stat[[j]] <- t.test(na.omit(sg1_measure), na.omit(sg2_measure))
}
measure_t_stat