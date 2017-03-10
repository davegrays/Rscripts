##############################
library(lavaan)
library(Hmisc)
library(psych)
library(corrplot)
library(FactMixtAnalysis)
library(mice)
##############################

#if in lab
setwd("/Users/grayson/Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")
#if from home
#setwd("Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")

rm(list=ls())
d <- read.csv("7_measures.csv")
cms_subj_list <- read.csv("SUBS_89ASDmaleids.txt",header=FALSE)
cms_subgroups <- read.csv("../../connectome_csvs/ASD_2DWI_RF_subtypes_89subs.csv",header=FALSE)

num_factors <- 5

#---------------------------------------------------------------------------------
#truncate the dataset:
#choose only males, time 1's
trunc <- d
trunc <- trunc[trunc$subject_gender=="Male",]
trunc <- trunc[trunc$visit=="1",]

#---------------------------------------------------------------------------------
#only include subjects with connectome data; match to connectome ordering
#trunc_with_connectomes<-match(as.character(cms_subj_list[,]),as.character(trunc$Case))
#trunc<-trunc[trunc_with_connectomes,]

#---------------------------------------------------------------------------------
#choose only ASD #don't fuck up the ordering here
#cms_subgroups <- as.data.frame(cms_subgroups[trunc$app_diag=="ASD",])
#cms_subj_list <- as.data.frame(cms_subj_list[trunc$app_diag=="ASD",])
trunc <- trunc[trunc$app_diag=="ASD",]

#---------------------------------------------------------------------------------
#extract measures

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
#EOWPVT<-trunc$EOWPVT_SS * -1 #this is a censored measure

#extract RBS measures
ADOS_rrb <- trunc$rrb_total
RBS <- trunc$RBS_T_O
ADI_rbeh <- trunc$BEHT_CS

#and IQ and ADOS severity
VDQ <- trunc$VDQ * -1
NVDQ <- trunc$NVDQ * -1
#DQ <- cbind(VDQ,NVDQ)
DQ<-trunc$DQ * -1
#ADOS_severity <- trunc$ados_severity #also a censored measure

#short sensory profile
sensory<-trunc$TOT_RS * -1

#CBC
CBC_internalizing<-trunc$internalizing_t
CBC_externalizing<-trunc$externalizing_t
#CBC_somatic<-trunc$somatic_complaints_t
#CBC_attention<-trunc$attention_problems_t
#CBC_aggressive<-trunc$aggressive_bx_t
#CBC_reactive<-trunc$emotionally_reactive_t
#CBC_sleep<-trunc$sleep_problems_t
#CBC_anxiety<-trunc$anxious_depressed_t

#concatenate observed variables
CBC<-cbind(CBC_externalizing,CBC_internalizing)
social_obs <- cbind(SRS_AVGsocial,ADOS_cst,ADI_social)
rbs_obs <- cbind(SRS_manner,ADI_rbeh,RBS)
#obs <- cbind(social_obs,rbs_obs,CBC,sensory)
obs <- cbind(social_obs,rbs_obs,CBC,sensory)

#------------------------------------------------------------------------------------------------
#impute missing data
obs_complete<-complete(mice(obs),5)
obs <- as.data.frame(obs)

#visualize correlations between observed variables:
corrplot(cor(obs_complete), method = "number", order = "hclust")
#corrplot(cor(obs_complete), method = "shade", order = "hclust")

#------------------------------------------------------------------------------------------------
#exploratory factor analyze on the imputed or non-imuted data
facm <- factanal(obs_complete,factors = num_factors,rotation="varimax", scores="regression")
print(facm, digits = 2, cutoff = .2, sort = TRUE)
facscores <- facm$scores

#------------------------------------------------------------------------------------------------
#confirmatory factor analyze on the imputed or non-imputed data
cfamodel <- ' parent  =~  sensory + CBC_externalizing+ SRS_AVGsocial + RBS + CBC_internalizing 
              social =~ ADI_social + SRS_AVGsocial + ADOS_cst
              rbeh   =~ ADI_rbeh + RBS
              ADI_comm ~~ ADI_social
              ADI_comm ~~ ADI_rbeh'
              #SRS_AVGsocial ~~ SRS_manner'
cfafit <- cfa(cfamodel, data = as.data.frame((obs_complete)))




#------------------------------------------------------------------------------------------------
#match up the truncated dataset with the connectome data
trunc_with_connectomes<-match(as.character(cms_subj_list[,]),as.character(trunc$Case))
trunc<-trunc[trunc_with_connectomes,]
obs<-obs[trunc_with_connectomes,]
obs_complete<-obs_complete[trunc_with_connectomes,]
facscores<-facscores[trunc_with_connectomes,]





#------------------------------------------------------------------------------------------------
#write fac scores to file
#write.csv(facscores,"ASDmales_89subs_imputed_6factors_TRY2.csv")




#------------------------------------------------------------------------------------------------
#run tests of imaging-subgroup differences on factor scores:
subgroups <- list()
measure_t_stat <- list()
for (i in 1:max(cms_subgroups)){
  subgroups[[i]] <- facscores[cms_subgroups==i,]
}
for (j in 1:length(facscores[1,])){
  sg1_measure <- subgroups[[1]][,j]
  sg2_measure <- subgroups[[2]][,j]
  measure_t_stat[[j]] <- t.test(na.omit(sg1_measure), na.omit(sg2_measure))
}

#------------------------------------------------------------------------------------------------
#repeat using 1/2 test and replication datasets
tmp<-na.omit(obs_complete)
lt<-length(tmp[,1])
test<-tmp[1:(lt/2),]
replication<-tmp[(lt/2 + 1):lt,]
facm <- factanal(test,factors = 6,rotation="varimax")
print(facm, digits = 2, cutoff = .2, sort = TRUE)
facm <- factanal(replication,factors = 6,rotation="varimax")
print(facm, digits = 2, cutoff = .2, sort = TRUE)

#------------------------------------------------------------------------------------------------
#do factor mixture analysis:
fmam.2.2 <- fma(obs_complete,2,2, scaling = TRUE)
fmam.2.3 <- fma(obs_complete,2,3, scaling = TRUE)
fmam.3.2 <- fma(obs_complete,3,2, scaling = TRUE)
fmam.3.3 <- fma(obs_complete,3,3, scaling = TRUE)
fmam.2.4 <- fma(obs_complete,2,4, scaling = TRUE)
fmam.3.4 <- fma(obs_complete,3,4, scaling = TRUE)
fmam.1.4 <- fma(obs_complete,1,4, scaling = TRUE)
fmam.1.3 <- fma(obs_complete,1,3, scaling = TRUE)
fmam.1.2 <- fma(obs_complete,1,2, scaling = TRUE)
fmam.1.5 <- fma(obs_complete,1,5, scaling = TRUE)
fmam.2.5 <- fma(obs_complete,2,5, scaling = TRUE)
fmam.3.5 <- fma(obs_complete,3,5, scaling = TRUE)

fmam.3.5$aic;fmam.3.4$aic;fmam.3.3$aic;fmam.2.5$aic;fmam.2.4$aic;fmam.2.3$aic;fmam.1.5$aic;fmam.1.4$aic;fmam.1.3$aic
fmam.3.5$bic;fmam.3.4$bic;fmam.3.3$bic;fmam.2.5$bic;fmam.2.4$bic;fmam.2.3$bic;fmam.1.5$bic;fmam.1.4$bic;fmam.1.3$bic