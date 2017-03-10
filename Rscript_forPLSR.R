#-------------------
library(CCA)
library(lavaan)
library(Hmisc)
library(psych)
library(mice)
library(pls)
library(corrplot)
library(GenABEL)
library(spsl)
library(glmnet)
#-------------------

#if in lab
#setwd("/Users/grayson/Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")
#if from home
setwd("Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")

rm(list=ls())
d <- read.csv("7_measures.csv")
cms_subj_list <- read.csv("../../connectome_csvs/123connectomes_scale33parcellation/SUBS_123maleids.txt",header=FALSE)
connectomes <- read.csv("../../connectome_csvs/123connectomes_scale33parcellation/connectomes_2D_csvs/ranknormal_2Dcsv_123connectome_scale33_nuisanceregressed.csv",header=FALSE)
ROIpairs <- read.csv("../../connectome_csvs/123connectomes_scale33parcellation/connectomes_2D_csvs/connection_names_scale33.csv",header=FALSE)
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

#------------------------------------------------------------------------------------------------

##choose only ASD #don't switch the order of cms_subj_list/connectomes/trunc here
cms_subj_list <- as.data.frame(cms_subj_list[trunc$app_diag=="ASD",])
#cms_subj_list <- cms_subj_list[trunc$version=="ADOS-G",]
connectomes <- as.data.frame(connectomes[trunc$app_diag=="ASD",])
#connectomes <- connectomes[trunc$version=="ADOS-G",]
trunc <- trunc[trunc$app_diag=="ASD",]
#trunc <- trunc[trunc$version=="ADOS-G",]

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
ADOS_severity <- trunc$ados_severity

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
social_obs <- cbind(SRS_AVGsocial)
rbs_obs <- cbind(SRS_manner,RBS)

#obs <- cbind(social_obs,rbs_obs,ADOS_severity)
#ordering here makes sure ADOS and ADI are first variables for vizualization
obs <- cbind(ADOS_cst,ADI_social,ADI_comm,ADOS_rrb,ADI_rbeh,social_obs,rbs_obs,PPVT,DQ,CBC,sensory)
#obs <- cbind(ADOS_cst,ADOS_rrb)
#obs <- cbind(ADI_social,ADI_comm,ADI_rbeh)
#obs <- cbind(social_obs,rbs_obs,CBC,sensory)
#obs <- cbind(PPVT,DQ)

#------------------------------------------------------------------------------------------------
#set here to test factor scores instead of using individual measures
obs <- read.csv("ASDmales_89subs_imputed_6factors_try2.csv",header=FALSE)
obs <- cbind(obs[,5]) #rbeh
#obs <- cbind(obs[,3]) #social

#------------------------------------------------------------------------------------------------
#transform connectome data using "log" or "ranknormal"
#source("transform_connectomes.R")
#trconnectomes <- transform_connectomes(connectomes, transform_type="ranknormal")
#plot(trconnectomes[[2]],trconnectomes[[3]],xlab="norm test pre-transformation; -log pval",ylab="norm test post-transformation;-log pval")
#connectomes <- trconnectomes[[1]]

#---------------------------------------------------------

#run PLSR!!
Nsamples <- length(obs[,1]) #sample size
Nmeasures <- length(obs[1,]) #number of behavioral measures
number_of_comps <- 10 #how many pls components to estimate

#impute missing behav data
if (sum(is.na(obs)) > 0) {
  obs_complete<-complete(mice(obs),5)
} else {
  obs_complete <-  as.data.frame(obs)
}
obs <- as.data.frame(obs)
#rank transform behav data
#trobs_complete <- transform_connectomes(obs_complete, transform_type="ranknormal")
#obs_complete <- trobs_complete[[1]]

####### classic plsr #########
#concatenate connectome and (imputed) behavior data into a single dataframe
#mydata <- as.data.frame(matrix(,nrow = Nsamples,ncol = 0))
#mydata$behavior <- I(as.matrix(obs_complete))
#mydata$connectomes <- I(as.matrix(connectomes))

#pls <- plsr(behavior ~ connectomes, ncomp = number_of_comps, data = mydata, scale = "TRUE", validation = "LOO")
#summary(pls)
#plot(RMSEP(pls),legendpos="topright")
#plot(pls, ncomp = 1, asp=1, line = TRUE)
#plot(pls, plottype = "scores", comps = 1:5)

####### sparse plsr #########
#cv <- cv.spls( connectomes, obs, eta = seq(0.001,0.901,.1), K = c(1:3) )
#f <- spls( connectomes, obs, eta = 0.5, K = 1 )
#print(f)
#plot.spls(f)

####### sparse regularized multiple regression #########

x_pred <- array(0,Nsamples)
for (i in seq(Nsamples)){ 
  x_test <- as.matrix(connectomes[i,])
  y_test <- as.matrix(obs[i,])
  x_train <- as.matrix(connectomes[-i,])
  y_train <- as.matrix(obs[-i,])
  
  fit = glmnet(x_train, y_train)
  x_pred[i] <- predict(fit, x_test, type = "response", s = 0.25)
}
cor(x_pred,as.matrix(obs))


x<- as.matrix(connectomes)
y<-as.matrix(obs)
#fit = glmnet(x,y)
cvfit = cv.glmnet(x,y, type.measure = "mse", nfolds = Nsamples)
#predict(cvfit,x,s=0.25)
plot(cvfit)
#cvfit$lambda.1se
#cvfit$lambda.min
