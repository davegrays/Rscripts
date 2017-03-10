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
setwd("/Users/grayson/Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")
#if from home
#setwd("Box Sync/WORK/Amaral/APP/Behavior_stuff/my_factor_analysis/")

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
obs <- cbind(ADOS_cst,ADI_social,ADI_comm,ADOS_rrb,ADI_rbeh,social_obs,rbs_obs,PPVT,DQ,CBC,sensory)

#obs <- cbind(ADI_social,ADI_comm,ADI_rbeh)
#obs <- cbind(social_obs,rbs_obs,CBC,sensory)
#obs <- cbind(PPVT,DQ)

#------------------------------------------------------------------------------------------------

##save imputed behavior (no missing data) matched to connectome order for Jeremy
#obs_complete<-complete(mice(obs),5)
#write.csv(obs_complete,"behavior_data_ASDkids.csv")

#------------------------------------------------------------------------------------------------

#transform connectome data using "log" or "ranknormal"
source("transform_connectomes.R")
trconnectomes <- transform_connectomes(connectomes, transform_type="ranknormal")
plot(trconnectomes[[2]],trconnectomes[[3]],xlab="norm test pre-transformation; -log pval",ylab="norm test post-transformation;-log pval")
abline(4.6,0)
connectomes <- trconnectomes[[1]]

#---------------------------------------------------------

#impute missing behavioral data
obs_complete<-complete(mice(obs),5)

#rank-normal transform behavioral data
#source("transform_connectomes.R")
#trobs <- transform_connectomes(obs_complete, transform_type="ranknormal")
#obs_complete <- trobs[[1]]

#run CCA using different amounts of regularization on the connectome data
#cc1 <- rcc(obs_complete,connectomes,0,1)
#cc0.5 <- rcc(obs_complete,connectomes,0,0.5)
#cc0.1 <- rcc(obs_complete,connectomes,0,0.1)
#cc0.01 <- rcc(obs_complete,connectomes,0,0.01)
#cc0.001 <- rcc(obs_complete,connectomes,0,0.001)
#cc0.0001 <- rcc(obs_complete,connectomes,0,0.0001)

#for each regularization parameter, visualize the loadings of the behavioral measures onto each CV
#corrplot(cc1$xcoef, method = "shade")
#corrplot(cc0.1$xcoef, method = "shade")
#corrplot(cc0.01$xcoef, method = "shade")
#corrplot(cc0.001$xcoef, method = "shade")
#corrplot(cc0.0001$xcoef, method = "shade")
#you can use 0.001 as a default

#------------------------------------------------------------------------
#
# I/O for CCA analysis
#

obs_missing <- obs_complete
Nsamples <- length(obs_complete[,1]) #sample size
Nmeasures <- length(obs_complete[1,]) #number of behavioral measures
number_of_pcs <- 15 #how many top principal components of connectome data to use
number_of_perms <- 1000 #how many permutations for signifcance testing
  
# GO!

#run CCA using PCA for dim-reduction on connectome data first
pcs <- prcomp(connectomes)
trunc_pcs <- pcs$x[,1:number_of_pcs]
#run cca between brain and behav (using the incomplete behav data)
ccs_missing <- CCA::cc(obs_missing,trunc_pcs)

#generate null distribution of maximal CV correlations by permuting the rows
#in obs and trunc_pcs, rerunning CCA, and extracting the maximal CV correlation
#(i.e. the correlation for the first CV) per null run
permvec <- seq(1,Nsamples)
corrs_null <- matrix(, nrow = number_of_perms, ncol = 1)
for (i in 1:number_of_perms){
  nullvec <- permute(permvec)
  nullvec2 <- permute(permvec)
  ccs_null <- CCA::cc(obs_missing[nullvec,],trunc_pcs[nullvec2,])
  corrs_null[i] <- ccs_null$cor[1]
}

#estimate significance of observed CV correlations
number_of_CVs <- min(number_of_pcs,Nmeasures)
corrected_p_CV <- matrix(, nrow = 1, ncol = number_of_CVs)
for (i in 1:number_of_CVs){
  corrected_p_CV[i] <- mean(ccs_missing$cor[i] < corrs_null)
}
corrected_p_CV
ccs_missing$cor

#visualize the loadings of the behavioral measures on to each CV
CV_loadings <- ccs_missing$scores
corrplot::corrplot(CV_loadings$corr.X.yscores, method = "shade")

