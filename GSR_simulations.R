library(MBESS)
library(MASS)
library(gdata)
library(psych)
library(plotrix)

rm(list=ls())

#function to return residuals of signal y (region) from signal x (GS)
regress <- function(x,y){
  m<-lm(y~x)
  return(resid(m))
}

#function to compute change in variances pre and post, with and without GSR
pre_post_vars_GSR <- function(ts,ts_post){
  vars_pre <- apply(ts,2,var)
  vars_post <- apply(ts_post,2,var)
  GS <- rowSums(ts)
  GS_post <- rowSums(ts_post)
  ts_GSR <- apply(ts,2,function(x) regress(GS,x))
  ts_post_GSR <- apply(ts_post,2,function(x) regress(GS_post,x))
  vars_pre_GSR <- apply(ts_GSR,2,var)
  vars_post_GSR <- apply(ts_post_GSR,2,var)
  return(list(vars_pre,vars_pre_GSR,vars_post,vars_post_GSR))
}

#function to generate synthetic independent data
#return vars pre and post GSR, before and after manipulation
varchange_independence <- function(vars, delta_vars){
  timepoints<-300
  ts<-matrix(data = NA,nrow=timepoints,ncol=length(vars))
  ts_post<-matrix(data = NA,nrow=timepoints,ncol=length(vars))
  for (i in 1:length(vars)){
    ts[,i]<-rnorm(timepoints, mean = 0, sd = sqrt(vars[i]))
    ts_post[,i]<-rnorm(timepoints, mean = 0, sd = sqrt((vars[i]+delta_vars[i])/vars[i]))
  }
  return(pre_post_vars_GSR(ts,ts_post))
}

#function to generate synthetic timeseries from asymmetric correlation mat + variances
timeseries <- function(cormat,vars){
  timepoints<-300
  #symmetrize with upper triangular
  lowerTriangle(cormat) <- upperTriangle(cormat, byrow=TRUE)
  diag(cormat)<-1
  #convert to covariance matrix
  covmat <- cor2cov(cormat,sqrt(vars))
  #generate simulated data
  #set.seed(1)
  ts <- mvrnorm(timepoints,replicate(length(vars),0),covmat)
  return(ts)
}

#function to generate synthetic covariances data
#returns vars pre and post GSR, before and after manipulation
synth_cormat <- function(rmat,deltar,vars,delta_vars){
  size<-length(vars)
  R <- matrix(runif(size^2), ncol=size,byrow=TRUE) 
  RtR <- R %*% t(R)
  cors_pre <- cov2cor(RtR)
  diag(cors_pre)<-0
  cors_pre<-fisherz2r(fisherz(cors_pre)/2)
  diag(cors_pre)<-1
  #replace amygdala's corvalues
  cors_pre[rmat>0]<-rmat[rmat>0]
  cors_post<-cors_pre+deltar
  #generate timeseries data for pre and post mats
  ts <- timeseries(cors_pre,vars)
  ts_post <- timeseries(cors_post,(vars+delta_vars)/vars)
  return(pre_post_vars_GSR(ts,ts_post))
}

##############
## set your variances and delta variances here
vars<-c(1,1,1,1,1,1)
delta_vars<-c(-0.3,-0.2,-0.15,-0.1,-0.05,0)

#have a zero-matrix of size 6 around
#use this instead of rmat and deltar arguments to generate random correlational structure
null<-replicate(6,replicate(6,0))

#z-values of 1, .75, .5, .25, 0 for amygdalo-connections
rmat<-matrix(c(
  0,0.7616,0.6351,0.4621,0.2449,0.0000001,
  0,0,0.7616,0.6351,0.4621,0.2449,
  0,0,0,0.7616,0.6351,0.4261,
  0,0,0,0,0.7616,0.6351,
  0,0,0,0,0,0.7616,
  0,0,0,0,0,0),
  nrow=length(vars),ncol=length(vars),byrow=TRUE)
lowerTriangle(rmat) <- upperTriangle(rmat, byrow=TRUE)

#reduce r^2 by .7
deltar<-matrix(c(
  0,-0.1244,-0.1037,-0.0755,-0.04,0,
  0,0,0,0,0,0,
  0,0,0,0,0,0,
  0,0,0,0,0,0,
  0,0,0,0,0,0,
  0,0,0,0,0,0),
  nrow=length(vars),ncol=length(vars),byrow=TRUE)
lowerTriangle(deltar) <- upperTriangle(deltar, byrow=TRUE)

perms<-500

##############
## GO!
##############

#condition 1: data is independent; only variances altered
#condition 2: data has positively skewed but random correlational structure; only variances altered
#condition 3: data has positively skewed, non-random correlational structure; only variances altered
#condition 4: data has positively skewed, non-random correlational structure; variances and correlations altered

#compute vars before and after manipulation pre and post GSR
vars_pre<-matrix(data = NA,nrow=length(vars),ncol=perms)
vars_pre_GSR<-matrix(data = NA,nrow=length(vars),ncol=perms)
vars_post<-matrix(data = NA,nrow=length(vars),ncol=perms)
vars_post_GSR<-matrix(data = NA,nrow=length(vars),ncol=perms)
delta_percvars<-matrix(data = NA,nrow=length(vars),ncol=perms)
delta_percvars_GSR<-matrix(data = NA,nrow=length(vars),ncol=perms)
for (i in 1:perms){
  
  #condition1
  #tmp4<-varchange_independence(vars, delta_vars)
  #condition2 
  #tmp4<-synth_cormat(null,null,vars,delta_vars)
  #condition3 
  #tmp4<-synth_cormat(rmat,null,vars,delta_vars)
  #condition4
  tmp4<-synth_cormat(rmat,deltar,vars,delta_vars)
  
  vars_pre[,i]<-tmp4[[1]]
  vars_pre_GSR[,i]<-tmp4[[2]]
  vars_post[,i]<-tmp4[[3]]
  vars_post_GSR[,i]<-tmp4[[4]]
  
  delta_percvars[,i]<-(vars_post[,i]-vars_pre[,i])/vars_pre[,i]
  delta_percvars_GSR[,i]<-(vars_post_GSR[,i]-vars_pre_GSR[,i])/vars_pre_GSR[,i]
}
vpre_ms<-apply(vars_pre,1,mean)
vpre_sds<-apply(vars_pre,1,std.error)
vpreGSR_ms<-apply(vars_pre_GSR,1,mean)
vpreGSR_sds<-apply(vars_pre_GSR,1,std.error)
vp_ms<-apply(vars_post,1,mean)
vp_sds<-apply(vars_post,1,std.error)
vpGSR_ms<-apply(vars_post_GSR,1,mean)
vpGSR_sds<-apply(vars_post_GSR,1,std.error)
dpv_ms<-apply(delta_percvars,1,mean)
dpv_sds<-apply(delta_percvars,1,std.error)
dpvGSR_ms<-apply(delta_percvars_GSR,1,mean)
dpvGSR_sds<-apply(delta_percvars_GSR,1,std.error)

plotstderrors <- function(x,y,sd){
  segments(x, y-sd,x, y+sd)
  epsilon = 0.02
  segments(x-epsilon,y-sd,x+epsilon,y-sd)
  segments(x-epsilon,y+sd,x+epsilon,y+sd)
}

plotme <- function(y1_m,y1_sd,y2_m,y2_sd,presetlim,ylabel){
  x<-seq(length(y1_m))
  plot (x, y1_m, ylim=presetlim,type="l",col="black",ylab=ylabel,xlab="node")
  plotstderrors(x,y1_m,y1_sd)
  #abline(1,0,col = "gray", lty = 2)
  lines(x, y2_m,col="red")
  plotstderrors(x,y2_m,y2_sd)
}

par(mfrow = c(1,2))
plotme(vp_ms,vp_sds,vpGSR_ms,vpGSR_sds,c(0.1,1.05),"raw variances")
lines(vpre_ms,col="black",lty=2)
lines(vpreGSR_ms,col="red",lty=2)
plotme(dpv_ms,dpv_sds,dpvGSR_ms,dpvGSR_sds,c(-0.35,0.1),"%∆ variance (pre vs. post)")

### visualize network in igraph ###
#g<-graph.adjacency(rmat,weighted=TRUE,mode="undirected",diag=FALSE)
#lo<-layout.spring(g,dim=2)
#plot(g, layout=lo, edge.width=E(g)$weight*5, edge.label=round(E(g)$weight, 3),
#     vertex.size=30, edge.arrow.size=2,rescale=FALSE,xlim=range(lo[,1]),
#     ylim=range(lo[,2]), vertex.label.dist=1, vertex.label.color="red")