.libPaths("/gpfs/gibbs/project/wang_zuoheng/gx28/R_lib")
options(bitmapType='cairo')
library(jomo)
library(xtable)
library(gtools)

setwd("/gpfs/gibbs/project/wang_zuoheng/gx28/missing_indicator/results_noMinMI_2k2k_checkconverg")
#simu_results <- read.csv("simulation_combined.csv",header = T)
para_args <- read.table("/gpfs/gibbs/project/wang_zuoheng/gx28/missing_indicator/param_args_valid_rmoutlier_v2.txt",header=F)
rdslist=list.files(path=".", pattern=".csv")
rdslist=  mixedsort(sort(rdslist))
length(rdslist)
nsim <- 2000
nscena <- 4704

#nscena <- dim(simu_results)[1]/2000
#simu_results <- simu_results[simu_results$pr_cmiss>0.05,]

#table(simu_results$rr_e)
#simu_results <- simu_results[simu_results$rr_e<1,]
#table(simu_results$rr_e)

#Case 1, bias=(beta-log(rr_e))/log(rr_e) if rr_e!=1, bias=beta otherwise, keep all
missing_cutoff <- 1999

sample_size <-  rep(0,nscena)
mean_beta <- matrix(0,nscena,5)
mean_capture <- matrix(0,nscena,4)
missing_sum <- rep(0,nscena)
missing_inc_ind <- rep(0,nscena)
missing_rate <- rep(0,4)
var_pro <- NULL
i=1
missing_rate_table=NULL
for(i in 1:nscena){
  index_i <- (i-1)*nsim+1:nsim
  sub_simu <- read.csv(paste0("./",rdslist[[i]]))
  missing_index <- c(which(is.na(sub_simu$beta_nm)),which(is.na(sub_simu$beta_mcim)),
                     which(is.na(sub_simu$beta_mimp)),which(is.na(sub_simu$beta_cc)))
  missing_ratei=c(length(which(is.na(sub_simu$beta_nm))),length(which(is.na(sub_simu$beta_mcim))),
                  length(which(is.na(sub_simu$beta_mimp))),length(which(is.na(sub_simu$beta_cc))))
  missing_rate <- missing_rate + c(length(which(is.na(sub_simu$beta_nm))),length(which(is.na(sub_simu$beta_mcim))),
                                   length(which(is.na(sub_simu$beta_mimp))),length(which(is.na(sub_simu$beta_cc))))
  missing_rate_table=rbind(missing_rate_table, c(missing_ratei,para_args[i,]))
  missing_index <- unique(missing_index)
  missing_sum[i] <- length(missing_index)
  if(missing_sum[i]<=missing_cutoff){
    missing_inc_ind[i] <- 1
    if(missing_sum[i]>0) {sub_simu <- sub_simu[-missing_index,]}
    sub_beta <- cbind(sub_simu$beta_nm,sub_simu$beta_mcim,sub_simu$beta_mimp,sub_simu$beta_cc)
    capture <- cbind(sub_simu$capture_nm,sub_simu$capture_mcim,sub_simu$capture_mimp,sub_simu$capture_cc)
    mean_beta[i,] <- c(apply(sub_beta,2,mean),sub_simu$rr_e[1])
    mean_capture[i,] <- apply(capture,2,mean)
    sample_size[i] <- dim(sub_beta)[1]
    #    var_pro <- rbind(var_pro,cbind((sub_simu$beta_nm_se/sub_simu$beta_mcim_se)^2,
    #                                   (sub_simu$beta_nm_se/sub_simu$beta_mimp_se)^2,
    #                                   (sub_simu$beta_nm_se/sub_simu$beta_cc_se)^2))
    #var_pro <- rbind(var_pro,cbind((sub_simu$beta_nm_se)^2,(sub_simu$beta_mcim_se)^2,
    #                               (sub_simu$beta_mimp_se)^2,
    #                               (sub_simu$beta_cc_se)^2))
  }
  
}
missing_rate_table[1:40,]

missing_rate/(nsim*nscena)

sum(missing_inc_ind)
mean_beta_sub <- mean_beta[missing_inc_ind==1,]
mean_capture_sub <- mean_capture[missing_inc_ind==1,]
sample_size <- sample_size[missing_inc_ind==1]
uni_rre <- unique(mean_beta_sub[,5])
bias_sub <- mean_beta_sub[,1:4]

for(i in 1:dim(mean_beta_sub)[1]){
  if(log(mean_beta_sub[i,5])!=0) {bias_sub[i,] <- (bias_sub[i,]-log(mean_beta_sub[i,5]))/log(mean_beta_sub[i,5])
  }else{bias_sub[i,] <- bias_sub[i,]
  }
}

beta_table <- cbind(apply(bias_sub,2,mean),
                    apply(bias_sub,2,median),
                    apply(bias_sub,2,function(x) length(which(abs(x)>0.05))/length(x)),
                    apply(bias_sub,2,function(x) length(which(abs(x)>0.1))/length(x)),
                    apply(bias_sub,2,function(x) length(which(abs(x)>0.5))/length(x)))
beta_table
length(bias_sub[,1])
var_pro_prop <- var_pro[,1]/var_pro[,-1]
var_table <- cbind(apply(var_pro_prop,2,mean),
                   apply(var_pro_prop,2,median))
var_table
var_pro_prop2 <- var_pro[,-1]/var_pro[,1]
var_table2 <- cbind(apply(var_pro_prop2,2,mean),
                   apply(var_pro_prop,2,median))
var_table2
#Case 2, bias=(beta-log(rr_e))/log(rr_e) if rr_e!=1, bias=beta otherwise, keep missing_cutoff < 1000
missing_cutoff <- 2500

sample_size <-  rep(0,nscena)
mean_beta <- matrix(0,nscena,5)
mean_capture <- matrix(0,nscena,4)
missing_sum <- rep(0,nscena)
missing_inc_ind <- rep(0,nscena)
var_pro2 <- NULL
missing_rate <- rep(0,4)
for(i in 1:nscena){
  index_i <- (i-1)*nsim+1:nsim
  sub_simu <- read.csv(paste0("./",rdslist[[i]]))
  missing_index <- c(which(is.na(sub_simu$beta_nm)),which(is.na(sub_simu$beta_mcim)),
                     which(is.na(sub_simu$beta_mimp)),which(is.na(sub_simu$beta_cc)))
  missing_index <- unique(missing_index)
  missing_sum[i] <- length(missing_index)
  if(missing_sum[i]<=missing_cutoff){
    missing_inc_ind[i] <- 1
    if(missing_sum[i]>0) {sub_simu <- sub_simu[-missing_index,]}
    #sub_simu <- sub_simu[1:2000,]
    sub_beta <- cbind(sub_simu$beta_nm,sub_simu$beta_mcim,sub_simu$beta_mimp,sub_simu$beta_cc)
    capture <- cbind(sub_simu$capture_nm,sub_simu$capture_mcim,sub_simu$capture_mimp,sub_simu$capture_cc)
    mean_beta[i,] <- c(apply(sub_beta,2,mean),sub_simu$rr_e[1])
    mean_capture[i,] <- apply(capture,2,mean)
    sample_size[i] <- dim(sub_beta)[1]
    #    var_pro <- rbind(var_pro,cbind((sub_simu$beta_nm_se/sub_simu$beta_mcim_se)^2,
    #                                   (sub_simu$beta_nm_se/sub_simu$beta_mimp_se)^2,
    #                                   (sub_simu$beta_nm_se/sub_simu$beta_cc_se)^2))
    var_pro2 <- rbind(var_pro2,cbind((sub_simu$beta_nm_se)^2,(sub_simu$beta_mcim_se)^2,
                                     (sub_simu$beta_mimp_se)^2,
                                     (sub_simu$beta_cc_se)^2))
  }
}
missing_rate/(nsim*nscena)

sum(missing_inc_ind)
mean_beta_sub <- mean_beta[missing_inc_ind==1,]
mean_capture_sub <- mean_capture[missing_inc_ind==1,]
sample_size <- sample_size[missing_inc_ind==1]
uni_rre <- unique(mean_beta_sub[,5])
bias_sub <- mean_beta_sub[,1:4]

for(i in 1:dim(mean_beta_sub)[1]){
  if(log(mean_beta_sub[i,5])!=0) {bias_sub[i,] <- (bias_sub[i,]-log(mean_beta_sub[i,5]))/log(mean_beta_sub[i,5])
  }else{bias_sub[i,] <- bias_sub[i,]
  }
}


beta_table2 <- cbind(apply(bias_sub,2,mean),
                     apply(bias_sub,2,median),
                     apply(bias_sub,2,function(x) length(which(abs(x)>0.05))/length(x)),
                     apply(bias_sub,2,function(x) length(which(abs(x)>0.1))/length(x)))
beta_table2
var_pro3 <- var_pro2[,1]/var_pro2[,-1]
var_table2 <- cbind(apply(var_pro3,2,mean),
                    apply(var_pro3,2,median))
var_table2
var_pro4 <- var_pro2[,-1]/var_pro2[,1]
var_table3 <- cbind(apply(var_pro4,2,mean),
                    apply(var_pro4,2,median))
var_table3

beta_table_all <- cbind(beta_table,beta_table2)
rownames(beta_table_all) <- c("NM","MCIM","MI","CC")
beta_table_all[,c(1,2,5,6)]<- beta_table_all[,c(1,2,5,6)]*100
xtable(beta_table_all)

var_table_all <- cbind(var_table,var_table2) 
rownames(var_table_all) <- c("MCIM","MI","CC")
xtable(var_table_all,digits=2)
range(var_pro[,3])

apply(var_pro,2,median)
summary(var_pro[,1])
summary(var_pro[,2])
summary(var_pro[,3])
summary(var_pro[,4])
boxplot(var_pro[,1])
which(var_pro[,1]>100000)
subvar_pro <- var_pro[-which(var_pro[,1]>100000),]
apply(subvar_pro,2,mean)
boxplot(subvar_pro[,1])

summary(subvar_pro[,1])
summary(subvar_pro[,2])
summary(subvar_pro[,3])
summary(subvar_pro[,4])
length(which(var_pro2[,1]>100000))




.libPaths("/gpfs/gibbs/project/wang_zuoheng/gx28/R_lib")
options(bitmapType='cairo')
library(jomo)

setwd("/gpfs/gibbs/project/wang_zuoheng/gx28/missing_indicator/results_noMinMI_2k2k_checkconverg")
#simu_results <- read.csv("simulation_combined.csv",header = T)
rdslist=list.files(path=".", pattern=".csv")
length(rdslist)
nsim <- 2000
nscena <- 4704

#nscena <- dim(simu_results)[1]/2000
#simu_results <- simu_results[simu_results$pr_cmiss>0.05,]

#table(simu_results$rr_e)
#simu_results <- simu_results[simu_results$rr_e<1,]
#table(simu_results$rr_e)

#Case 1, bias=(beta-log(rr_e))/log(rr_e) if rr_e!=1, bias=beta otherwise, keep all
missing_cutoff <- 1999


missing_rate <- rep(0,4)
nscena_sum <- rep(0,4)
bias_nm <- NULL
bias_mcim <- NULL
bias_mimp <- NULL
bias_cc <- NULL
var_nm <- NULL

i=156
for(i in 1:nscena){
  index_i <- (i-1)*nsim+1:nsim
  sub_simu <- read.csv(paste0("./",rdslist[[i]]))
  missing_index_nm <- which(is.na(sub_simu$beta_nm))
  missing_index_mcim <-  which(is.na(sub_simu$beta_mcim))
  missing_index_mimp <-  which(is.na(sub_simu$beta_mimp))
  missing_index_cc <-  which(is.na(sub_simu$beta_cc))
  missing_rate <- missing_rate + c(length(which(is.na(sub_simu$beta_nm))),length(which(is.na(sub_simu$beta_mcim))),
                                   length(which(is.na(sub_simu$beta_mimp))),length(which(is.na(sub_simu$beta_cc))))

  if(length(missing_index_nm)<=missing_cutoff){
    nscena_sum[1] <- nscena_sum[1] + 1
    if(sub_simu$rr_e[1]!=1) {bias_nm <- c(bias_nm,(mean(sub_simu$beta_nm,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_nm <- c(bias_nm,mean(sub_simu$beta_nm,na.rm=T))
    }
  }
  if(length(missing_index_mcim)<=missing_cutoff){
    nscena_sum[2] <- nscena_sum[2] + 1
    if(sub_simu$rr_e[1]!=1) {bias_mcim <- c(bias_mcim,(mean(sub_simu$beta_mcim,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_mcim <- c(bias_mcim,mean(sub_simu$beta_mcim,na.rm=T))
    }
  }
  if(length(missing_index_mimp)<=missing_cutoff){
    nscena_sum[3] <- nscena_sum[3] + 1
    if(sub_simu$rr_e[1]!=1) {bias_mimp <- c(bias_mimp,(mean(sub_simu$beta_mimp,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_mimp <- c(bias_mimp,mean(sub_simu$beta_mimp,na.rm=T))
    }
  }
  if(length(missing_index_cc)<=missing_cutoff){
    nscena_sum[4] <- nscena_sum[4] + 1
    if(sub_simu$rr_e[1]!=1) {bias_cc <- c(bias_cc,(mean(sub_simu$beta_cc,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_cc <- c(bias_cc,mean(sub_simu$beta_cc,na.rm=T))
    }
  }
}

hist(bias_nm)
mean_bias <- c(mean(bias_nm),mean(bias_mcim),mean(bias_mimp),mean(bias_cc))
median_bias <- c(median(bias_nm),median(bias_mcim),median(bias_mimp),median(bias_cc))
tail5_bias <- c(length(which(abs(bias_nm)>0.05))/length(bias_nm),
                 length(which(abs(bias_mcim)>0.05))/length(bias_mcim), 
                 length(which(abs(bias_mimp)>0.05))/length(bias_mimp),
                length(which(abs(bias_cc)>0.05))/length(bias_cc))
tail10_bias <- c(length(which(abs(bias_nm)>0.1))/length(bias_nm),
                  length(which(abs(bias_mcim)>0.1))/length(bias_mcim), 
                 length(which(abs(bias_mimp)>0.1))/length(bias_mimp),
                  length(which(abs(bias_cc)>0.1))/length(bias_cc))
beta_table <- cbind(mean_bias,median_bias,tail5_bias,tail10_bias)
beta_table


var_pro_prop <- var_pro[,1]/var_pro[,-1]
var_table <- cbind(apply(var_pro_prop,2,mean),
                   apply(var_pro_prop,2,median))
var_table
var_pro_prop2 <- var_pro[,-1]/var_pro[,1]
var_table2 <- cbind(apply(var_pro_prop2,2,mean),
                    apply(var_pro_prop,2,median))
var_table2

#Case 2, bias=(beta-log(rr_e))/log(rr_e) if rr_e!=1, bias=beta otherwise, keep missing_cutoff < 1000
missing_cutoff <- 1000

missing_rate <- rep(0,4)
bias_nm <- NULL
bias_mcim <- NULL
bias_mimp <- NULL
bias_cc <- NULL
var_nm <- NULL

i=156
for(i in 1:nscena){
  index_i <- (i-1)*nsim+1:nsim
  sub_simu <- read.csv(paste0("./",rdslist[[i]]))
  missing_index_nm <- which(is.na(sub_simu$beta_nm))
  missing_index_mcim <-  which(is.na(sub_simu$beta_mcim))
  missing_index_mimp <-  which(is.na(sub_simu$beta_mimp))
  missing_index_cc <-  which(is.na(sub_simu$beta_cc))
  missing_rate <- missing_rate + c(length(which(is.na(sub_simu$beta_nm))),length(which(is.na(sub_simu$beta_mcim))),
                                   length(which(is.na(sub_simu$beta_mimp))),length(which(is.na(sub_simu$beta_cc))))
  
  if(length(missing_index_nm)<=missing_cutoff){
    if(sub_simu$rr_e[1]!=1) {bias_nm <- c(bias_nm,(mean(sub_simu$beta_nm,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_nm <- c(bias_nm,mean(sub_simu$beta_nm,na.rm=T))
    }
  }
  if(length(missing_index_mcim)<=missing_cutoff){
    if(sub_simu$rr_e[1]!=1) {bias_mcim <- c(bias_mcim,(mean(sub_simu$beta_mcim,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_mcim <- c(bias_mcim,mean(sub_simu$beta_mcim,na.rm=T))
    }
  }
  if(length(missing_index_mimp)<=missing_cutoff){
    if(sub_simu$rr_e[1]!=1) {bias_mimp <- c(bias_mimp,(mean(sub_simu$beta_mimp,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_mimp <- c(bias_mimp,mean(sub_simu$beta_mimp,na.rm=T))
    }
  }
  if(length(missing_index_cc)<=missing_cutoff){
    if(sub_simu$rr_e[1]!=1) {bias_cc <- c(bias_cc,(mean(sub_simu$beta_cc,na.rm=T)-log(sub_simu$rr_e[1]))/log(sub_simu$rr_e[1]))
    }else{ bias_cc <- c(bias_cc,mean(sub_simu$beta_cc,na.rm=T))
    }
  }
}

hist(bias_nm)
mean_bias <- c(mean(bias_nm),mean(bias_mcim),mean(bias_mimp),mean(bias_cc))
median_bias <- c(median(bias_nm),median(bias_mcim),median(bias_mimp),median(bias_cc))
tail5_bias <- c(length(which(abs(bias_nm)>0.05))/length(bias_nm),
                length(which(abs(bias_mcim)>0.05))/length(bias_mcim), 
                length(which(abs(bias_mimp)>0.05))/length(bias_mimp),
                length(which(abs(bias_cc)>0.05))/length(bias_cc))
tail10_bias <- c(length(which(abs(bias_nm)>0.1))/length(bias_nm),
                 length(which(abs(bias_mcim)>0.1))/length(bias_mcim), 
                 length(which(abs(bias_mimp)>0.1))/length(bias_mimp),
                 length(which(abs(bias_cc)>0.1))/length(bias_cc))
beta_table <- cbind(mean_bias,median_bias,tail5_bias,tail10_bias)
beta_table

var_pro3 <- var_pro2[,1]/var_pro2[,-1]
var_table2 <- cbind(apply(var_pro3,2,mean),
                    apply(var_pro,2,median))
var_table2
var_pro4 <- var_pro2[,-1]/var_pro2[,1]
var_table3 <- cbind(apply(var_pro4,2,mean),
                    apply(var_pro4,2,median))
var_table3

beta_table_all <- cbind(beta_table,beta_table2)
rownames(beta_table_all) <- c("NM","MCIM","MI","CC")
beta_table_all

var_table_all <- cbind(var_table,var_table2) 
rownames(var_table_all) <- c("MCIM","MI","CC")
var_table_all
range(var_pro[,3])

apply(var_pro,2,median)
summary(var_pro[,1])
summary(var_pro[,2])
summary(var_pro[,3])
summary(var_pro[,4])
boxplot(var_pro[,1])
which(var_pro[,1]>100000)
subvar_pro <- var_pro[-which(var_pro[,1]>100000),]
apply(subvar_pro,2,mean)
boxplot(subvar_pro[,1])

summary(subvar_pro[,1])
summary(subvar_pro[,2])
summary(subvar_pro[,3])
summary(subvar_pro[,4])
length(which(var_pro2[,1]>100000))







hist(mean_capture_sub[,1],breaks=40)
par(mfrow=c(2,2))
plot(sample_size,mean_capture_sub[,1],ylab="95% coverage probability",main="No missing")
lines(x=seq(0,2000,by=20),y=rep(0.99,101),col="red")
plot(sample_size,mean_capture_sub[,2],ylab="95% coverage probability",main="MCIM")
lines(x=seq(0,2000,by=20),y=rep(0.99,101),col="red")
plot(sample_size,mean_capture_sub[,3],ylab="95% coverage probability",main="MI")
lines(x=seq(0,2000,by=20),y=rep(0.99,101),col="red")
plot(sample_size,mean_capture_sub[,4],ylab="95% coverage probability",main="CC")
lines(x=seq(0,2000,by=20),y=rep(0.99,101),col="red")
sub_simu[1:20,]
dev.off()


A <- which(is.na(sub_simu$beta_nm))
B <- which(is.na(sub_simu$beta_mcim))
C <- which(is.na(sub_simu$beta_mi))
D <- which(is.na(sub_simu$beta_cc))
all(A==B)
all(A==C)
all(A==D)
