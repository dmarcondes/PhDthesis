######################################################################################
#####A data-driven systematic, consistent and feasible approach to Model Selection####
#####4.1 Learning via the Partition Lattice Learning Space                        ####
#####PhD Thesis, Diego Marcondes                                                  ####
#####Universisty of São Paulo, 2022                                               ####
######################################################################################

#Install U-curve package
#library(devtools)
#install_github() #PREENCHER

#Libraries
library(partitionUcurve)
library(tidyverse)
library(xtable)

#Source utils.R from partitionUcurve package
source("utils.R")

#####Auxiliary functions####
#Get test error
minError <- function(h,xtest,ytest){
  tab <- prop.table(table(xtest,ytest))
  e <- 0
  for(i in 1:nrow(tab))
    e <- e + tab[i,2-h$fx[i]]
  return(e)
}

#Get real error from joint distribution
realError <- function(h,joint,prob){
  j <- data.frame(data.table::rbindlist(joint))
  j$p <- prob
  h$X <- levels(j$Var2)[as.numeric(h$X)]
  e <- sum(apply(h,1,function(x) j$p[j$Var2 == x[1] & j$Var1 != x[2]]))
  return(e)
}

#Simulate optimal search on the PLLS for given joint distribution
simulate_PLLS <- function(joint,prob,n,m,stop,exhaust = 100,rep,cores,verbose = F,path,Xd,Yd){
  #joint: list with joint values of X and Y
  #prob: vector of probabilities of same order as joint
  #n: vector of sample sizes to consider
  #m: number of points in X domain
  #stop: number of models to exhaust before stopping optimal algorithm and performing exhaustive search
  #exhaust: number of models to exhaust before stopping suboptimal algorithm
  #rep: number of repetitions for each sample size
  #cores: number of cores to use
  #verbose: trace ucurve algorithm
  #path: prefix of files with simulation results
  #X: vector with domain variable x
  #Y: vector with domain variable y

  #Read lattice
  Lh <- fst::read.fst(paste("part_",m,".fst",sep = ""))

  #Tab to store solutions
  tab <- data.frame("Optimal" = NA,"n" = NA,"exhausted" = NA,"remain" = NA,"number solutions" = NA,"val error" = NA,"real error" = NA,"ERM error" = NA,"Mhat error" = NA,
                    "time" = NA)
  trace <- NULL
  info <- NULL

  #Joint as tab
  tab_joint <- data.table::rbindlist(joint) %>% data.frame()
  tab_joint$prob <- prob
  tab_joint <- tab_joint %>% spread("Var1","prob")
  tab_joint <- cbind(tab_joint[,2],tab_joint[,3])

  #For each sample size
  for(size in n){
    cat("\n")
    cat(paste("Sample size =",size))
    cat("\n")

    #Train, validation and independent sample size
    ntrain <- floor(0.5*size)
    nval <- floor(0.25*size)
    nindependent <- size - ntrain - nval

    #Vectors to store results optimal search
    nsolutions <- vector()
    val_error <- vector()
    real_error <- vector()
    erm_error <- vector()
    Mhat_error <- vector()
    remain <- vector()
    exhausted <- vector()
    time <- vector()

    #Vectors to store results suboptimal search
    est_nsolutions <- vector()
    est_val_error <- vector()
    est_real_error <- vector()
    est_erm_error <- vector()
    est_Mhat_error <- vector()
    est_remain <- vector()
    est_exhausted <- vector()
    est_time <- vector()

    #Repeat
    for(i in 1:rep){
      cat("\n")
      cat(i)
      cat("\n")

      #Sample train data
      train <- data.frame(data.table::rbindlist(sample(x = joint,size = ntrain,replace = T,prob = prob)))
      names(train) <- c("y","x")

      #Sample validation data
      val <- data.frame(data.table::rbindlist(sample(x = joint,size = nval,replace = T,prob = prob)))
      names(val) <- c("y","x")

      #Sample independent data
      independent <- data.frame(data.table::rbindlist(sample(x = joint,size = nindependent,replace = T,prob = prob)))
      names(independent) <- c("y","x")

      #Domains
      independent$x <- factor(independent$x,Xd)
      independent$y <- factor(independent$y,Yd)

      #Run U-curve Optimal
      u <- ucurve(xtrain = train$x,ytrain = train$y,xval = val$x,yval = val$y,Lh = Lh,cores = cores,stop = stop,verbose = verbose,X = Xd,
                  Y = Yd)

      #Get estimated hypotheses of each returned partition
      est_h <- unique(lapply(rbind(u$partitions),function(x) optimalHyp(part = x,jtrain = prop.table(table(independent$x,independent$y)),
                                                                        increasing = F,n = length(unique(independent$x)))))

      #Get target hypotheses of returned partitions
      hstar <- unique(lapply(rbind(u$partitions),function(x) optimalHyp(part = x,jtrain = tab_joint,
                                                                        increasing = F,n = length(unique(independent$x)))))

      #Get ERM hypotheses of learning directly in H
      erm_hyp <- optimalHyp(part = namePartition(makePartition(as.list(1:m))),
                            jtrain = prop.table(table(c(train$x,val$x,independent$x),c(train$y,val$y,as.character(independent$y)))),increasing = F,
                            n = length(unique(c(train$x,val$x,independent$x))))

      #Store results
      exhausted[i] <- u$exhausted #Number of models exhausted
      remain[i] <- u$remain #Models remaining when
      nsolutions[i] <- length(unique(est_h)) #Number of solutions
      val_error[i] <- u$error #Validation error
      real_error[i] <- min(unlist(lapply(est_h,function(x) realError(x,joint,prob)))) #Real error of best estimated hypotheses
      Mhat_error[i] <- min(unlist(lapply(hstar,function(x) realError(x,joint,prob)))) #Real error of Mhat
      time[i] <- u$time #Running time
      erm_error[i] <- realError(erm_hyp,joint,prob) #Real error of best erm hypotheses

      #Store trace of algorithm
      if(is.null(trace))
        trace <- data.frame("optimal" = "Yes","n" = size,"id" = i,"ex" = u$exhausted_until_prune,"rem" = u$remain_after_prune)
      else
        trace <- rbind.data.frame(trace,data.frame("optimal" = "Yes","n" = size,"id" = i,"ex" = u$exhausted_until_prune,"rem" = u$remain_after_prune))

      #Run Suboptimal U-curve un
      u <- ucurve(xtrain = train$x,ytrain = train$y,xval = val$x,yval = val$y,Lh = Lh,cores = cores,stop = stop,verbose = verbose,X = Xd,
                  Y = Yd,optimal = F,exhaust = exhaust)

      #Get estimated hypotheses of each returned partition
      est_h <- unique(lapply(rbind(u$partitions),function(x) optimalHyp(part = x,jtrain = prop.table(table(independent$x,independent$y)),
                                                                        increasing = F,n = length(unique(independent$x)))))

      #Get target hypotheses of returned partitions
      hstar <- unique(lapply(rbind(u$partitions),function(x) optimalHyp(part = x,jtrain = tab_joint,
                                                                        increasing = F,n = length(unique(independent$x)))))

      #Get ERM hypotheses of learning directly in H
      erm_hyp <- optimalHyp(part = namePartition(makePartition(as.list(1:m))),
                            jtrain = prop.table(table(c(train$x,val$x,independent$x),c(train$y,val$y,as.character(independent$y)))),increasing = F,
                            n = length(unique(c(train$x,val$x,independent$x))))

      #Store results
      est_exhausted[i] <- u$exhausted #Number of models exhausted
      est_remain[i] <- u$remain #Models remaining when
      est_nsolutions[i] <- length(unique(est_h)) #Number of solutions
      est_val_error[i] <- u$error #Validation error
      est_real_error[i] <- min(unlist(lapply(est_h,function(x) realError(x,joint,prob)))) #Real error of best estimated hypotheses
      est_Mhat_error[i] <- min(unlist(lapply(hstar,function(x) realError(x,joint,prob)))) #Real error of Mhat
      est_time[i] <- u$time #Running time
      est_erm_error[i] <- realError(erm_hyp,joint,prob) #Real error of best erm hypotheses

      #Store trace of algorithm
      trace <- rbind.data.frame(trace,data.frame("optimal" = "No","n" = size,"id" = i,"ex" = u$exhausted_until_prune,"rem" = u$remain_after_prune))
      rm(train,val,independent,u,est_h)
    }
    #Summarize optimal search with given sample size
    tab <- rbind.data.frame(tab,data.frame("Optimal" = "Yes","n" = size,
                                           "exhausted" = paste(round(median(exhausted),2)," (",
                                                               round(quantile(exhausted,0.025),2),",",
                                                               round(quantile(exhausted,0.975),2),")",sep = ""),
                                           "remain" = paste(round(median(remain),2)," (",
                                                            round(quantile(remain,0.025),2),",",
                                                            round(quantile(remain,0.975),2),")",sep = ""),
                                           "number solutions" = paste(median(nsolutions)," (",
                                                                      quantile(nsolutions,0.025),",",
                                                                      quantile(nsolutions,0.975),")",sep = ""),
                                           "val error" = paste(round(median(val_error),5)," (",
                                                               round(quantile(val_error,0.025),5),",",
                                                               round(quantile(val_error,0.975),5),")",sep = ""),
                                           "real error" = paste(round(median(real_error),5)," (",
                                                                round(quantile(real_error,0.025),5),",",
                                                                round(quantile(real_error,0.975),5),")",sep = ""),
                                           "ERM error" = paste(round(median(erm_error),5)," (",
                                                               round(quantile(erm_error,0.025),5),",",
                                                               round(quantile(erm_error,0.975),5),")",sep = ""),
                                           "Mhat error" = paste(round(median(Mhat_error),5)," (",
                                                                round(quantile(Mhat_error,0.025),5),",",
                                                                round(quantile(Mhat_error,0.975),5),")",sep = ""),
                                           "time" = paste(round(median(time),5)," (",
                                                          round(quantile(time,0.025),5),",",
                                                          round(quantile(time,0.975),5),")",sep = "")))

    #Summarize suboptimal search with given sample size
    tab <- rbind.data.frame(tab,data.frame("Optimal" = "No","n" = size,
                                           "exhausted" = paste(round(median(est_exhausted),2)," (",
                                                               round(quantile(est_exhausted,0.025),2),",",
                                                               round(quantile(est_exhausted,0.975),2),")",sep = ""),
                                           "remain" = NA,
                                           "number solutions" = paste(median(est_nsolutions)," (",
                                                                      quantile(est_nsolutions,0.025),",",
                                                                      quantile(est_nsolutions,0.975),")",sep = ""),
                                           "val error" = paste(round(median(est_val_error),5)," (",
                                                               round(quantile(est_val_error,0.025),5),",",
                                                               round(quantile(est_val_error,0.975),5),")",sep = ""),
                                           "real error" = paste(round(median(est_real_error),5)," (",
                                                                round(quantile(est_real_error,0.025),5),",",
                                                                round(quantile(est_real_error,0.975),5),")",sep = ""),
                                           "ERM error" = paste(round(median(est_erm_error),5)," (",
                                                               round(quantile(est_erm_error,0.025),5),",",
                                                               round(quantile(est_erm_error,0.975),5),")",sep = ""),
                                           "Mhat error" = paste(round(median(est_Mhat_error),5)," (",
                                                                round(quantile(est_Mhat_error,0.025),5),",",
                                                                round(quantile(est_Mhat_error,0.975),5),")",sep = ""),
                                           "time" = paste(round(median(est_time),5)," (",
                                                          round(quantile(est_time,0.025),5),",",
                                                          round(quantile(est_time,0.975),5),")",sep = "")))

    #Store all information about simulations
    if(is.null(info))
      info <- data.frame("id"= 1:rep,"Optimal" = "Yes","n" = size,"exhausted" = exhausted,"remain" = remain,"number solutions" = nsolutions,
                         "val error" = val_error,"real error" = real_error,"ERM error" = erm_error,
                         "Mhat error" = Mhat_error,
                         "time" = time)
    else
      info <- rbind.data.frame(info,data.frame("id"= 1:rep,"Optimal" = "Yes","n" = size,"exhausted" = exhausted,"remain" = remain,"number solutions" = nsolutions,
                                               "val error" = val_error,"real error" = real_error,"ERM error" = erm_error,
                                               "Mhat error" = Mhat_error,
                                               "time" = time))
    info <- rbind.data.frame(info,data.frame("id"= 1:rep,"Optimal" = "No","n" = size,"exhausted" = est_exhausted,"remain" = est_remain,"number solutions" = est_nsolutions,
                                             "val error" = est_val_error,"real error" = est_real_error,"ERM error" = est_erm_error,
                                             "Mhat error" = est_Mhat_error,
                                             "time" = est_time))

    #Save information stored so far
    write.table(x = tab,file = paste(path,".txt",sep = ""))
    write_rds(x = list("table" = tab,"trace" = trace,"info" = info),file = paste(path,".rds",sep = ""))
  }
  #Return:
  #tabe: table with the summary of simulations
  #trace: the trace of all U-curve calls
  #info: information about all simulated samples
  tab <- tab[-1,]
  return(list("tabe" = tab,"trace" = trace,"info" = info))
}

#####Examples####

#Example 1
joint <- split(expand.grid(0:1,c("01","02","03","04","05","06","07","08")),seq(16))
prob <- c(0.05,0.95,
          0.02,0.98,
          0.51,0.49,
          0.51,0.49,
          0.51,0.49,
          0.90,0.10,
          0.8,0.2,
          0.99,0.01)*(1/8)
set.seed(93893475)
ex1 <- simulate_PLLS(joint = joint,prob = prob,n = c(64,96,128,160,192,224,256),m = 8,stop = 1000,rep = 100,cores = 6,
                     path = "example1",verbose = F,
                     Xd = c("01","02","03","04","05","06","07","08"),Yd = c("0","1"))

#Example 2
joint <- split(expand.grid(0:1,c("01","02","03","04","05","06","07","08")),seq(16))
prob <- c(0.05,0.95,
          0.02,0.98,
          0.52,0.48,
          0.52,0.48,
          0.52,0.48,
          0.90,0.10,
          0.8,0.2,
          0.96,0.04)*(1/8)
set.seed(645645)
ex2 <- simulate_PLLS(joint = joint,prob = prob,n = c(64,96,128,160,192,224,256),m = 8,stop = 1000,rep = 100,cores = 6,path = "example2",verbose = F,
                        X = c("01","02","03","04","05","06","07","08"),Y = c("0","1"))

#Example 3
joint <- split(expand.grid(0:1,c("01","02","03","04","05","06","07","08")),seq(16))
prob <- c(0.05,0.95,
          0.02,0.98,
          0.56,0.44,
          0.56,0.44,
          0.56,0.44,
          0.90,0.10,
          0.8,0.2,
          0.84,0.16)*(1/8)
set.seed(453457)
ex3 <- simulate_PLLS(joint = joint,prob = prob,n = c(64,96,128,160,192,224,256),m = 8,stop = 1000,rep = 100,cores = 6,path = "example3",verbose = F,
                        X = c("01","02","03","04","05","06","07","08"),Y = c("0","1"))

#Example 4
joint <- split(expand.grid(0:1,c("01","02","03","04","05","06","07","08")),seq(16))
prob <- c(0.05,0.95,
          0.02,0.98,
          0.7,0.3,
          0.7,0.3,
          0.7,0.3,
          0.90,0.10,
          0.52,0.48,
          0.7,0.3)*(1/8)
set.seed(64567)
ex4 <- simulate_PLLS(joint = joint,prob = prob,n = c(64,96,128,160,192,224,256),m = 8,stop = 1000,rep = 100,cores = 6,path = "example4",verbose = F,
                        X = c("01","02","03","04","05","06","07","08"),Y = c("0","1"))

#####Results of simulations####

#Theme for ggplot
titles <- theme(strip.text = element_text(size = 20), axis.text = element_text(size = 14,
                                                                               color = "black"),
                axis.title = element_text(size = 20), legend.text = element_text(size = 20),
                legend.title = element_text(size = 14), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), panel.border = element_blank(),
                panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                legend.position="bottom",legend.spacing.x = unit(0.5, 'cm'))

#Oprn RDS results
results <- list()
for(i in 1:4)
  results[[i]] <- readRDS(paste("example",i,".rds",sep = ""))

#Table of results
tab <- data.frame("Example" = 1,results[[1]]$table)
for(i in 2:4)
  tab <- rbind.data.frame(tab,data.frame("Example" = i,results[[i]]$table))

#Print table in latex
print(xtable(tab[,-c(5,6)]),include.rownames = F)

#Comparing via PLLS x via erm
info <- data.frame("ex" = 1,results[[1]]$info)
for(i in 2:4)
  info <- rbind.data.frame(info,data.frame("ex" = i,results[[i]]$info))
View(info)

tmp <- info %>% select(ERM.error,real.error,Optimal,n,id,ex)
tmp <- tmp %>% group_by(Optimal,n,ex) %>% summarise("Learning via Learning Spaces better than learning via ERM" = 100*sum(real.error < ERM.error)/length(real.error),
                                                    "Learning via Learning Spaces as good as learning via ERM" = 100*sum(real.error == ERM.error)/length(real.error),
                                                    "Learning via Learning Spaces worse than learning via ERM" = 100*sum(real.error > ERM.error)/length(real.error)) %>%
  data.frame()
tmp <- tmp %>% gather("cat","porc",-Optimal,-n,-ex)
tmp$cat <- factor(tmp$cat,c("Learning.via.Learning.Spaces.worse.than.learning.via.ERM",
                            "Learning.via.Learning.Spaces.as.good.as.learning.via.ERM",
                            "Learning.via.Learning.Spaces.better.than.learning.via.ERM"))
tmp$cat <- plyr::mapvalues(tmp$cat,levels(tmp$cat),c("Learning via Learning Spaces worse than learning via ERM",
                                                     "Learning via Learning Spaces as good as learning via ERM",
                                                     "Learning via Learning Spaces better than learning via ERM"))
tmp$n <- factor(tmp$n)
tmp$Optimal <- factor(tmp$Optimal)
tmp$Optimal <- plyr::mapvalues(tmp$Optimal,c("Yes","No"),c("Optimal","Suboptimal"))
tmp$ex <- factor(tmp$ex)
tmp$ex <- plyr::mapvalues(tmp$ex,1:4,c("Example 1","Example 2","Example 3","Example 4"))

p <- ggplot(tmp,aes(x = n,y = porc,fill = cat)) + theme_linedraw() + titles + facet_grid(ex~Optimal) + geom_bar(stat = "identity") +
  scale_fill_manual("",values = c("red","gold","green")) + ylab("Percentage of simulated samples (%)") + xlab("Sample size") +
  geom_text(aes(label = paste(round(porc),"%",sep = "")),position = "stack",vjust = 1.2,col = "black",size = 5) +
  guides(fill = guide_legend(nrow = 3,byrow = TRUE))
pdf("SimulationPartition.pdf",width = 0.9*20,height = 0.9*15)
p
dev.off()

#Comparing optimal x Suboptimal algorithm
tmp <- info %>% select(real.error,Optimal,n,id,ex)
tmp$key <- paste(tmp$ex,tmp$id,tmp$n)
tmp1 <- tmp %>% filter(Optimal == "Yes") %>% select(key,real.error,n,id,ex)
colnames(tmp1)[2] <- c("erro.opt")
tmp2 <- tmp %>% filter(Optimal == "No") %>% select(key,real.error)
colnames(tmp2)[2] <- c("erro.sub")
tmp <- merge(tmp1,tmp2)

tmp <- tmp %>% group_by(ex,n) %>% summarise("SuboptimalBetter" = 100*sum(erro.sub < erro.opt)/length(erro.sub),
                                            "Equal" = 100*sum(erro.sub == erro.opt)/length(erro.sub),
                                            "OptimalBetter" = 100*sum(erro.sub > erro.opt)/length(erro.sub),) %>%
  data.frame()
tmp <- data.frame(tmp[1:14,],tmp[15:28,])
print(xtable(tmp),include.rownames = F)

#Joint distributions of simulations
prob1 <- matrix(c(0.05,0.95,
                  0.02,0.98,
                  0.51,0.49,
                  0.51,0.49,
                  0.51,0.49,
                  0.90,0.10,
                  0.8,0.2,
                  0.99,0.01),nrow = 8,ncol = 2,byrow = T)*(1/8)
prob2 <- matrix(c(0.05,0.95,
                  0.02,0.98,
                  0.52,0.48,
                  0.52,0.48,
                  0.52,0.48,
                  0.90,0.10,
                  0.8,0.2,
                  0.96,0.04),nrow = 8,ncol = 2,byrow = T)*(1/8)
prob3 <- matrix(c(0.05,0.95,
                  0.02,0.98,
                  0.56,0.44,
                  0.56,0.44,
                  0.56,0.44,
                  0.90,0.10,
                  0.8,0.2,
                  0.84,0.16),nrow = 8,ncol = 2,byrow = T)*(1/8)
prob4 <- matrix(c(0.05,0.95,
                  0.02,0.98,
                  0.7,0.3,
                  0.7,0.3,
                  0.7,0.3,
                  0.90,0.10,
                  0.52,0.48,
                  0.7,0.3),nrow = 8,ncol = 2,byrow = T)*(1/8)
probs <- cbind(prob1,prob2,prob3,prob4)
rownames(probs) <- 1:8

#Calculate epsilonstar, minimum error and conditional entropy of joint distributions
epsilonStar <- c(min(apply(prob1,1,function(x) abs(x[1] - x[2]))),NA,
                 min(apply(prob2,1,function(x) abs(x[1] - x[2]))),NA,
                 min(apply(prob3,1,function(x) abs(x[1] - x[2]))),NA,
                 min(apply(prob4,1,function(x) abs(x[1] - x[2]))),NA)
entropy <- c(-sum(prob1*log(prob1/(1/8))),NA,
             -sum(prob2*log(prob2/(1/8))),NA,
             -sum(prob3*log(prob3/(1/8))),NA,
             -sum(prob4*log(prob4/(1/8))),NA)
error <- c(sum(apply(prob1,1,min)),NA,sum(apply(prob2,1,min)),NA,sum(apply(prob3,1,min)),NA,sum(apply(prob4,1,min)),NA)
probs <- rbind(probs,epsilonStar,entropy,error)

#Print joint distributions
print(xtable(data.frame(probs),digits = 5),include.rownames = T)

