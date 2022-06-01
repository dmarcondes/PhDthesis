######################################################################################
#####A data-driven systematic, consistent and feasible approach to Model Selection####
#####4.3 Multilayer W-operator                                                    ####
#####PhD Thesis, Diego Marcondes                                                  ####
#####Universisty of São Paulo, 2022                                               ####
######################################################################################

#Digits classification MNIST
library(keras)
library(tidyverse)
library(ggplot)
library(ggpubr)
library(xtable)

#Read MNIST from keras
c(c(xtrain, ytrain), c(xtest, ytest)) %<-% dataset_mnist()
for(i in 1:nrow(xtrain))
  xtrain[i,,] <- as.numeric(xtrain[i,,] > 0)
for(i in 1:nrow(xtest))
  xtest[i,,] <- as.numeric(xtest[i,,] > 0)

#Divide train into train and validation
set.seed(759)
s <- sample(x = 1:nrow(xtrain),size = floor(0.75*nrow(xtrain)),replace = F)
yval <- ytrain[!(1:nrow(xtrain) %in% s)]
ytrain <- ytrain[(1:nrow(xtrain) %in% s)]
xval <- xtrain[!(1:nrow(xtrain) %in% s),,]
xtrain <- xtrain[(1:nrow(xtrain) %in% s),,]

#Save each image in a file (already saved in MNIST folder)
#for(i in 1:nrow(xtrain))
#  write.table(x = xtrain[i,,],file = paste("train",formatC(x = i,digits = 5,flag = "0"),".txt",sep = ""),sep = " ",row.names = F,col.names = F)
#write.table(x = cbind(ytrain),file = "ytrain.txt",sep = " ",row.names = F,col.names = F)

#for(i in 1:nrow(xval))
#  write.table(x = xval[i,,],file = paste("val",formatC(x = i,digits = 5,flag = "0"),".txt",sep = ""),sep = " ",row.names = F,col.names = F)
#write.table(x = cbind(yval),file = "yval.txt",sep = " ",row.names = F,col.names = F)


#####Process result to calculate validation error####

#Read windows and operators (joint)
W <- read.table("W_Estimated.txt")
joint <- read.table("Joint_Estimated.txt",colClasses = c("factor","factor","factor"))

#Calculate prediction for each image in validation sample
res <- vector()
for(s in 1:nrow(xval)){
  if(s %% 500 == 0)
    print(s)
  M <- rbind(rep(0,29),cbind(rep(0,28),xval[s,,]))
  for(l in 0:6){
    Wtmp <- W %>% filter(V1 == l)
    jointtmp <- joint %>% filter(V1 == l)
    M1 <- matrix(nrow = nrow(M) - 4,ncol = ncol(M) - 4)
    for(i in 3:(nrow(M) - 2))
      for(j in 3:(ncol(M) - 2)){
        val <- ""
        for(k in 1:nrow(Wtmp))
          val <- paste(val,M[i + Wtmp$V2[k],j + Wtmp$V3[k]],sep = "")
        M1[i-2,j-2] <- as.numeric(as.character(jointtmp$V3[jointtmp$V2 == val]))
      }
    M <- M1
  }
  res[s] <- M[1,1]
}

#Save prediction on validation sample
write_rds(x = res,file = "prediction_val.rds")

#Read prediction on validation sample
res <- read_rds("prediction_val.rds")

#Calculate confusion matrix
c <- table(yval == 0,res)

#Validation error
1-sum(diag(prop.table(c)))

#####Test error####

#Read windows and operators (joint)
W <- read.table("W_Estimated.txt")
joint <- read.table("Joint_Estimated.txt",colClasses = c("factor","factor","factor"))

#Calculate prediction for each image in test sample
res <- vector()
for(s in 1:nrow(xtest)){
  if(s %% 500 == 0)
    print(s)
  M <- rbind(rep(0,29),cbind(rep(0,28),xtest[s,,]))
  for(l in 0:6){
    Wtmp <- W %>% filter(V1 == l)
    jointtmp <- joint %>% filter(V1 == l)
    M1 <- matrix(nrow = nrow(M) - 4,ncol = ncol(M) - 4)
    for(i in 3:(nrow(M) - 2))
      for(j in 3:(ncol(M) - 2)){
        test <- ""
        for(k in 1:nrow(Wtmp))
          test <- paste(test,M[i + Wtmp$V2[k],j + Wtmp$V3[k]],sep = "")
        M1[i-2,j-2] <- as.numeric(as.character(jointtmp$V3[jointtmp$V2 == test]))
      }
    M <- M1
  }
  res[s] <- M[1,1]
}

#Save prediction on test sample
write_rds(x = res,file = "prediction_test.rds")

#Read prediction on validation sample
res <- read_rds("prediction_test.rds")

#Calculate confusion matrix
c <- table(ytest == 0,res)
xtable(c)

#Validation error
1-sum(diag(prop.table(c)))

#####Process windows####

#Open estimated windows
W <- read.table("W_Estimated.txt")

#Plot estimated windows
p <- list()
lab <- c(expression(W[1]),expression(W[2]),expression(W[3]),expression(W[4]),expression(W[5]),expression(W[6]),expression(W[7]))

for(l in 0:6){
  w <- W[W$V1 == l,]
  p[[l+1]] <- ggplot(w, aes(V2, V3)) + theme_void() +
    geom_tile(fill = "black") + geom_tile(data = expand.grid(-2:2,-2:2),mapping = aes(x = Var1,y = Var2),inherit.aes = F,color = "black",
                                           fill = NA) + ggtitle(lab[l+1]) + 
    theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
}

#Save plot
pl <- ggarrange(plotlist = p[-7])
tmp <-  ggarrange(plotlist = list(ggplot() + theme_void(),p[[7]],ggplot() + theme_void()),nrow = 1)
pl <- ggarrange(plotlist = list(pl,tmp),nrow = 2,heights = c(2,1))
pdf(file = "windows.pdf",width = 10,height = 10)
pl
dev.off()
  
