######################################################################################
#####A data-driven systematic, consistent and feasible approach to Model Selection####
#####4.2 Forecasting variable order Markov chains                                 ####
#####PhD Thesis, Diego Marcondes                                                  ####  
#####Universisty of São Paulo, 2022                                               ####
######################################################################################

#Install MarkovLS package
library(devtools)
install_github("dmarcondes/PhDthesis/MarkovLS")

#Packages
library(tidyverse)
library(ggplot2)
library(lubridate)
library(MarkovLS)
library(xtable)
library(RColorBrewer)
library(ggpubr)

titles <- theme(strip.text = element_text(size = 18), axis.text = element_text(size = 14,
                                                                               color = "black"),
                axis.title = element_text(size = 18), legend.text = element_text(size = 18),
                legend.title = element_text(size = 20), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), panel.border = element_blank(),
                panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                legend.position="bottom",legend.spacing.x = unit(0.5, 'cm'))

#####Train models####
set.seed(10)

#Read data
dados <- readRDS("data_bitcoin.rds")
summary(dados)

#Divide data
train <- dados %>% filter(sample == "Train")
val <- dados %>% filter(sample == "Validation")
test <- dados %>% filter(sample == "Test")

#Train models via MarkovLS package
mod_bit8 <- search_LS(train = train,val = val,k = 30,max_leaves = 8,max_init = 3,verbose = T)

mod_bit16 <- search_LS(train = train,val = val,k = 30,max_leaves = 16,max_init = 4,verbose = T)

mod_bit32 <- search_LS(train = train,val = val,k = 30,max_leaves = 32,max_init = 5,verbose = T)

mod_bit64 <- search_LS(train = train,val = val,k = 30,max_leaves = 64,max_init = 6,verbose = T)

mod_bit128 <- search_LS(train = train,val = val,k = 30,max_leaves = 128,max_init = 7,verbose = T)

mod_bit256 <- search_LS(train = train,val = val,k = 30,max_leaves = 256,max_init = 8,verbose = T)

mod_bit512 <- search_LS(train = train,val = val,k = 30,max_leaves = 512,max_init = 9,verbose = T)
l <- list("bit8" = mod_bit8,"bit_16" = mod_bit16,"bit_32" = mod_bit32,"bit_64" = mod_bit64,"bit_128" = mod_bit128,
          "bit_256" = mod_bit256,"bit_512" = mod_bit512)
write_rds(l,"results_bitcoin.rds")

#####Results####
#Data
dados <- readRDS("data_bitcoin.rds")
train <- dados %>% filter(sample == "Train")
val <- dados %>% filter(sample == "Validation")
test <- dados %>% filter(sample == "Test")

#Results
l <- readRDS("results_bitcoin.rds")
mod_bit8 <- l$bit8
mod_bit16 <- l$bit_16
mod_bit32 <- l$bit_32
mod_bit64 <- l$bit_64
mod_bit128 <- l$bit_128
mod_bit256 <- l$bit_256
mod_bit512 <- l$bit_512

xtable(mod_bit8$hhat)

#Predict
pred <- list(mod_bit8$predict,mod_bit16$predict,mod_bit32$predict,mod_bit64$predict,mod_bit128$predict,
             mod_bit256$predict,mod_bit512$predict)
names <- c("LS8","LS16","LS32","LS64","LS128","LS256","LS512")

#Compare models
compare_bit <- error_models(sample = test,predict = pred,names = names,simulate = T,init_value = 1000)

#Table with results
tab <- data.frame("Model" = c("LS8","LS16","LS32","LS64","LS128","LS256","LS512"),
                  "Time" = as.numeric(c(mod_bit8$time,mod_bit16$time,mod_bit32$time,mod_bit64$time,
                                        mod_bit128$time,mod_bit256$time,mod_bit512$time))/60,
                  "Order" = c(max(nchar(mod_bit8$hhat$x)),max(nchar(mod_bit16$hhat$x)),max(nchar(mod_bit32$hhat$x)),max(nchar(mod_bit64$hhat$x)),
                              max(nchar(mod_bit128$hhat$x)),max(nchar(mod_bit256$hhat$x)),max(nchar(mod_bit512$hhat$x))),
                  "VCd" = c(nrow(mod_bit8$hhat),nrow(mod_bit16$hhat),nrow(mod_bit32$hhat),nrow(mod_bit64$hhat),
                            nrow(mod_bit128$hhat),nrow(mod_bit256$hhat),nrow(mod_bit512$hhat)),
                  "ValError" = c(mod_bit8$error,mod_bit16$error,mod_bit32$error,mod_bit64$error,
                                 mod_bit128$error,mod_bit256$error,mod_bit512$error),
                  "TestError" = compare_bit$result$Error,
                  "MinSpread" = c(min(compare_bit$sample$sim_LS8/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS16/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS32/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS64/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS128/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS256/compare_bit$sample$sim_Observed)-1,
                                  min(compare_bit$sample$sim_LS512/compare_bit$sample$sim_Observed)-1),
                  "MaxSpread" = c(max(compare_bit$sample$sim_LS8/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS16/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS32/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS64/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS128/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS256/compare_bit$sample$sim_Observed)-1,
                                  max(compare_bit$sample$sim_LS512/compare_bit$sample$sim_Observed)-1),
                  "FinalSpread" = c((compare_bit$sample$sim_LS8/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS16/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS32/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS64/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS128/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS256/compare_bit$sample$sim_Observed)[nrow(test)]-1,
                                    (compare_bit$sample$sim_LS512/compare_bit$sample$sim_Observed)[nrow(test)]-1),
                  "TrainSize" = nrow(train %>% filter(nchar(x) >= 30)),
                  "ValSize" = nrow(val),
                  "TestSize" = nrow(test))
xtable(tab,digits = 3)
tapply(dados$date,dados$sample,summary)

#Plot
plot <- compare_bit$data_plot
tmp <- plot %>% filter(Model == "Observed") %>% select(date,value)
names(tmp)[2] <- "observed"
plot <- merge(plot %>% filter(Model != "Observed"),tmp,by = "date")
plot$Model <- factor(plot$Model)
plot$Model <- plyr::mapvalues(plot$Model,levels(plot$Model),paste("d =",c(8,16,32,64,128,256,512)))
dados$sample <- plyr::mapvalues(dados$sample,"Train","Training")

p0 <- ggplot(dados,aes(x = date,y = close,color = sample,group = 1)) + theme_linedraw() + titles + geom_line() +
  scale_y_log10(breaks = (2^(0:10))*100,labels = scales::comma) + xlab("Date") +
  ylab("Bitcoin value at end of day in US Dollars") +
  scale_colour_manual("Sample",values = c("dodgerblue2","black","chocolate1"),breaks = c("Training","Validation","Test")) +
  scale_x_date(breaks =ymd("2013-07-01","2014-07-01","2015-07-01","2016-07-01","2017-07-01","2018-07-01","2019-07-01","2020-07-01",
                           "2021-07-01","2022-04-06"),labels = c(paste("Jul-",2013:2021,sep = ""),"Apr-22"))

p1 <- ggplot(plot %>% filter(Model != "d = 512"),aes(x = date,y = value,group = Model)) + theme_linedraw() + titles +
  geom_line(aes(x = date,y = observed,group = 1),color = "red") +
  geom_line(color = "chartreuse3") + facet_wrap(~Model,nrow = 2) +
  geom_hline(yintercept = 1000,linetype = "dashed") + xlab("") + ylab("Account balance") +
  scale_y_continuous(breaks = seq(1000,10000,500),labels = scales::comma) +
  scale_x_date(breaks = ymd(c("2021-02-01","2021-06-01","2021-11-01","2022-03-01")),
               labels = c("Feb-2021","Jun-2021","Nov-2021","Mar-2022"))

p2 <- ggplot(plot %>% filter(Model == "d = 512"),aes(x = date,y = value,group = Model)) + theme_linedraw() + titles +
  geom_line(aes(x = date,y = observed,group = 1,color = "Stay on the market")) +
  geom_line(aes(color = "Predicted strategy")) + facet_wrap(~Model) +
  geom_hline(yintercept = 1000,linetype = "dashed") +xlab("Date") + scale_color_manual("",breaks = c("Stay on the market",
                                                                                                     "Predicted strategy"),
                                                                                       values = c("red","chartreuse3")) +
  ylab("") +
  scale_y_continuous(breaks = seq(1000,10000,500),labels = scales::comma) +
  scale_x_date(breaks = ymd(c("2021-02-01","2021-06-01","2021-11-01","2022-03-01")),
               labels = c("Feb-2021","Jun-2021","Nov-2021","Mar-2022"))
pdf(file = "example_bitcoin.pdf",width = 15,height = 15)
ggarrange(ggarrange(p0,labels = "(A)"),ggarrange(p1,labels = "(B)"),ggarrange(ggplot() + theme_void(),p2,ggplot() + theme_void(),nrow = 1),
          nrow=3,heights = c(7.5,7.5,7.5/2))
dev.off()
