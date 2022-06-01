#' @import tidyverse
#' @import ggplot2
#' @export
#' @title Predict output from multiple models
#'
#' @description Predict the output for each observation in a sample and calculate the error for multiple models.
#'
#' @details Given a list of functions from models estimated via search_LS function, compare them by calculating errors and ploting the results.
#'
#' @param sample A data frame with the sample data. Should have columns named x and y representing the input and output strings, and columns variation and date containing the daily variation and date for simulation.
#' @param predict A list of functions which, given an input string, returns a predicted output.
#' @param names A vector containing the names of the models which generated the predict functions.
#' @param simulate Logical. Wheter to simulate the data and plot it.
#' @param init_value Initial value of the account for simulation.
#'
#' @return A list containing the sample data frame with new columns for prediction and simulation, the error of each model, a plot with the simularion and the data used in the plot.
#'
#' @examples
#' mod1 <- search_LS(train = bitcoin[bitcoin$sample == "Train",],
#' val = bitcoin[bitcoin$sample == "Validation",],
#' k = 30,max_leaves = 8,max_init = log(8,2),verbose = FALSE)
#'
#' mod2 <- search_LS(train = bitcoin[bitcoin$sample == "Train",],
#' val = bitcoin[bitcoin$sample == "Validation",],
#' k = 30,max_leaves = 16,max_init = log(16,2),verbose = FALSE)
#'
#' mod3 <- search_LS(train = bitcoin[bitcoin$sample == "Train",],
#' val = bitcoin[bitcoin$sample == "Validation",],
#' k = 30,max_leaves = 32,max_init = log(32,2),verbose = FALSE)
#'
#' e <- error_models(sample = bitcoin[bitcoin$sample == "Test",],
#' predict = list(mod1$predict,mod2$predict,mod3$predict),names = c("d8","d16","d32"),
#' simulate = TRUE)
#'
#' @references Diego Marcondes. A data-driven systematic, consistent and feasible approach to Model Selection. \emph{PhD thesis}. Universidade de São Paulo, 2022

error_models <- function(sample,predict,names,simulate = F,init_value = 1000){
  #Tema para os gráficos
  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 12,
                                                                                 color = "black"), axis.title = element_text(size = 14), legend.text = element_text(size = 14),
                  legend.title = element_text(size = 14), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  legend.position="bottom",legend.spacing.x = unit(0.5, 'cm'))

  result <- NULL
  for(i in 1:length(predict)){
    sample[[gsub(pattern = " ",replacement = "",x = names[i])]] <- NA
    for(j in 1:nrow(sample))
      sample[[gsub(pattern = " ",replacement = "",x = names[i])]][j] <- predict[[i]](sample$x[j])

    #Error
    e <- sum(sample$y != sample[[gsub(pattern = " ",replacement = "",x = names[i])]])/nrow(sample)
    if(is.null(result))
      result <- data.frame("Model" = names[i],"Error" = e)
    else
      result <- rbind.data.frame(result,data.frame("Model" = names[i],"Error" = e))
  }

  #Simulate
  if(simulate){
    for(i in 1:length(predict)){
      sample[[paste("sim_",gsub(pattern = " ",replacement = "",x = names[i]),sep = "")]][1] <- init_value
      for(j in 1:nrow(sample)){
        sample[[paste("sim_",gsub(pattern = " ",replacement = "",x = names[i]),sep = "")]][j] <-
          sample[[paste("sim_",gsub(pattern = " ",replacement = "",x = names[i]),sep = "")]][j]*(1 +
                                            as.numeric(sample[[gsub(pattern = " ",replacement = "",x = names[i])]][j])*sample$variation[j])
        if(j < nrow(sample)){
          sample[[paste("sim_",gsub(pattern = " ",replacement = "",x = names[i]),sep = "")]][j+1] <-
            sample[[paste("sim_",gsub(pattern = " ",replacement = "",x = names[i]),sep = "")]][j]
        }
      }
    }

    #Obseverd
    sample$sim_Observed[1] <- init_value
    for(j in 1:nrow(sample)){
      sample$sim_Observed[j] <- sample$sim_Observed[j]*(1 + sample$variation[j])
      if(j < nrow(sample)){
        sample$sim_Observed[j+1] <- sample$sim_Observed[j]
      }
    }

    #Plot
    tmp <- sample %>% select(date,paste("sim_",gsub(pattern = " ",replacement = "",x = names),sep = ""),"sim_Observed") %>% gather("Model","value",-date)
    tmp$Model <- factor(gsub(pattern = "sim_",replacement = "",x = tmp$Model),c("Observed",names))
    p <- ggplot(tmp,aes(x = date,y = value,group = Model,color = Model)) + theme_linedraw() + titles + geom_line() + xlab("Date") + ylab("Simulated value")
  }
  else{
    tmp <- NULL
    p <- NULL
  }

  return(list("sample" = sample,"result" = result,"plot" = p,"data_plot" = tmp))
}
