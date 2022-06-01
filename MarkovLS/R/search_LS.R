#' @import tidyverse
#' @export
#' @title Search the Markov Boolean Learning Space
#'
#' @description Performs an stochastic search on the Markov Boolean Learning Space.
#'
#' @details
#'
#' @param train A data frame with the training data. Should have collumns named x and y representing the input and output strings.
#' @param val A data frame with the validation data. Should have collumns named x and y representing the input and output strings.
#' @param k Maximum length of Markov Chain.
#' @param max_leaves Maximum number of leaves in a node.
#' @param max_init Maximum length of initial node to try.
#' @param verbose Logical. Whether to print algorithm progress.
#' @param error Wheter to consider the simple error or the weighted error
#'
#' @return A list containing the validation error of the estimated node, a table with the estimated hypothesis, a function to predict the output of a given input and the algorithm time.

search_LS <- function(train,val,k,max_leaves = 50,max_init,verbose = T,error = "simple"){
  start_time <- Sys.time()

  #Lost sample
  train <- train %>% filter(nchar(x) >= k)

  #All sample
  all_sample <- rbind.data.frame(train,val)

  #x that matter
  train$x <- substr(train$x,1,k)
  val$x <- substr(val$x,1,k)
  all_sample$x <- substr(all_sample$x,1,k)

  #Frequency tables
  if(error == "simple"){
    tab_train <- as.matrix(table(train$x,train$y))
    tab_val <- as.matrix(table(val$x,val$y))
  }
  else if(error == "weighted"){
    #Train
    train$key <- paste(train$x,train$y,sep = "-")
    tmp <- train %>% group_by(key) %>% summarise("weight" = sum(abs(weight))) %>% as.matrix()
    tmp <- cbind(tmp,unlist(lapply(strsplit(tmp[,1],split = "-"),function(x) x[1])))
    tmp <- cbind(tmp,unlist(lapply(strsplit(tmp[,1],split = "-"),function(x) x[2])))
    tmp <- tmp[,-1]
    tab_train <- matrix(data = 0,nrow = length(unique(tmp[,2])),ncol = 2)
    i <- 1
    for(x in unique(tmp[,2])){
      if("1" %in% tmp[tmp[,2] == x,3])
        tab_train[i,2] <- as.numeric(tmp[tmp[,2] == x & tmp[,3] == "1",1])
      if("0" %in% tmp[tmp[,2] == x,3])
        tab_train[i,1] <- as.numeric(tmp[tmp[,2] == x & tmp[,3] == "0",1])
      i <- i + 1
    }
    rownames(tab_train) <- unique(tmp[,2])

    #Val
    val$key <- paste(val$x,val$y,sep = "-")
    tmp <- val %>% group_by(key) %>% summarise("weight" = sum(abs(weight))) %>% as.matrix()
    tmp <- cbind(tmp,unlist(lapply(strsplit(tmp[,1],split = "-"),function(x) x[1])))
    tmp <- cbind(tmp,unlist(lapply(strsplit(tmp[,1],split = "-"),function(x) x[2])))
    tmp <- tmp[,-1]
    tab_val <- matrix(data = 0,nrow = length(unique(tmp[,2])),ncol = 2)
    i <- 1
    for(x in unique(tmp[,2])){
      if("1" %in% tmp[tmp[,2] == x,3])
        tab_val[i,2] <- as.numeric(tmp[tmp[,2] == x & tmp[,3] == "1",1])
      if("0" %in% tmp[tmp[,2] == x,3])
        tab_val[i,1] <- as.numeric(tmp[tmp[,2] == x & tmp[,3] == "0",1])
      i <- i + 1
    }
    rownames(tab_val) <- unique(tmp[,2])
  }

  #Current error
  current_error <- 1
  current_node <- NULL

  #For each initial value
  for(init in 1:max_init){
    #First node
    node <- clean_node(as.list(apply(expand.grid(replicate(init, 0:1, simplify = FALSE)),1,function(x) paste(x,collapse = ""))),tab_train)
    error <- get_error_node(tab_train,tab_val,node,nrow(val))

    while(length(node) < max_leaves){
      if(verbose){
          cat(paste("Init: ",init," Leaves: ",length(node)," Error: ",round(error,3),
                    " Current error: ",round(min(error,current_error),3),sep = ""))
          cat("\n")
      }
      keep <- NULL
      keep_error <- error
      for(i in 1:length(node)){
        #Add leaves
        new_node <- node[-i]
        new_node <- append(new_node,list(paste(node[[i]],"0",sep = ""),paste(node[[i]],"1",sep = "")))
        new_error <- get_error_node(tab_train,tab_val,new_node,nrow(val))
        if(new_error <= keep_error){
          keep <- new_node
          keep_error <- new_error
        }

        #Remove leaves
        if(nchar(node[[i]]) > 1){
            new_node <- node[unlist(lapply(node,function(x) substr(x,1,nchar(node[[i]])-1) != substr(node[[i]],1,nchar(node[[i]])-1)))]
          new_node <- append(new_node,list(substr(node[[i]],1,nchar(node[[i]])-1)))
          new_error <- get_error_node(tab_train,tab_val,new_node,nrow(val))
          if(new_error < keep_error){
            keep <- new_node
            keep_error <- new_error
          }
        }
      }

      if(is.null(keep))
        break
      else if(error == keep_error){
        node <- keep
        node <- node[sample(1:length(node),length(node),F)]
      }
      else{
        node <- keep
        error <- keep_error
        node <- node[sample(1:length(node),length(node),F)]
      }
    }
    node <- clean_node(node,tab_train)

    if(is.null(current_node)){
      current_node <- node
      current_error <- error
    }
    else if(error < current_error){
      current_node <- node
      current_error <- error
    }
  }

  node <- current_node
  error <- get_error_node(tab_train,tab_val,node,nrow(val))
  hhat <- get_hypothesis(tab_train,node)

  #Predict
  predict <- function(string){
    y <- hhat$hx[unlist(lapply(as.list(hhat$x),function(x) string_start(string = string,start = x)))]
    return(ifelse(y == "Tie",0,y))
  }

  end_time <- Sys.time()

  #Return
  return(list("error" = error,"hhat" = hhat,"predict" = predict,"time" = end_time - start_time))
}
