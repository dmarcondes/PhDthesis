#' @import plyr
#' @import rlist
#' @import fst
#' @import parallel
#' @import numbers
#' @import stringi
#' @import ggplot2
#' @import ggthemes
#' @import utils
#' @import stats
#' @export
#' @title U-curve algorith on the Partition Lattice Learning Space
#'
#' @description TBD
#'
#' @details TBD
#'
#' @param xtrain Vector with the training sample of x values.
#' @param ytrain Vector with the training sample of y values.
#' @param xval Vector with the validation sample of x values.
#' @param yval Vector with the validation sample of y values.
#' @param X Domain of variable x.
#' @param Y Domain of variable y.
#' @param optimal Logical indicating if the algorithm should return an optimal solution.
#' @param increasing Logical indicating if target function is increasing.
#' @param earlystop Logical indicating if it should be possible to early stop the algorithm.
#' @param exhaust Number of points to exhaust before stopping the algorithm.
#' @param sampleNeigh Either false to consider all neighboors, or the maximum number of neighboors to sample at each exhaustion. If a number,then optimal should be false.
#' @param verbose Logical to print a trace of the algorithm.
#' @param stop Number of nodes yet to evaluate to trigger exhaustive search.
#' @param path Path to preprocessed partition files.
#' @param Lh A data frame with the partition lattice.
#' @param cores Number of cores for parallel computing.
#' @param save Logical. Whether to save the nodes exhausted.
#' @return \item{hypotheses}{The estimated hypothesis of the global minimums with least VC dimension.}
#' @return \item{partitions}{Partitions of the global minimums with least VC dimension.}
#' @return \item{error}{Validation error of the global minimums.}
#' @return \item{exhausted}{Number of nodes exhausted during algorithm.}
#' @return \item{remain}{Number of nodes remaining after algorithm stopped.}
#' @return \item{finished}{If the algorithm was finished or ended after not finding any Strong Local Minimum.}
#' @return \item{SLMvis}{Number of nodes exhasuted until the last Strong Local Minimum was found.}
#' @return \item{remain_after_prune}{Number of nodes remaining after finding each Strong Local Minimum.}
#' @return \item{exhausted_until_prune}{Number of nodes exhausted until finding each Strong Local Minimum.}
#' @return \item{optimal}{Wheter an optimal solution was returned.}
#' @return \item{plot}{Plot with the error of each Strong Local Minimum.}
#' @examples
#' set.seed(1)
#' x <- sample(x = c("01","02","03","04","05","06","07","08","09","10"),size = 50,replace = TRUE)
#' y <- as.factor(ifelse(as.numeric(x)-5+rnorm(50,0,5/3) > 0,1,0))
#' x <- factor(x)
#' train <- sample(1:50,35,FALSE)
#' xtrain <- x[train]
#' ytrain <- y[train]
#' xval <- x[!(c(1:50) %in% train)]
#' yval <- y[!(c(1:50) %in% train)]
#' u <- ucurve(xtrain,ytrain,xval,yval,optimal = FALSE,sampleNeigh = 5000,exhaust = 10,cores = 1)

ucurve <- function(xtrain,ytrain,xval,yval,X,Y,optimal = T,exhaust = 1000,sampleNeigh = F,increasing = F,earlystop = F,
                   verbose = T,stop = 0,path = "~/GDrive/Doutorado/Códigos/Particoes/",Lh = NULL,cores = 4,save = F){

  #Test if parameters correct
  if(!!sampleNeigh){
    if(optimal)
      stop("Can only sample neighboors for suboptimal algorithm.")
  }

  if(increasing){
    if(optimal)
      stop("When target increasing can only return suboptimal solution.")
    if(!!sampleNeigh)
      stop("When target increasing all neighboors are considered.")
  }

  #Timing
  start_time <- Sys.time()

  #Get sample info
  #X <- unique(c(as.character(xtrain),as.character(xval))) #Domain
  #X <- X[order(X)] #Order domain
  n <- length(X) #Length of domain
  if(optimal) #Maximum domain size for optimal solution
    if(length(X) > 12)
      stop("Sorry, but this algorithm is not scalable! It works for at most 12 points in X domain.")
  #Y <- unique(c(as.character(ytrain),as.character(yval))) #Image
  #Y <- Y[order(Y)] #Order Image
  if(length(Y) != 2) #Only binary classification problems
    stop("The algorithm only work for binary classification problems.")

  #Turn sample into factor
  xtrain <- factor(xtrain,X)
  xval <- factor(xval,X)
  ytrain <- factor(ytrain,Y)
  yval <- factor(yval,Y)
  x <- NULL
  y <- NULL

  #Information about the algorithm
  if(verbose){
    cat("------------------------------------------------------------------------------\n")
    cat(paste(length(unique(X)),"points in the domain\n"))
    cat(paste(ifelse(n <= 218,format(bell(n),big.mark = ",",scientific = F),"???"),
                            "nodes on the Partition Lattice Learning Space\n"))
    cat(paste("Training sample size:",length(xtrain),"\n"))
    cat(paste("Validation sample size:",length(xval),"\n"))
    cat(paste("Returning",ifelse(increasing,"increasing suboptimal",ifelse(optimal,"optimal","suboptimal")),
              "hypotheses\n"))
    cat("------------------------------------------------------------------------------\n")
    cat("\n")
  }

  #Delete file with visited nodes
  suppressWarnings(system("rm -f visited.dat"))

  #Get parameters
  search <- 1 #Start search
  jtrain <- jointDistribution(xtrain,ytrain) #Train joint distribution
  jval <- jointDistribution(xval,yval) #Validation joint distribution
  part <- makePartition(list(1:n)) #First part to start algorithm
  strongMinimums <- vector() #Declare vector of strong minimums
  vis <- 0 #Start number of exhasted nodes
  SLMvis <- 0 #Exhausted nodes until last Strong Local Minimum found
  exhausted_until_prune <- vector() #Nodes exhausted until each Strong Local Minimum was found
  min_error <- NULL
  Plot <- NULL
  if(optimal)
    remain_after_prune <- vector() #Nodes remaining after prunning when found each Strong Local Minimum
  else
    remain_after_prune <- NA

  #If optimal, get partitions
  if(optimal){
    if(is.null(Lh)){
      cat("\n")
      cat("Reading partitions...")
      cat("\n")
      Lh <- read.fst(paste(path,"part_",n,".fst",sep = ""))
      Lstore <- Lh
    }
    else
      Lstore <- Lh
  }

  if(verbose){
    cat("------------------------------------------------------------------------------\n")
    cat("Starting U-curve algorithm...\n")
    cat("------------------------------------------------------------------------------\n")
  }

  if(verbose)
    cat("*")

  #Start algorithm
  while(search){
    if(optimal){
      #If there is no more nodes to restart from end algorithm
      if(length(colnames(Lh)) == 0)
        break
    }

    #Correct possible errors
    part <- makePartition(part)

    #Calculate error
    errorPart <- getError(part,jtrain,jval,increasing,n)

    #Add to nodes exhausted
    vis <- vis + 1

    #Save as visited node
    if(save)
      addNode(part,errorPart)

    #Get neighbours
    if(increasing)
      N <- findNeighborsIncreasing(part,n,cores)
    else{
      if(!sampleNeigh)
        N <- findNeighbors(part,n,cores)
      else
        N <- as.vector(sampleNeighbors(sample_size = sampleNeigh,part = part,n = n,cores = cores))
    }

    #Error of neighboors
    err <- unlist(mclapply(data.frame(rbind(N)),function(x) getError(getPartition(x),jtrain,jval,increasing,n),mc.cores = cores))
    if(optimal)
      Lh <- Lh[,!(colnames(Lh) %in% N[err > errorPart])]

    #If the node is a Strong Local Minimum
    if(errorPart <= min(err)){
      #Plot
      if(is.null(min_error))
        min_error <- errorPart
      else
        min_error <- c(min_error,errorPart)
      p <- data.frame("x" = 1:(length(strongMinimums) + 1),y = min_error)
      if(length(strongMinimums) > 0 & verbose){
        Plot <- ggplot(p,aes(x = x,y = y,group = 1)) + theme_linedraw() + geom_point(color = "white") +
          geom_line(color = "white") +
          xlab("Strong Local Minimum") + ylab("Validation Error") + ylim(c(0.95*min(min_error),1.05*max(min_error))) +
          scale_x_continuous(breaks = 1:1e4) + geom_hline(yintercept = min(min_error),linetype = "dashed",color = "white") +
          theme_solarized(light = FALSE) + theme(strip.background = element_blank(),
                                                 strip.text = element_text(size = 20,face = "bold",color = "white")) +
          theme(legend.title = element_text(face = "bold"),legend.position = "none") +
          theme(plot.title = element_text(face = "bold",size = 25,color = "white",hjust = 0.5),
                axis.text.x = element_text(size = 15,face = "bold",color = "white"),
                axis.text.y = element_text(size = 15,face = "bold",color = "white"),
                legend.box.margin = unit(x=c(20,0,0,0),units="mm"),
                legend.key.width=unit(3.5,"cm"),panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.title = element_text(color = "white",size = 20),
                plot.caption = element_text(face = "bold",color = "white",hjust = 0,size = 15)) +
          theme(plot.margin = unit(c(1,1,1,1), "lines")) +
          ggtitle(paste(ifelse(optimal,"Optimal","Subotimal"),"U-curve Algorithm \nPartition Lattice Learning Space"))
        print(Plot)
      }

      #Nodes exhausted until found the Strong Local Minimum
      SLMvis <- vis

      #Store node as Strong Local Minimum
      strongMinimums <- c(strongMinimums,namePartition(part))

      if(optimal){
        if(is.data.frame(Lh)){
          #Erase nodes associated to it
          num_part <- Lstore[,colnames(Lstore) == namePartition(part)]
          Lh <- Lh[,!unlist(mclapply(Lh,function(x) is.related(x,num_part),mc.cores = cores))]
        }
      }

      if(verbose){
        if(optimal){
          cat("\n")
          cat("Found a SLM. There are ",ncol(Lh)," (",ifelse(n <= 218,round(100*ncol(Lh)/bell(n),5),0),"%) nodes remaining after ",
              vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%) exhausted",sep = "")
          cat("\n")
        }
        else{
          cat("\n")
          cat("Found a SLM after exhausting ",vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = "")
          cat("\n")
        }
      }

      if(optimal){
        #Nodes remaining after pruning
        remain_after_prune <- c(remain_after_prune,ncol(Lh))
      }

      #Nodes exhausted until found Strong Local Minimum
      exhausted_until_prune <- c(exhausted_until_prune,vis)

      if(optimal){
        if(length(colnames(Lh)) <= stop){
          #If less nodes remaining than stop, stop algorithm
          search <- 0
        }
        else{
          #Choose a new node to restart
          part <- sample(colnames(Lh),1) #Sample new part
          part <- getPartition(part) #Get partition
        }
      }
      else{
        #Sample new partition to restart
        part <- samplePartition(n,increasing)
      }
      #If exhausted more points than desired
      if(!optimal){
        if(vis >= exhaust)
          break
      }

      if(earlystop){
        cont <- waitForKey(message = "Stop?",time = 3)
        if(cont == "key")
          break
      }
      if(verbose)
        cat("*") #SLM symbol
    }
    else{
      N1 <-  N[err < errorPart] #Get neighboors with less loss

      if(optimal){
        #Erase last node
        Lh <- Lh[,colnames(Lh) != namePartition(part)]

        #Neighboors with lesser loss still on space
        N1 <- N1[N1 %in% colnames(Lh)] #Get neighboors still on the space
      }

      #If some neighboor is available to continue, sample next node from it
      if(length(N1) > 0){
          part <- getPartition(sample(N1,1))
          if(verbose)
            cat("-") #Next node symbol
      }
      else{
        if(optimal){
          if(length(colnames(Lh)) > 0){
            #Choose a new node to restart
            part <- sample(colnames(Lh),1) #Sample new part
            part <- getPartition(part) #Get partition
          }
          else
            break #End algorithm
        }
        else{
          #Sample new partition to restart
          part <- samplePartition(n,increasing)
        }

        #If exhausted more points than desired
        if(!optimal){
          if(vis >= exhaust)
            break
        }
        if(verbose)
          cat("*") #SLM symbol
      }
    }
  }

  if(optimal){
    remain <- length(colnames(Lh)) #Number of nodes remaining

    #If nodes remaining start exhaustion
    if(remain > 0 & remain <= stop){
      finished <- 1
      if(verbose){
        cat("\n")
        cat(paste("Ending U-curve algorithm and starting exhaustion of ",ncol(Lh)," (",
                  ifelse(n <= 218,round(100*ncol(Lh)/bell(n),5),0),
                  "%) nodes remaining...",sep = ""))
        cat("\n")
      }
      #Get error of all remaining
      err <- unlist(mclapply(data.frame(rbind(colnames(Lh))),function(x) getError(getPartition(x),jtrain,jval,increasing,n),mc.cores = cores))

      #Get minimum error of strong local minimums
      m <- min(unlist(mclapply(data.frame(rbind(strongMinimums)),function(x) getError(getPartition(x),jtrain,jval,increasing,n),mc.cores = cores)))

      #Add nodes with less error to strong local minimums
      strongMinimums <- c(strongMinimums,colnames(Lh)[err <= m])

      if(verbose){
        cat("\n")
        cat("---------------------------------------------------------------------------------")
        cat("\n")
        cat(paste("There was ",remain," (",ifelse(n <= 218,round(100*remain/bell(n),5),0),"%)"," nodes to visit after exhausting ",
                  vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = ""))
        cat("\n")
        cat("---------------------------------------------------------------------------------")
      }
    }
    else if(remain == 0){
      finished <- 1
      if(verbose){
        cat("\n")
        cat("---------------------------------------------------------------------------------")
        cat("\n")
        cat(paste("There was ",remain," (",ifelse(n <= 218,round(100*remain/bell(n),5),0),"%)"," nodes to visit after exhausting ",
                  vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = ""))
        cat("\n")
        cat("---------------------------------------------------------------------------------")
      }
    }
    else{
      finished <- 0
      if(verbose){
        cat("\n")
        cat("---------------------------------------------------------------------------------")
        cat("\n")
        cat(paste("There was ",remain," (",ifelse(n <= 218,round(100*remain/bell(n),5),0),"%)"," nodes not visited early stopping afeter exhausting ",
                  vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = ""))
        cat("\n")
        cat("---------------------------------------------------------------------------------")
      }
    }
  }
  else{
    finished <- 1
    if(verbose){
      cat("\n")
      cat(paste("U-curve algorithm returned a suboptimal solution after the exhaustion of ",vis," (",
                ifelse(n <= 218,round(100*vis/bell(n),5),0),
                "%) nodes...",sep = ""))
      cat("\n")
    }
  }

  #Solution
  e <- apply(rbind(strongMinimums),2,function(x) getError(getPartition(x),jtrain,jval,increasing,n)) #Error of SLM
  global <- strongMinimums[e == min(e)] #Minimums of SLMs
  global <- apply(rbind(global),2,function(x) getPartition(x)) #Get partition of minimums
  solution <- unique(unlist(lapply(global,function(x) namePartition(makePartition(x))))) #Get name of minimums partition

  #Solution
  oh <- tapply(cbind(solution),solution,function(x) optimalHyp(part = x,jtrain = jtrain,increasing = increasing,n)) #Get hypotheses estimated of solutions
  names(oh) <- NULL
  if(!increasing)
    oh <- lapply(oh,function(x) data.frame("X" = X,"fx" = plyr::mapvalues(factor(x[,2]),c("0","1"),Y)))
  else
    oh <- unlist(oh,recursive = F)
  error <- unique(tapply(cbind(solution),solution,function(x) getError(part = getPartition(x),jtrain = jtrain,jval = jval,increasing,n))) #Error of solution
  names(error) <- NULL

  #Timing
  end_time <- Sys.time()

  #Classifier
  classifier <- function(x){
    y <- unlist(lapply(unique(oh),function(y) y[y[,1] == x,2]))
    y <- table(y)
    y <- names(y)[y == max(y)]
    if(length(y) > 1)
      y <- sample(y,1)
    return(y)
  }

  return(list("classifier" = classifier,"hypotheses" = unique(oh),"partitions" = solution,"error" = error,
              "exhausted" = vis,"remain" = ifelse(optimal,remain,NA),"finished" = finished,"SLMvis" = ifelse(!finished,SLMvis,NA),
              "remain_after_prune" = remain_after_prune,"exhausted_until_prune" = exhausted_until_prune,
              "time" = difftime(end_time,start_time,units = "mins"),"optimal" = optimal,"plot" = Plot))
}
