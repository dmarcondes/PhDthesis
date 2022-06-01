####Utility functions for partitionUcurve package#####

#Test if partion and order parts by minimum
makePartition <- function(part){
  #Is partition
  v <- unlist(part)
  if(!(min(v) == 1 & max(v) == length(v) & length(unique(v)) == length(v)))
    stop("It is not a partition. Please review the code.")

  #Order partition
  part <- part[order(unlist(lapply(part,min)))]
  part <- lapply(part,function(x) x[order(x,decreasing = F)])

  #Name partition
  names(part) <- unlist(lapply(part,function(x) paste(x,collapse = ",")))

  return(part)
}

#Get name of partition
namePartition <- function(part){
  return(paste(names(part),collapse = "-"))
}

#Get name from index
nameFromIndex <- function(index){
  n <- length(index)
  name <- namePartition(makePartition(tapply(c(1:n),index,function(x) x)))

  return(name)
}

#Get partition from name
getPartition <- function(name){
  part <- sapply(unlist(strsplit(x = name,split = "-")),function(x) lapply(strsplit(x,","),as.numeric))
  return(part)
}

#Test if two partitions are related
is.related <- function(part1,part2){
  if(all.equal(part1,part2) == T)
    return(T)
  m1 <- max(part1)
  m2 <- max(part2)
  if(m1 == m2)
    return(F)
  if(m1 > m2)
    return(sum(!unlist(lapply(tapply(part2,part1,function(x) x),function(x) length(unique(x)) == 1))) == 0)
  if(m1 < m2)
    return(sum(!unlist(lapply(tapply(part1,part2,function(x) x),function(x) length(unique(x)) == 1))) == 0)
}

numberNeighboors <- function(part){
  #Number of blocks unions
  #numb_union <- (length(part)*(length(part) - 1))/2

  #Number of blocks breaks
  part <- lapply(part,length)
  part <- part[part > 1]
  if(length(part) > 0)
    numb_break <- sum(unlist(lapply(part,function(x) ifelse(x > 2,ifelse(x %% 2 == 0,sum(choose(x,1:(x/2 - 1))) + choose(x,x/2)/2,
                                  sum(choose(x,1:floor(x/2)))),1))))
  else
    numb_break <- 0
  return(list("total" = numb_break))
}

numberNeighboorsBreakBlock <- function(size){
  #Number of block's breaks
  if(size > 1)
    numb_break <- ifelse(size > 2,ifelse(size %% 2 == 0,sum(choose(size,1:(size/2 - 1))) + choose(size,size/2)/2,
                                                                         sum(choose(size,1:floor(size/2)))),1)
  else
    numb_break <- 0
  return(numb_break)
}

#Find neighbors of a partition
findNeighbors <- function(part,n,cores){

  #All breaks of partition
  m <- matrix(ncol = n)
  for(i in 1:length(part))
    m[1,part[[i]]] <- i
  size <- max(m)
  nsize <- size + 1

  if(max(m) < n){
    for(i in 1:max(m)){
      if(length(m[1,m[1,] == i]) > 1){
        mtemp <- expand.grid(rlist::list.append(c(1),replicate(length(m[1,m[1,] == i])-1, c(0,1), simplify=FALSE)))
        mtemp <- mtemp[-nrow(mtemp),]
        mtemp <- as.matrix(nsize*mtemp)
        mtemp[mtemp == 0] <- i
        m2 <- matrix(0,nrow = nrow(mtemp),ncol = n)
        m2[,m[1,] == i] <- as.matrix(mtemp)
        m3 <- t(replicate(nrow(m2),m[1,]))
        m3[,m[1,] == i] <- 0
        m2 <- m2 + m3
        m <- rbind(m,m2)
      }
    }
    breaks <- unlist(mclapply(data.frame(t(m)),nameFromIndex,mc.cores = cores))[-1]
  }
  else
    breaks <- NULL

  #All unions of partition
  # if(length(part) > 1){
  #   nameUnion <- combn(1:length(part),2)
  #   if(ncol(nameUnion) > 1){
  #     union <- plyr::alply(nameUnion,2,function(x) c(x[[1]],x[[2]],part[[x[1]]],part[[x[2]]]))
  #     union <- mclapply(union,function(x) rlist::list.append(list(x[-c(1,2)]),part[-c(x[1:2])]),mc.cores = cores)
  #     union <- mclapply(union,function(y) lapply(rapply(y, enquote, how="unlist"), eval),mc.cores = cores)
  #     union <- mclapply(union,makePartition,mc.cores = cores)
  #     union <- unlist(lapply(union,namePartition))
  #     names(union) <- NULL
  #   }
  #   else
  #     union <- namePartition(makePartition(list(c(part[[1]],part[[2]]))))
  #   if(is.null(breaks))
  #       return(union)
  #   else
  #     return(c(breaks,union))
  # }
  # else
  #   return(breaks)
  return(breaks)
}

#Sample neighbors of a partition
sampleNeighbors <- function(sample_size,part,n,cores){
  numb <- numberNeighboors(part)
  sbreak <- sample_size
  if(numb$total < sample_size)
    return(findNeighbors(part,n,cores))

  #Sample sizes
  # if(numb$breaks < numb$union){
  #   sbreak <- floor(sample_size * numb$breaks/numb$total)
  #   sbreak <- ifelse(sbreak <= numb$breaks,sbreak,numb$breaks)
  #   sunion <- sample_size - sbreak
  # }
  # else{
  #   sunion <- floor(sample_size * numb$union/numb$total)
  #   sunion <- ifelse(sunion <= numb$union,sunion,numb$union)
  #   sbreak <- sample_size - sunion
  # }

  #All breaks of partition
  if(sbreak == numb$total){
    m <- matrix(ncol = n)
    for(i in 1:length(part))
      m[1,part[[i]]] <- i
    size <- max(m)
    nsize <- size + 1

    if(max(m) < n){
      for(i in 1:max(m)){
        if(length(m[1,m[1,] == i]) > 1){
          mtemp <- expand.grid(rlist::list.append(c(1),replicate(length(m[1,m[1,] == i])-1, c(0,1), simplify=FALSE)))
          mtemp <- mtemp[-nrow(mtemp),]
          mtemp <- as.matrix(nsize*mtemp)
          mtemp[mtemp == 0] <- i
          m2 <- matrix(0,nrow = nrow(mtemp),ncol = n)
          m2[,m[1,] == i] <- as.matrix(mtemp)
          m3 <- t(replicate(nrow(m2),m[1,]))
          m3[,m[1,] == i] <- 0
          m2 <- m2 + m3
          m <- rbind(m,m2)
        }
      }
      breaks <- unlist(mclapply(data.frame(t(m)),nameFromIndex,mc.cores = cores))[-1]
    }
    else
      breaks <- NULL
  }
  else{
    m <- matrix(ncol = n)
    for(i in 1:length(part))
      m[1,part[[i]]] <- i
    size <- max(m)
    nsize <- size + 1
    avai_blocks <- names(table(m[1,]))[table(m[1,]) > 1]
    wheight <- table(m[1,])[table(m[1,]) > 1]
    wheight <- apply(cbind(wheight),1,numberNeighboorsBreakBlock)

    remains <- sbreak
    sampled <- 0
    mtemp2 <- NULL
    while(sampled < sbreak){
      mtemp <- mclapply(data.frame(rbind(c(1:remains))),function(x){
        #Sample block
        block <- sample(x = avai_blocks,size = 1,prob = wheight)
        bsize <- sum(m[1,] == block)

        #Size new
        if(bsize %% 2 == 0 & bsize > 2)
          p <- c(choose(bsize,1:(bsize/2 - 1)),choose(bsize,bsize/2)/2)
        else if(bsize %% 2 != 0 & bsize > 2)
          p <- choose(bsize,1:floor(bsize/2))
        else if(bsize == 2)
          p <- 1
        size_new <- sample(x = 1:floor(bsize/2),size = 1,prob = p)

        #Which to turn
        turn <- sample(x = c(1:ncol(m))[m[1,] == block],size = size_new)
        if(min(turn) == min(c(1:ncol(m))[m[1,] == block]))
          turn <- c(1:ncol(m))[m[1,] == block][!(c(1:ncol(m))[m[1,] == block] %in% turn)]

        #New block
        mtemp <- m[1,]
        mtemp[turn] <- nsize
        return(as.list(mtemp))},mc.cores = cores)
      mtemp2 <- unique(rbind(mtemp2,as.matrix(data.table::rbindlist(mtemp))))
      sampled <- nrow(mtemp2)
      remains <- sbreak - sampled
    }
    m <- mtemp2

    if(nrow(m) > 0)
      breaks <- unlist(mclapply(data.frame(t(m)),nameFromIndex,mc.cores = cores))
    else
      breaks <- NULL
  }

  #All unions of partition
  # if(length(part) > 1 & sunion > 1){
  #   nameUnion <- combn(1:length(part),2)
  #   nameUnion <- cbind(nameUnion[,sample(x = 1:ncol(nameUnion),size = sunion,replace = F)])
  #   if(ncol(nameUnion) > 1){
  #     union <- plyr::alply(nameUnion,2,function(x) c(x[[1]],x[[2]],part[[x[1]]],part[[x[2]]]))
  #     union <- mclapply(union,function(x) rlist::list.append(list(x[-c(1,2)]),part[-c(x[1:2])]),mc.cores = cores)
  #     union <- mclapply(union,function(y) lapply(rapply(y, enquote, how="unlist"), eval),mc.cores = cores)
  #     union <- mclapply(union,makePartition,mc.cores = cores)
  #     union <- unlist(lapply(union,namePartition))
  #     names(union) <- NULL
  #   }
  #   else{
  #     union <- namePartition(makePartition(rlist::list.append(part[-nameUnion[,1]],unlist(part[nameUnion[,1]],use.names = F))))
  #   }
  #
  #   if(is.null(breaks))
  #     return(union)
  #   else
  #     return(c(breaks,union))
  # }
  # else
  return(breaks)
}

#Joint distribution
jointDistribution <- function(x,y){
  return(prop.table(table(x,y)))
}

#Calculate error of a partition
getError <- function(part,jtrain,jval,increasing,n){

  if(increasing){
    if(!is.increasing(part))
      stop("Partition is not increasing.")
    M <- unlist(lapply(part,max))
    M <- c(0,M)
    image <- lapply(data.frame(rbind(M)),function(x) c(rep(2,x),rep(1,n-x)))
    error <- unlist(lapply(image,function(x) sum(unlist(lapply(data.frame(rbind(1:nrow(jtrain))),function(y) jtrain[y,x[y]])))))
    image <- image[error == min(error)]
    error <- unlist(lapply(image,function(x) sum(unlist(lapply(data.frame(rbind(1:nrow(jval))),function(y) jval[y,x[y]])))))
    error <- min(error)

    return(error)
  }

  else{
    #Estimated hypothesis
    parJoint <- lapply(part,function(x) colSums(rbind(jtrain[x,],c(0,0))))
    optimalComplement <- unlist(lapply(parJoint,function(x) c(1:2)[x == min(x)][1]))

    #Error
    parVal <- unlist(lapply(part,function(x) colSums(rbind(jval[x,],c(0,0)))))[optimalComplement + seq(0,2*(length(part)-1),2)]
    error <- sum(parVal)

    return(error)
  }
}

#Add a visited node to file
addNode <- function(part,error){
  name <- paste(namePartition(part),"e",error,sep = "")
  write(name,file = "visited.dat",append = T)
}

#Test if node has been visited and return its error
visitedNode <- function(part){
  name <- namePartition(part)
  vis <- suppressWarnings(system(paste("grep ",name," visited.dat"),intern = T))
  if(length(vis) > 0)
    vis <- as.numeric(unlist(strsplit(vis,"e"))[2])
  else
    vis <- NULL

  return(vis)
}

#Get optimal hypothesis of partition
optimalHyp <- function(part,jtrain,increasing,n){
  part <- getPartition(part)

  if(increasing){
    if(!is.increasing(part))
      stop("Partition is not increasing.")
    M <- unlist(lapply(part,max))
    M <- c(0,M)
    image <- lapply(data.frame(rbind(M)),function(x) c(rep(2,x),rep(1,n-x)))
    error <- unlist(lapply(image,function(x) sum(unlist(lapply(data.frame(rbind(1:nrow(jtrain))),function(y) jtrain[y,x[y]])))))
    image <- lapply(image,function(x) 2-x)
    image <- image[error == min(error)]
    f <- lapply(image,function(x) data.frame("X" = rownames(jtrain),"fx" = x))
    names(f) <- NULL
  }
  else{
    #Estimated hypothesis
    parJoint <- lapply(part,function(x) colSums(rbind(jtrain[x,],c(0,0))))
    optimal <- unlist(lapply(parJoint,function(x) c(0:1)[x == max(x)][1]))

    #Hypothesis
    domain <- unlist(strsplit(names(optimal),","))
    domain <- domain[order(as.numeric(domain))]
    f <- vector()
    for(i in domain)
      f[domain == i] <- optimal[unlist(lapply(strsplit(names(optimal),","),function(x) i %in% x))]
    f <- data.frame("X" = domain,"fx" = f)
  }

  return(f)
}

#Sample a partition
samplePartition <- function(n,increasing){
  #Size of new partion
  size <- 1 + rbinom(1,n-1,1/2)

  if(size > 1){
    #Sample an ordenation of the domain
    if(increasing)
      o <- 1:n
    else
      o <- sample(1:n,n,F)

    #Sample last position of each partition
    pos <- sample(1:(n-1),size-1,F)
    pos <- pos[order(pos,decreasing = F)]

    #Creating part
    part <- list()
    pos[size] <- n
    for(i in 1:size)
      part[[i]] <- o[ifelse(i == 1,1,pos[i-1]+1):pos[i]]
    part <- makePartition(part)
  }
  else
    part <- makePartition(list(1:n))
  return(part)
}

#Test if partition is increasing
is.increasing <- function(x){
  if(length(x) > 1){
    m <- unlist(lapply(x,min))
    m <- m[2:length(m)]
    M <- unlist(lapply(x,max))
    M <- M[1:(length(M)-1)]
    return(sum(M < m) == length(M))
  }
  return(T)
}

#Break partition after a value
breakAfter <- function(part,x,cores){
  p <- unlist(mclapply(part,function(y) x %in% y,mc.cores = cores))
  pos <- c(1:length(part[p][[1]]))[part[p][[1]] == x]
  p1 <- part[p][[1]][1:pos]
  p2 <- part[p][[1]][(pos+1):length(part[p][[1]])]
  npart <- makePartition(list.append(part[!p],p1,p2))

  return(namePartition(npart))
}

#Union part with next one
unionNext <- function(part,x,cores){
  p <- part[-c(x,x+1)]
  p0 <- unlist(part[c(x,x+1)])
  p <- makePartition(list.append(p,p0))

  return(namePartition(p))
}

#Find increasing neighbors of a partition
findNeighborsIncreasing <- function(part,n,cores){

  #All increasing breaks of partition
  if(length(part) < n){
    b <- as.list(unlist(mclapply(part,function(x) x[ifelse(length(x) == 1,0,1):(length(x) - 1)],mc.cores = cores)))
    b <- unlist(mclapply(b,function(x) breakAfter(part,x,cores)))
  }
  else
    b <- NULL

  #All increasing unions of partition
  # if(length(part) > 1){
  #   u <- unlist(mclapply(data.frame(rbind(1:(length(part)-1))),function(x) unionNext(part,x,cores)))
  # }
  # else
  #   u <- NULL
  #
  # if(is.null(b))
  #   return(u)
  # else if(is.null(u))
  return(b)
  # else
  #   return(c(b,u))
}

#Code from https://github.com/Dans-labs/R-package_EmilMisc
waitForKey <- function(message='Continuing in {n} seconds, or press any key.', time=10, counter=.5, precision=.01) {
  stopifnot(is.character(message), length(message)==1, !is.na(message),
            is.numeric(time), length(time)==1, !is.na(time), time>0,
            is.numeric(counter), length(counter)==1, !is.na(counter), counter>=0,
            is.numeric(precision), length(precision)==1, !is.na(precision), precision>0)
  my_in <- file('stdin')
  open(my_in, blocking=FALSE)
  ans <- readLines(my_in, n=1)
  if(counter==0) cat(message)
  while(time>0 && !length(ans)) {
    if(counter>0 && !is.null(message) && round(time/counter, digits = 3) %% 1L==0L) {
      cat(gsub('\\{n\\}', format(round(time, digits=4), width=5, scientific = F), message), '\r', sep = '')
    }
    Sys.sleep(precision)
    time <- time-precision
    ans <- readLines(my_in, n=1)
  }
  close(my_in)
  if(length(ans)) {
    return(invisible('key'))
  } else {
    return(invisible('timer'))
  }
}

#####W Operator####

#Apply filter to training data
# applyFilter <- function(w,pt,xtrain,d,wsize){
#   f <- function(img){
#     g <- Vectorize(function(i,j) pt[[paste(as.vector(apply(w,1,function(x) ifelse(i + x[1] <= d & j + x[2] <= d,img[i + x[1],j + x[2]],0))),collapse = "+")]])
#     return(outer(X = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),Y = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),
#           FUN = g))
#   }
#   f2 <- function(img){
#     g <- Vectorize(function(i,j) pt[[paste(as.vector(apply(w,1,function(x) ifelse(i + x[1] <= d & j + x[2] <= d,xtrain[img,i + x[1],j + x[2]],0))),collapse = "+")]])
#     return(outer(X = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),Y = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),
#                  FUN = g))
#   }
#   new_xtrain <- apply(xtrain,1,f)
#   tmp <- aperm(xtrain,c(2,3,1))
#   a <- mclapply(1:nrow(xtrain),f2,,mc.cores = 4)
#
#   return(new_xtrain)
# }
#
# #Get error of multi layer W operator
# getWError <- function(W,xtrain,xval,Cytrain,Cyval,nlayer,d,cor_dim,wsize){
#   if(length(W) > 1){
#     domain <- apply(expand.grid(replicate(nrow(W[[1]]), c(0,1), simplify=FALSE)),1,function(x) paste(x,collapse = "+"))
#     tmp <- expand.grid(replicate(length(domain), c(0,1), simplify=FALSE))
#     colnames(tmp) <- domain
#     domain <- tmp
#     for(i in 1:nrow(domain)){
#       pt <- domain[i,]
#       new_xtrain <- applyFilter(w = W[[1]],pt = pt,xtrain = xtrain,d = d,cor_dim = cor_dim,wsize = wsize)
#     }
#   }
#   else{
#
#   }
# }
#
# #Apply a filter until layer l
# applyFilter <- function(x,W,layer,cor_dim){
#   for(l in 1:layer){
#     dic <- lapply(strsplit(x = names(W[[l]])[-length(W[[l]])],split = "+"),function(x) ifelse(as.numeric(c(x[1],x[3])) == 2,-1,as.numeric(c(x[1],x[3]))))
#     W[[l]] <- cbind(W[[l]],apply(W[[l]][1:(ncol(W[[l]])-1)],1,function(x) paste(x,collapse = "")))
#     f <- function(img){
#       g <- Vectorize(function(i,j) W[[l]][W[[l]][ncol(W[[l]])] == paste(unlist(lapply(dic,function(pos) ifelse(i + pos[1] > 0 & i + pos[1] <= d & j + pos[2] <= d & j + pos[2] > 0,
#                                                                           x[img,i + pos[1],j + pos[2]],0))),collapse = ""),ncol(W[[l]])-1])
#       return(outer(X = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),Y = (1+floor(wsize/2)):(cor_dim - floor(wsize/2)),
#                FUN = g))
#       }
#     x <- lapply(1:nrow(x),f)
#     if(nrow(x[[1]]) > 1){
#       x <- aperm(simplify2array(x),c(3,1,2))
#       cor_dim <- nrow(x[1,,])
#     }
#     else
#       x <- unlist(x)
#   }
#   return(x)
# }
#
# #Get error of given window
# getWError <- function(x,y,w,cd){
#   ypred <- applyFilter(x = x,W = w,layer = length(w),cor_dim = cd)
#   e <- 1-sum(diag(prop.table(table(ypred,y))))
#   return(e)
# }
#
# #Update W
# updateW <- function(W,xtrain_batch,ytrain_batch,cor_dim){
#   err <- data.frame("layer" = NA,"pos" = NA,"erro" = NA)
#   for(l in length(W):1){
#     if(l > 1){
#       x <- applyFilter(xtrain_batch,W,l-1,cor_dim)
#       tmp_err <- vector()
#       for(t in 1:nrow(W[[l]])){
#         w <- W[l:(length(W))]
#         w[[1]][t,ncol(w[[1]])] <- ifelse(w[[1]][t,ncol(w[[1]])] == 0,1,0)
#         tmp_err[[t]] <- getWError(x,ytrain_batch,w,ncol(x[1,,]))
#       }
#       pos <-sample(x = c(1:length(tmp_err))[tmp_err == min(tmp_err)],1)
#       err <- rbind.data.frame(err,data.frame("layer" = l,"pos" = pos,"erro" = min(tmp_err)))
#     }
#   }
#   err <- na.omit(err)
#   pos <- sample(x = c(1:nrow(err))[err$erro == min(err$erro)],1)
#   err <- err[pos,]
#   W[[err$layer]]$h[err$pos] <- ifelse(W[[err$layer]]$h[err$pos]  == 0,1,0)
#
#   return(list("W" = W,"err" = err$erro))
# }
#
