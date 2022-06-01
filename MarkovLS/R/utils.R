#Get error
get_error_node <- function(tab_train,tab_val,node,total){
  part_train <- vector(length = nrow(tab_train))
  part_val <- vector(length = nrow(tab_val))
  for(i in 1:length(node)){
    part_train[substr(rownames(tab_train),1,nchar(node[[i]])) == node[[i]]] <- i
    part_val[substr(rownames(tab_val),1,nchar(node[[i]])) == node[[i]]] <- i
  }

  error <- 0
  for(i in 1:length(node)){
    s0_train <- sum(tab_train[part_train == i,1])
    s1_train <- sum(tab_train[part_train == i,2])
    if(s0_train >= s1_train)
      error <- error + sum(tab_val[part_val == i,2])
    else
      error <- error + sum(tab_val[part_val == i,1])
  }

  return(error/total)
}

#Get hypothesis
get_hypothesis <- function(tab,node){
  part <- vector(length = nrow(tab))
  for(i in 1:length(node))
    part[substr(rownames(tab),1,nchar(node[[i]])) == node[[i]]] <- i

  h <- data.frame("x" = NA,"hx" = NA)
  for(i in 1:length(node)){
    s0 <- sum(tab[part == i,1])
    s1 <- sum(tab[part == i,2])
    if(s0 == s1)
      h <- rbind.data.frame(h,data.frame("x" = node[[i]],"hx" = "Tie"))
    else if(s0 >= s1)
      h <- rbind.data.frame(h,data.frame("x" = node[[i]],"hx" = "0"))
    else
      h <- rbind.data.frame(h,data.frame("x" = node[[i]],"hx" = "1"))
  }

  h <- h[order(h$x),]
  return(na.omit(h))
}

#Test if string starts with some string
string_start <- function(string,start){
  return(start == substr(string,1,nchar(start)))
}

#Clean node
clean_node <- function(node,tab_train){
  new_node <- vector()
  node <- unlist(node)
  node <- node[order(nchar(node),decreasing = F)]
  s <- length(node) + 1
  while(length(node) < s){
    h <- get_hypothesis(tab = tab_train,node = node)
    s <- length(node)
    while(length(node) > 0){
      if(nchar(node[1]) > 1){
        x <- node[c(T,substr(node[-1],1,nchar(node[1])-1) == substr(node[1],1,nchar(node[1])-1))]
        y <- unique(h$hx[h$x %in% x])
        if(length(y) > 1 || length(x) == 1)
          new_node <- c(new_node,x)
        else
          new_node <- c(new_node,substr(node[1],1,nchar(node[1])-1))
        node <- node[!(node %in% x)]
        if(length(node) == 1){
          new_node <- c(new_node,node[1])
          node <- vector()
        }
      }
      else{
        new_node <- c(new_node,node[1])
        node <- node[-1]
        if(length(node) == 1){
          new_node <- c(new_node,node[1])
          node <- vector()
        }
      }
    }
    node <- new_node[order(nchar(new_node),decreasing = F)]
    new_node <- vector()
  }
  return(as.list(node))
}

#Is.valid node
is.valid_node <- function(node){
  return(max(unlist(lapply(node,function(x) sum(x == unlist(lapply(node,function(y) substr(y,1,nchar(x)))))))) == 1)
}


