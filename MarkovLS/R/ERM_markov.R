#' @export
#' @title ERM estimative of Markov model of length k
#'
#' @description Estimate the hypothesis that minimizes the empirical error under the Markov model of length k.
#'
#' @details
#'
#' @param train A data frame with the training data. Should have collumns named x and y representing the input and output strings.
#' @param k Length of Markov Chain.
#'
#' @return A list containing the table with the estimated hypothesis and a function to predict the output of a given input.

ERM_markov <- function(train,k){
  #ERM hypothesis
  tab <- table(factor(substr(train$x,1,k)),train$y)
  hhat <- data.frame("x" = rownames(tab),"hx" = unlist(apply(as.matrix(tab),1,function(x) ifelse(min(x) != max(x),
                                                                                                 c(0,1)[x == max(x)],
                                                                                                 "Tie"))))

  #Function to predict
  pred <- function(x){
    x <- substr(x,1,k)
    hx <- hhat$hx[hhat$x == x]
    if(length(x) == 0)
      hx <- "0"
    else if(hx == "Tie")
      hx <- "0"
    return(hx)
  }

  return(list("hhat" = hhat,"predict" = pred))
}
