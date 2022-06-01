#' @export
#' @title ERM estimative of Markov model of length k
#'
#' @description Estimate the hypothesis that minimizes the empirical classification error under the Markov model of length k.
#'
#' @details Receives the sequence values y and its past x, and returns the predictor of length k that minimizes the classification error in the sample x,y.
#'
#' @param train A data frame with the training data. Should have columns named x and y representing the input and output strings.
#' @param k Length of Markov Chain.
#'
#' @return A list containing the table with the estimated hypothesis and a function to predict the output of a given input.
#'
#' @examples
#' ERM_markov(bitcoin,2)
#'
#' @references Diego Marcondes. A data-driven systematic, consistent and feasible approach to Model Selection. \emph{PhD thesis}. Universidade de São Paulo, 2022

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
