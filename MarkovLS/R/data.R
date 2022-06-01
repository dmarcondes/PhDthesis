#' @title Bitcoin history
#'
#' @description Data set containing the bitcoin value history from April 30th 2013 to April 6th 2022.
#' @details The data is divided into training, validation and test samples.
#' @return \item{date}{The date}
#' @return \item{open}{Open bitcoin value in US dollars in the respective date}
#' @return \item{close}{Close bitcoin value in US dollars in the respective date}
#' @return \item{variation}{Percentage variation from open to close value}
#' @return \item{y}{Variation positive (1) or negative (0)}
#' @return \item{sample}{If the respetive date is part of the training, validation or test sample}
#' @return \item{x}{The past variation (positive or negative) on the prior 30 days}
#' @references Diego Marcondes. A data-driven systematic, consistent and feasible approach to Model Selection. \emph{PhD thesis}. Universidade de São Paulo, 2022
"bitcoin"
