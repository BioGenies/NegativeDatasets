#' Create grey model GM(1, 1)
#' 
#' @description Builds grey model with parameters 1 and 1
#' from the supplied data. Returns a numeric vector of
#' length 2 with \code{a} and \code{b} coefficients for 
#' supplied observation.
#' 
#' @param row
#'  Training data in numeric vector format without target
#'  value. Should contain data about one sequence.
#' 
#' @examples
#' X <- list(
#'   X1 = 1:5,
#'   X2 = c(10, 12, 24, 60),
#'   X3 = 9:2
#' )
#' lapply(X, grey_model_1_1)
#' 
#' @noRd
grey_model_1_1 <- function(row) {
  # First we compute first-order AGO
  ago <- cumsum(row)
  # Unfortunately, had to write for-loop, because of trying to
  #  iterate over rows of two matrices simultaneously.
  # -a is called "developing coefficient"
  # b is called "influence coefficient"
  # [a, b]^T = [B^T*B]^-1*B^T*Y
  B <- cbind(-0.5*(ago[-length(ago)] + ago[-1]), 1)
  Y <- row[-1]
  # ^ has higher precedence than %*%, so no parentheses necessary
  ret <- as.numeric((t(B) %*% B)^-1 %*% t(B) %*% Y)
  setNames(ret, c("a", "b"))
}
