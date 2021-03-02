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
#' # If we have a tibble with data and want to append columns
#' # created with grey model:
#' data <- tibble(Ala = c(3, 6, 2), Gly = c(3, 1, 3))
#' # Each element of a list is a vector of some statistic
#' # (e.g. hydrophobicity) for given sequence:
#' stats <- list(c(1, 6, 3, 7), c(2, 0, 0, 4, 0), 6:14)
#' # Transposition and vapply are the safest way to bind cols
#' # However, names have to be updated still
#' bind_cols(
#'   data,
#'   t(vapply(stats, grey_model_1_1, FUN.VALUE = numeric(2)))
#' )
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
