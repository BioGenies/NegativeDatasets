#' Create grey model GM(1, 1)
#' 
#' @description Builds grey model with parameters 1 and 1
#' from the supplied data. Returns a matrix with \code{a}
#' and \code{b} coefficients for each observation.
#' 
#' @param data
#'  Training data in numeric matrix format without target
#'  column. Each row contains data about one sequence.
#' 
#' @examples
#' X <- matrix(1:15, nrow = 5,
#'             dimnames = list(LETTERS[11:15], c("V1", "V2", "V3")))
#' grey_model_1_1(X)
#' 
#' @noRd
grey_model_1_1 <- function(data) {
  # First we compute first-order AGO
  ago <- t(apply(data, 1, cumsum))
  # Unfortunately, had to write for-loop, because of trying to
  #  iterate over rows of two matrices simultaneously.
  # -a is called "developing coefficient"
  # b is called "influence coefficient"
  # [a, b]^T = [B^T*B]^-1*B^T*Y
  ret <- matrix(nrow = nrow(data), ncol = 2,
                dimnames = list(rownames(data), c("a", "b")))
  for (i in 1:nrow(data)) {
    B <- cbind(-0.5*(ago[i, -ncol(ago)] + ago[i, -1]),
               1)
    Y <- data[i, -1]
    # ^ has higher precedence than %*%, so no parentheses necessary
    ret[i, ] <- (t(B) %*% B)^-1 %*% t(B) %*% Y
  }
  ret
}
