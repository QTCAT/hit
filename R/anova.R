#' @title Fast ANOVA
#' 
#' @description A fast sequential analysis of variance (ANOVA). Mainly 
#' developed for internal use.
#' 
#' @param x Design matrix of dimension \code{n * p}.
#' @param y Response vector of observations of length \code{n}.
#' @param assign Integer vector assigning columns to terms can be also given as 
#' \code{x} attribute in which case the argument is ignored. If an intercept 
#' exist it is expected to be coded by '0'.For details about assign see 
#' \code{\link[stats]{model.matrix}}.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{anova}}, and 
#' \code{\link[stats]{aov}}.
#' 
#' @examples
#' y <- rnorm(n = 100)
#' x <- matrix(data = rnorm(1000), nrow = 100)
#' a <- 1:10
#' fast.anova(x = x, y = y, assign = a)
#' 
#' @importFrom stats lm.fit pf
#' @export
fast.anova <- function(x, y, assign = NULL) {
  if (is.null(assign) && !is.null(attr(x, "assign"))) {
    assign <- attr(x, "assign")
  }
  if (is.null(assign)) {
    stop("either a 'x' attribute 'assign' or the 'assign' argument must be specified")
  }
  stopifnot(ncol(x) == length(assign))
  stopifnot(nrow(x) == length(y))
  # LM fit by pivoted QR decomposition
  fit <- lm.fit(x, y)
  if (assign[1L] == 0L) { 
    # with intercept
    full.rank <- 1L:(fit$rank - 1L)
    assign.pivot <- assign[fit$qr$pivot[full.rank + 1L]]
    var <- fit$effects[-1L]^2
  } else { 
    # without intercept
    full.rank <- 1L:fit$rank
    assign.pivot <- assign[fit$qr$pivot[full.rank]]
    var <- fit$effects^2
  }
  # Treatment: Sum Sq | Df | Mean Sq
  ss.treat <- tapply(var[full.rank], assign.pivot, "sum")
  df.treat <- table(assign.pivot)
  ms.treat <- ss.treat / df.treat
  # Residuals: Sum Sq | Df | Mean Sq
  ss.res <- sum(var[-full.rank])
  df.res <- fit$df.residual
  ms.res <- ss.res / df.res
  # F value
  f <- ms.treat / ms.res 
  # p value
  p <- rep(1, max(assign))
  p[unique(assign.pivot)] <- pf(f, df.treat, df.res, lower.tail = FALSE)
  p
}
