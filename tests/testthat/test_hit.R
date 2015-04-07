context("HIT testing")

test_that("hit testing", {
  set.seed(123)
  n <- 100
  p <- 150
  # x with correlated columns
  sigmaChol <- chol(toeplitz((p:1/p)^3), pivot = TRUE)
  sigmaChol <- sigmaChol[, order(attr(sigmaChol, "pivot"))]
  x <- matrix(rnorm(n * p), nrow = n) %*% sigmaChol
  x <- sweep(x, 2, rnorm(p, 10), "+")
  colnames(x) <- paste0("x", 1:p)
  # y
  y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
  # hierarchy
  dend <- as.dendrogram(hclust(dist(t(x))))
  hier <- hierarchy(dend)
  # HIT
  fit <- hit(x, y, hier)
  # checks
  expect_equal(class(fit), "hit")
  expect_equal(unname(unlist(lapply(fit, class))), 
               c("numeric", "numeric","hierarchy", "character"))
#   expect_equal(names(fit), 
#                c("pValues", "selectFreq", "hierarchy", "additionals"))
})
