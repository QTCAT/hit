context("hirarchy testing")

test_that("hirarchy testing", {
  set.seed(123)
  n <- 100
  p <- 150
  # x with correlated columns
  corMat <- toeplitz(p:1/p)
  corMatQ <- chol(corMat)
  x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
  colnames(x) <- paste0("x", 1:p)
  # y
  y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
  # hierarchy
  dend <- as.dendrogram(hclust(dist(t(x))))
  hier <- hierarchy(dend, max.height = 20)
  # check:
  expect_equal(class(hier), "hierarchy")
#   expect_equal(unname(unlist(lapply(hier, class))), 
#                c("list", "list", "character", "numeric"))
#   expect_equal(names(hier), 
#                c("hierarchyCluster", "clusterMembers", 
#                  "labels", "hierarchyLevel"))
})
