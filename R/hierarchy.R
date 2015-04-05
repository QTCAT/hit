
#' @title Hierachy indeces from desndrogram
#' @param x object of class 'dendrogram'.
#' @param height vector of cutting points.
#' @author Jonas Klasen
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 150
#' # x with correlated columns
#' sigmaChol <- chol(toeplitz((p:1/p)^3), pivot = TRUE)
#' sigmaChol <- sigmaChol[, order(attr(sigmaChol, "pivot"))]
#' x <- matrix(rnorm(n * p), nrow = n) %*% sigmaChol
#' x <- sweep(x, 2, rnorm(p, 10), "+")
#' colnames(x) <- paste0("x", 1:p)
#' # hierarchy
#' dend <- as.dendrogram(hclust(dist(t(x))))
#' hier <- hierarchy(dend)
#' @importFrom parallel mclapply
#' @export
hierarchy <- function (x, height = NULL) {
  make.hierarchy <- function(subtree, level, superset) {
    h <- attr(subtree, "height")
    newLevel <- length(height) - findInterval(h, vec = hInterval)
    if (is.leaf(subtree) && newLevel > level) {
      CLUSTERS[[COUNTER]] <<- match(labels(subtree), varLabel)
      attr(CLUSTERS[[COUNTER]], "superset") <<- superset
      subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
      attr(CLUSTERS[[superset]], "subset") <<- subset
      LEVELS[[newLevel]] <<- c(LEVELS[[newLevel]], COUNTER)
      COUNTER <<- COUNTER + 1L
    } else {
      if (newLevel > level) {
        CLUSTERS[[COUNTER]] <<- match(labels(subtree), varLabel)
        attr(CLUSTERS[[COUNTER]], "superset") <<- superset
        subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
        attr(CLUSTERS[[superset]], "subset") <<- subset
        LEVELS[[newLevel]] <<- c(LEVELS[[newLevel]], COUNTER)
        superset <- COUNTER
        COUNTER <<- COUNTER + 1L
      }
      lapply(subtree, make.hierarchy, level = newLevel, 
             superset = superset)
    }
    return(NULL)
  }
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  varLabel <- labels(x)
  if (is.null(height)) 
    height <- heightLevels(x)
  height <- sort(height)
  hInterval <- height[-1] - diff(height)/2
  CLUSTERS <- list(seq_along(varLabel))
  attr(CLUSTERS[[1]], "superset") <- 0L
  LEVELS <- rep(list(c()), length(height))
  LEVELS[[1L]] <- 1L
  COUNTER <- 2L
  lapply(x, make.hierarchy, level = 1L, superset = 1L)
  usedLevels <- !sapply(LEVELS, is.null)
  out <- list(hierarchyCluster = LEVELS[usedLevels], 
              clusterMembers = CLUSTERS, 
              labels = varLabel, 
              hierarchyLevel = sort(height,  decreasing = TRUE)[usedLevels])
  class(out) <- "hierarchy"
  out
}

#' @title all heights from a dendrogram 
#' @param x a dendrogram
#' @keywords internal
heightLevels <- function (x) {
  node.height <- function(d) {
    if (is.list(d)) {
      r <- attributes(d)$height
      return(c(r, node.height(d[[1L]]), node.height(d[[2L]])))
    }
    attributes(d)$height
  }
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  sort(unique(node.height(x)), decreasing = TRUE)
}
