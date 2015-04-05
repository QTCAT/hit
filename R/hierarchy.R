
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
    newLevel <- sum(attr(subtree, "height") <= height)
    if (is.leaf(subtree) && newLevel > level) {
      CLUSTERS[[COUNTER]] <<- labels(subtree)
      attr(CLUSTERS[[COUNTER]], "superset") <<- superset
      attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
      subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
      attr(CLUSTERS[[superset]], "subset") <<- subset      
      COUNTER <<- COUNTER + 1L
    } else {
      if (newLevel > level) {
        CLUSTERS[[COUNTER]] <<- labels(subtree)
        attr(CLUSTERS[[COUNTER]], "superset") <<- superset
        attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
        subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
        attr(CLUSTERS[[superset]], "subset") <<- subset
        superset <- COUNTER
        COUNTER <<- COUNTER + 1L
      }
      lapply(subtree, make.hierarchy, level = newLevel, superset = superset)
    }
    return(NULL)
  }
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  if (is.null(height)) 
    height <- heightLevels(x)
  height <- sort(height, decreasing = TRUE)
  CLUSTERS <- list(labels(x))
  attr(CLUSTERS[[1L]], "height") <- height[1L]
  COUNTER <- 2L
  lapply(x, make.hierarchy, level = 1L, superset = 1L)
  class(CLUSTERS) <- "hierarchy"
  CLUSTERS
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

#' @title root of thr hierarchy 
#' @param x a hierarchy
#' @keywords internal
topSet <- function(x) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  which(sapply(x, function(x) is.null(attr(x, which = "superset"))))
}

#' @title leaf of thr hierarchy 
#' @param x a hierarchy
#' @keywords internal
bottomSet <- function(x) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  which(sapply(x, function(x) is.null(attr(x, which = "subset"))))
}