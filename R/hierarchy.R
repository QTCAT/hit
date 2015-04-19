
#' @title Hierachy indeces from desndrogram
#' @param x object of class 'dendrogram'.
#' @param height vector of heghts at which nodes are grouped.
#' @param max.height is the maximal heigt at which the testing starts.
#' @param names names in order.
#' @author Jonas Klasen
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 150
#' # x with correlated columns
#' corMat <- toeplitz(p:1/p)
#' corMatQ <- chol(corMat)
#' x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
#' colnames(x) <- paste0("x", 1:p)
#' # hierarchy
#' dend <- as.dendrogram(hclust(dist(t(x))))
#' hier <- hierarchy(dend, max.height = 20)
#' @importFrom parallel mclapply
#' @export
hierarchy <- function (x, height, max.height, names) {
  make.hierarchy <- function(subtree, level, superset) {
    newLevel <- sum(attr(subtree, "height") <= height)
    if (is.leaf(subtree) && newLevel > level) {
      CLUSTERS[[COUNTER]] <<- match(labels(subtree), names)
      attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
      attr(CLUSTERS[[COUNTER]], "superset") <<- superset
      subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
      attr(CLUSTERS[[superset]], "subset") <<- subset      
      COUNTER <<- COUNTER + 1L
    } else {
      if (newLevel > level) {
        CLUSTERS[[COUNTER]] <<- match(labels(subtree), names)
        attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
        attr(CLUSTERS[[COUNTER]], "superset") <<- superset
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
  if (missing(height)) 
    height <- heightDendrogram(x)
  height <- sort(height, decreasing = TRUE)
  if (missing(max.height))
    max.height <- attr(x, "height")
  height <- height[height <= max.height]
  if (missing(names)) {
    names <- labels(x)
  } else if (length(setdiff(labels(x), names))) {
    stop("'x' includs variabels not in 'names'")
  }
  CLUSTERS <- list(seq_along(names))
  attr(CLUSTERS[[1L]], "names") <- names
  attr(CLUSTERS[[1L]], "height") <- height[1L]
  COUNTER <- 2L
  lapply(x, make.hierarchy, level = 1L, superset = 1L)
  class(CLUSTERS) <- "hierarchy"
  CLUSTERS
}

#' @title All heights from a dendrogram 
#' @param x a dendrogram
#' @keywords internal
heightDendrogram <- function (x) {
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

#' @title All heights from a hierarchy 
#' @param x a hierarchy
#' @keywords internal
heightHierarchy <- function(x) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  sort(unique(sapply(x, attr, which = "height")))
}

# #' @title Leaf of the hierarchy 
# #' @param x a hierarchy
# #' @keywords internal
# bottomNodeIndex <- function(x) {
#   if (!inherits(x, "hierarchy")) 
#     stop("'x' is not a hierarchy")
#   which(sapply(x, function(x) is.null(attr(x, which = "subset"))))
# }

#' @title names of variables in hierarchy
#' @param x a hierarchy.
#' @method names hierarchy 
#' @export
names.hierarchy <- function(x) {
  names(x[[1]])
}

#' @title reorder hierarchy according to names vector
#' @param x a hierarchy.
#' @param names names in new order.
#' @param ... further arguments passed to or from other methods (not used).
#' @importFrom stats reorder
#' @method reorder hierarchy 
#' @export
reorder.hierarchy <- function(x, names, ...) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  if (length(setdiff(names(x[[1]]), names)))
    stop("'x' includs variabels not in 'names'")
  newOrder <- match(names(x[[1]]), names)
  out <- lapply(x, function(x, newOrder){x[] <- sort(newOrder[x]); x}, newOrder)
  names(out[[1]]) <- names[out[[1]]]
  class(out) <- "hierarchy"
  out
}
