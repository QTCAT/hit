#' @title Clustering hierachy
#' @description Stores variable indexes of clustering hierarchies in a fast 
#' accessible manner.
#' @param x a \code{\link[stats]{dendrogram}}.
#' @param height vector of heights at which nodes are grouped.
#' @param max.height is the maximal height below the global node height which 
#' is considered.
#' @param names variable names in the order in which the indexes shut be given 
#' to the variables.
#' @author Jonas Klasen
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 150
#' # x with correlated columns
#' corMat <- toeplitz((p:1/p)^3)
#' corMatQ <- chol(corMat)
#' x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
#' colnames(x) <- paste0("x", 1:p)
#' # hierarchy
#' dend <- as.dendrogram(hclust(dist(t(x))))
#' hier <- hierarchy(dend, max.height = 20)
#' @importFrom parallel mclapply
#' @export
hierarchy <- function (x, height, max.height, names) {
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  if (missing(height)) 
    height <- heightDendrogram(x)
  height <- sort(height, decreasing = TRUE)
  if (missing(max.height))
    max.height <- attr(x, "height")
  height <- height[height <= max.height]
  if (attr(x, "height") > max.height)
    height <- c(attr(x, "height"), height)
  if (missing(names)) {
    names <- labels(x)
  } else if (length(setdiff(labels(x), names))) {
    stop("'x' includs variabels not in 'names'")
  }
  out <- unname(dend2hier(x, height, names))
  ordAll <- order(out[[1L]])
  out[[1L]][] <- out[[1L]][ordAll]
  names(out[[1L]]) <- names(out[[1L]])[ordAll]
  class(out) <- "hierarchy"
  out
}

#' @title Heights of dendrogram
#' @description All heights from a dendrogram. 
#' @param x a \code{\link[stats]{dendrogram}}.
#' @author Jonas Klasen
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

#' @title Heights of Hierarchy
#' @description All heights from a hierarchy.
#' @param x a \code{\link{hierarchy}}.
#' @author Jonas Klasen
#' @keywords internal
heightHierarchy <- function(x) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  sort(unique(sapply(x, attr, which = "height")))
}

# #' @title Leaf of hierarchy 
# #' @description All the leafs of the hierarchy.
# #' @param x a \code{\link{hierarchy}}.
# #' @keywords internal
# bottomNodeIndex <- function(x) {
#   if (!inherits(x, "hierarchy")) 
#     stop("'x' is not a hierarchy")
#   which(sapply(x, function(x) is.null(attr(x, which = "subset"))))
# }

#' @title Names  of hierarchy
#' @description Names of variables of an hierarchy.
#' @param x a \code{\link{hierarchy}}.
#' @author Jonas Klasen
#' @method names hierarchy 
#' @export
names.hierarchy <- function(x) {
  names(x[[1L]])
}

#' @title Reorder hierarchy
#' @description Reorder indexes according to a vector of names.
#' @param x a \code{\link{hierarchy}}.
#' @param names variable names in the order in which the indexes shut be given 
#' to the variables.
#' @param ... further arguments passed to or from other methods (not used).
#' @author Jonas Klasen
#' @importFrom stats reorder
#' @method reorder hierarchy 
#' @export
reorder.hierarchy <- function(x, names, ...) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  if (length(setdiff(names(x[[1L]]), names)))
    stop("'x' includs variabels not in 'names'")
  newOrder <- match(names(x[[1L]]), names)
  out <- lapply(x, function(x, newOrder) {
    x[] <- sort(newOrder[x]) 
    x
  }, newOrder)
  names(out[[1L]]) <- names[out[[1L]]]
  class(out) <- "hierarchy"
  out
}

# #' A pure R version of hierarchy
# #' @export
# hierarchy2 <- function (x, height, max.height, names) {
#   make.hierarchy <- function(subtree, level, superset) {
#     newLevel <- sum(attr(subtree, "height") <= height)
#     if (is.leaf(subtree) && newLevel > level) {
#       CLUSTERS[[COUNTER]] <<- sort(match(labels(subtree), names))
#       attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
#       attr(CLUSTERS[[COUNTER]], "superset") <<- superset
#       subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
#       attr(CLUSTERS[[superset]], "subset") <<- subset      
#       COUNTER <<- COUNTER + 1L
#     } else {
#       if (newLevel > level) {
#         CLUSTERS[[COUNTER]] <<- sort(match(labels(subtree), names)) 
#         attr(CLUSTERS[[COUNTER]], "height") <<- height[newLevel]
#         attr(CLUSTERS[[COUNTER]], "superset") <<- superset
#         subset <- c(attr(CLUSTERS[[superset]], "subset"), COUNTER)
#         attr(CLUSTERS[[superset]], "subset") <<- subset
#         superset <- COUNTER
#         COUNTER <<- COUNTER + 1L
#       }
#       lapply(subtree, make.hierarchy, level = newLevel, superset = superset)
#     }
#     return(NULL)
#   }
#   if (!inherits(x, "dendrogram")) 
#     stop("'x' is not a dendrogram")
#   if (missing(height)) 
#     height <- heightDendrogram(x)
#   height <- sort(height, decreasing = TRUE)
#   if (missing(max.height))
#     max.height <- attr(x, "height")
#   height <- height[height <= max.height]
#   if (attr(x, "height") > max.height)
#     height <- c(attr(x, "height"), height)
#   if (missing(names)) {
#     names <- labels(x)
#   } else if (length(setdiff(labels(x), names))) {
#     stop("'x' includs variabels not in 'names'")
#   }
#   CLUSTERS <- list(seq_along(names))
#   attr(CLUSTERS[[1L]], "names") <- names
#   attr(CLUSTERS[[1L]], "height") <- height[1L]
#   COUNTER <- 2L
#   lapply(x, make.hierarchy, level = 1L, superset = 1L)
#   class(CLUSTERS) <- "hierarchy"
#   CLUSTERS
# }
