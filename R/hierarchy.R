#' @title Hierarchy Structure
#' 
#' @description Stores variable indexes of clustering hierarchies in a fast 
#' accessible manner.
#' 
#' @param x A S3 object.
#' @param max.height Is the maximal height below the global node height which 
#' is considered.
#' @param height A vector of heights at which nodes are grouped.
#' @param names Variable names in the order in which the indexes shut be given 
#' to the variables.
#' @param ... Further arguments.
#' 
#' @details For the HIT algorithm it is important to have the hierarchical 
#' clustering structure in a fast accessible format. This is provided by the 
#' hierarchy object generated with this function.
#' 
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
#' hc <- hclust(dist(t(x)))
#' hier <- as.hierarchy(hc, max.height = 20)
#' 
#' @export
as.hierarchy <- function(x, max.height, height, names) 
  UseMethod("as.hierarchy")


#' @export
as.hierarchy.hierarchy <- function (x, max.height, height, names, ...) x


#' @export
as.hierarchy.hclust <- function (x, max.height, height, names, ...) {
    as.hierarchy(as.dendrogram(x), max.height, height, names, ...)
}


#' @export
as.hierarchy.dendrogram <- function (x, max.height, height, names, ...) {
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
  out <- unname(dend2hier(x, as.numeric(height), as.character(names)))
  ordAll <- order(out[[1L]])
  out[[1L]][] <- out[[1L]][ordAll]
  names(out[[1L]]) <- names(out[[1L]])[ordAll]
  class(out) <- "hierarchy"
  out
}


#' @title Heights of Dendrogram
#' 
#' @description All heights from a dendrogram. 
#' 
#' @param x a \code{\link[stats]{dendrogram}}.
#' 
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
#' 
#' @description All heights from a hierarchy.
#' 
#' @param x a \code{\link{as.hierarchy}}.
#' 
#' @keywords internal
heightHierarchy <- function(x) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  sort(unique(sapply(x, attr, which = "height")))
}


#' @title Names of Hierarchy
#' 
#' @description Names of variables of an hierarchy.
#' 
#' @param x a \code{\link{as.hierarchy}}.
#' 
#' @export
names.hierarchy <- function(x) {
  names(x[[1L]])
}


#' @title Reorder Hierarchy
#' 
#' @description Reorder indexes according to a vector of names.
#' 
#' @param x a \code{\link{as.hierarchy}}.
#' @param names variable names in the order in which the indexes shut be given 
#' to the variables.
#' @param ... further arguments passed to or from other methods (not used).
#' 
#' @importFrom stats reorder
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


# #' @title Leaf of Hierarchy 
# #'
# #' @description All the leafs of the hierarchy.
# #' @param x a \code{\link{as.hierarchy}}.
# #'
# #' @keywords internal
# bottomNodeIndex <- function(x) {
#   if (!inherits(x, "hierarchy")) 
#     stop("'x' is not a hierarchy")
#   which(sapply(x, function(x) is.null(attr(x, which = "subset"))))
# }


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
