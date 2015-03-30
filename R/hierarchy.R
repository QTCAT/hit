
#' @title Hierachy indeces from desndrogram
#' @param dendrogram object of class 'dendrogram'.
#' @param min.dist min. distence to be tested.
#' @param max.dist max. distence to be tested.
#' @param h.levels number of hierarchy levels to be tested.
#' @param h.cut vector of cutting points. If it is specified min.dist, 
#' max.dist, and h.levels are ignored.
#' @param mc.cores number of cores for parallelising, 
#' see \code{\link{mclapply}}.
#' @author Jonas Klasen
#' @examples
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
hierarchy <- function (dendrogram, min.dist=0, max.dist=max(dist), 
                       h.levels=100L, h.cut=NULL, mc.cores=1L) {
  stopifnot(class(dendrogram) == "dendrogram")
  dend.names <- labels(dendrogram)
  p <- length(dend.names)
  if (is.null(h.cut)) {
    dist <- dend.heights(dendrogram)
    cut.min <- which.min(abs(dist-min.dist))
    cut.max <- which.min(abs(dist-max.dist))
    if ((cut.min-cut.max) > h.levels) {
      h.cut <- dist[seq(cut.max, cut.min, length.out=h.levels)]
    } else {
      h.cut <- dist[cut.max:cut.min]
    }
  } else {
    h.cut <- sort(h.cut, decreasing=TRUE)
  } # if (is.null(h.cut))
  clust <- function(h.cut, dendrogram) {    
    cluster <- clust.h(dendrogram, h.cut)
    cluster
  } # clust()
  cluster <- mclapply(h.cut, clust, dendrogram,
                      mc.cores=mc.cores, mc.cleanup=TRUE)
  if (length(unique(na.exclude(cluster[[1]]))) > 1) {
    h.cut <- c(Inf, h.cut)
    rone <- rep(1L, p)
    cluster <- c(list(rone), cluster)
  } 
  cluster <- do.call("rbind", cluster)
  hierarchy <- clusterInx(cluster)
  out <- list("hierarchyCluster"=hierarchy$hierarchyCluster,
              "clusterMembers"=hierarchy$clusterMembers,
              "labels"=dend.names,
              "hierarchyLevel"=h.cut)
  class(out) <- "hierarchy"
  out
} # hierarchy

#' @title all heights from a dendrogram 
#' @param x a dendrogram
#' @keywords internal
dend.heights <- function(x) {
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  node.height <- function(d) {
    if (is.list(d)) {
      r <- attributes(d)$height
      return(c(r, node.height(d[[1L]]), node.height(d[[2L]])))
    }
    attributes(d)$height
  } # node.height()
  sort(unique(round(node.height(x), 8L)), decreasing=TRUE)
} # dend.heights

#' @title cluster dendrogram
#' @param x dendrogram
#' @param h cutting hight 
#' @keywords internal
clust.h <- function(x, h) {
  stopifnot(!missing(h))
  if (h >= attr(x, "height")) {
    names.clust <- labels(x)
    cluster <- rep(1, length(names.clust))
    names(cluster) <- names.clust
  } else {
    cut.x <- cut(x, h=h)$lower
    clust.member <- function(i, x) {
      names.clust <- labels(x[[i]])
      clust <- rep(i, length(names.clust))
      names(clust) <- names.clust
      return(clust)
    }
    cluster <- lapply(1:length(cut.x), clust.member, cut.x)
  }  
  unlist(cluster)
} # clust.h
