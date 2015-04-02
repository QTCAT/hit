
#' @title Hierarchical inference testing
#' @description Hierarchical inference testing for linear models with 
#' high-dimensional and/or correlated covariates by repeated sample splitting.
#' @param x design matrix of dimension \code{n * p}, without intercept.
#' Variables not part of the dendrogram are added to the HO-model, see Details 
#' below.
#' @param y quantitative response variable dimension \code{n}.
#' @param hierarchy object of class \code{\link{hierarchy}}. Must include all 
#' variables of \code{x} which should be trested.
#' @param B number of sample-splits.
#' @param p.samp1 fraction of data used for the LASSO. The ANOVA uses 
#' \code{1 - p.samp1}.
#' @param gamma vector of gamma-values.
#' @param max.p.esti maximum alpha level. All p-values above this value are set 
#' to one. Small \code{max.p.esti} values reduce computing time.
#' @param mc.cores number of cores for parallelising. Theoretical maximum is 
#' 'B'. For details see \code{\link[parallel]{mclapply}}.
#' @param trace if TRUE it prints current status of the program.
#' @param ... additional arguments for \code{\link[glmnet]{cv.glmnet}}.
#' @details The H0-model contains variables, with are not tested, like 
#' experimental-design variables. These variables are not penalised in the 
#' LASSO model selection and are always include in the reduced ANOVA model.
#' @author Jonas Klasen
#' @references
#'   Mandozzi, J. and Buehlmann, P. (2013). \emph{Hierarchical testing in the 
#'   high-dimensional setting with correlated variables.} To appear in the 
#'   Journal of the American Statistical Association. Preprint arXiv:1312.5556
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
#' # y
#' y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
#' # hierarchy
#' dend <- as.dendrogram(hclust(dist(t(x))))
#' hier <- hierarchy(dend)
#' # HIT
#' out <- hit(x, y, hier)
#' summary(out)
#' @importFrom parallel mclapply
#' @export 
hit <- function(x, y, hierarchy, B=50, p.samp1=0.5, 
                gamma=seq(0.05, 0.99, by=0.01), max.p.esti=1, 
                mc.cores=1L, trace=FALSE, ...) {
  ## Mandozzi and Buehlmann, 2013
  # see chapter 2 Description of method
  if (is.null(colnames(x)))
    stop ("column names of 'x' are missing")
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(class(hierarchy) == "hierarchy")
  # checks order of variables
  x.names <- colnames(x)
  hier.names <- hierarchy$labels
  if (any(!(hier.names %in% x.names)))
    stop("'hierarchy' includs variabels not in 'x'")
  if (length(hier.names) == length(x.names) && all(hier.names == x.names)) {
    x.nonTested <- NULL
    additionals <- character(0L)
  } else if (setequal(hier.names, x.names)) {
    x <- x[, hier.names]
    x.nonTested <- NULL
    additionals <- character(0L)
  } else {
    x <- x[, union(hier.names, x.names)]
    x.nonTested <- which(!(colnames(x) %in% hier.names))
    additionals <- colnames(x)[x.nonTested]
  }
  gamma <- sort(gamma)
  if (p.samp1 < .1 | p.samp1 > .5)
    stop("'p.samp1' must be between .1 and 0.5")
  n.samp1 <- as.integer(n*p.samp1)
  n.samp2 <- n-n.samp1
  penalty.factor <- rep(1L, p)
  penalty.factor[x.nonTested] <- 0L
  ## subsamples
  allSamp1.ids <- replicate(B, sample.int(n, n.samp1), simplify=FALSE)
  ## 2.2 Screening
  if (isTRUE(trace))
    cat("LASSO has started at:\n\t", as.character(Sys.time()), "\n")
  allActSet.ids <- mclapply(allSamp1.ids, samp1.lasso,
                            x, y, n.samp2, penalty.factor, ...,
                            mc.cores=mc.cores, mc.cleanup=TRUE)
  ## 2.3 Testing and multiplicity adjustmen; and 
  ## 2.4 Aggregating and Hierarchical adjustment
  if (isTRUE(trace))
    cat("Significance testing has started at:\n\t", 
        as.character(Sys.time()), "\n\t Loop:\n")
  p.values <- rep(NA_real_, length(hierarchy$clusterMembers))
  for (level in seq_along(hierarchy$hierarchyCluster)) {
    if (isTRUE(trace))
      cat(level, " ")
    clust.ids <- hierarchy$hierarchyCluster[[level]]
    if (level > 1L) {
      max.upper <- max(hierarchy$hierarchyCluster[[level-1L]])
      clust.ids <- clust.ids[clust.ids > max.upper]
    }
    p.values[clust.ids] <- unlist(mclapply(clust.ids, samp2.testing, 
                                           x, y, p.values, hierarchy, 
                                           level, x.nonTested,
                                           allSamp1.ids, allActSet.ids, 
                                           max.p.esti, B, gamma, 
                                           mc.preschedule=TRUE, 
                                           mc.cores=mc.cores, 
                                           mc.cleanup=TRUE))
  }
  # make output
  if (isTRUE(trace))
    cat("\nHIT has finished at:\n\t", as.character(Sys.time()), "\n")
  asi <- sort(unlist(allActSet.ids))
  sel.tab <- rep(0, p)
  sel.tab[unique(asi)] <- table(asi)/B
  # Output
  out <- list("pValues"=p.values,
              "selectFreq"=sel.tab,
              "hierarchy"=hierarchy,
              "additionals"=additionals)
  class(out) <- "hit"
  out
} # hit

#' @title LASSO screening
#' @param samp1 list of index for subsample (mclapply index).
#' @param x design matrix, of dimension n x p.
#' @param y vector of quantitative response variable.
#' @param n.samp2 number of individuals in samp2 which is the max. 
#' for non zero coefficients.
#' @param penalty.factor see glmnet.
#' @param ... aditional agruments
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @keywords internal
samp1.lasso <- function (samp1, x, y, n.samp2, penalty.factor, ...) {
  ## Mandozzi and Buehlmann, 2013
  # see chapter 2 Description of method
  ## 2.2 Screening
  x <- x[samp1, ]
  y <- y[samp1]
  lasso.fit <- cv.glmnet(x, y, penalty.factor=penalty.factor, 
                         dfmax = n.samp2-2L, ...)
  beta <- coef(lasso.fit)[-1L]
  actSet <- which(beta != 0 & penalty.factor == 1L)
  return(actSet)
} # samp1.lasso

#' @title ANOVA testing, multiplicity adjustment, aggregating and hierarchical 
#' adjustment
#' @param clust.id index for cluster-ID (mclapply index).
#' @param x design matrix, of dimension n x p.
#' @param y vector of quantitative response variable.
#' @param p.values vector of p value per clusters.
#' @param hierarchy a hierarchy object
#' @param level hierarchy level
#' @param x.nonTested  vector of indeces of non tested variabels.
#' @param allSamp1.ids  list of subsampels.
#' @param allActSet.ids list of active sets.
#' @param max.p.esti maximum alpha level. All p-values above this value are set 
#' to one. Small max.p.esti values reduce computing time.
#' @param B number of sample-splits.
#' @param gamma vector of gamma-values.
#' @keywords internal
samp2.testing <- function(clust.id, x, y, p.values, hierarchy, level, 
                          x.nonTested, allSamp1.ids, allActSet.ids, 
                          max.p.esti, B, gamma) { 
  clust.inx <- hierarchy$clusterMembers[[clust.id]]
  if (level == 1) {
    upper.p <- 0
  } else {
    upper.level <- hierarchy$hierarchyCluster[[level-1L]]
    upper.inx <- hierarchy$clusterMembers[upper.level]
    upper.hier <- sapply(upper.inx, function(i, j) any(j %in% i), clust.inx)
    upper.clust <- upper.level[upper.hier]
    upper.p <- max(p.values[upper.clust])
  }
  if (upper.p >= max.p.esti) {
    p.value <- 1
  } else {
    ## 2.3 Testing and multiplicity adjustment
    p.clust <- sapply(1L:B, samp2.anova, 
                 x, y, clust.inx, x.nonTested, 
                 allSamp1.ids, allActSet.ids)
    ## 2.4-1 Aggregating 
    q.aggre <- sapply(gamma, 
                      function(i, x) { min(1, quantile(x/i, i)) }, 
                      x = p.clust)
    p.aggre <- min(1, (1-log(min(gamma)))*min(q.aggre))
    ## 2.4-2 Hierarchical adjustment
    p.value <- max(p.aggre, upper.p)
  }
  p.value
} # samp2.testing

#' @title ANOVA testing and multiplicity adjustment
#' @param k index for subsample (mclapply index).
#' @param x design matrix, of dimension n x p.
#' @param y vector of quantitative response variable.
#' @param cluster clusters to be tested.
#' @param x.nonTested vector of indeces of non tested variabels.
#' @param allSamp1.ids  list of subsampels.
#' @param allActSet.ids list of active sets.
#' @keywords internal
samp2.anova <- function (k, x, y, cluster, x.nonTested, 
                         allSamp1.ids, allActSet.ids) {
  n <- nrow(x)
  ## Mandozzi and Buehlmann, 2013
  # see chapter 2 Description of method
  actClust <- intersect(allActSet.ids[[k]], cluster)
  nonActClust <- setdiff(allActSet.ids[[k]], cluster)
  ## 2.3-1 Testing
  p.cluster <- 1
  if (l.actClust <- length(actClust)) {
    # ANOVA between active set and active set minus cluster
    samp2 <- (1L:n)[-allSamp1.ids[[k]]]
    y <- y[samp2]
    if (l.nonActClust <- length(nonActClust)) {
      if (l.nonTested <- length(x.nonTested)) {
        x  <- cbind(1, x[samp2, c(x.nonTested, nonActClust, actClust)])
        assign <- rep(0L:3L, c(1L, l.nonTested, l.nonActClust, l.actClust))
        get.p <- 3L
      } else {
        x  <- cbind(1, x[samp2, c(nonActClust, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonActClust, l.actClust))
        get.p <- 2L
      }
    } else {
      if (l.nonTested <- length(x.nonTested)) {
        x  <- cbind(1, x[samp2, c(x.nonTested, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonTested, l.actClust))
        get.p <- 2L
      } else {
        x  <- cbind(1, x[samp2, actClust])
        assign <- rep(0L:1L, c(1L, l.actClust))
        get.p <- 1L
      }
    } # if (l.nonActClust <- length(nonActClust))
    p.test <- fast.anova(x, y, assign)[get.p]
    ## 2.3-2 Multiplicity adjustment
    p.cluster <- min(1, p.test*(l.actClust+l.nonActClust)/l.actClust)
  } # if (l.actClust <- length(actClust))
  return(p.cluster) 
} # samp2.anova

#' Significant clusters at alpha threshold
#'
#' @param object a hit object.
#' @param alpha alpha level.
#' @param min.dist minimal distance for significant clusters.
#' @param max.dist maximal distance for significant clusters.
#' @param ... further arguments passed to or from other methods (not used).
#' @method summary hit
#' @export 
summary.hit <- function(object, alpha=.05, 
                        min.dist=0, max.dist=max(dist), ...) {
  stopifnot(is(object, "hit"))
  dist <- object$hierarchy$hierarchyLevel
  min.inx <- sum(dist >= min.dist)
  max.inx <- sum(dist > max.dist | dist == Inf)
  out <- sigCluster(object$pValues, 
                    object$hierarchy$hierarchyCluster, 
                    object$hierarchy$clusterMember, 
                    length(object$hierarchy$labels),
                    alpha, min.inx, max.inx)
  colnames(out) <- object$hierarchy$labels
  rownames(out) <- c("clusters", "pValues")
  out <- out[, out[1L, ] != 0, drop=FALSE]
  out
} # summary.hit

#' @title Significant hierarchy
#' @description Significant hierarchy
#' @param x a hit object
#' @details makes a matrix of p-values for image(hit:::p.hierarchy(x))
#' @keywords internal
p.hierarchy <- function (x) {
  out <- matrix(NA_real_, 
                length(x$hierarchy$hierarchyCluster), 
                length(x$hierarchy$labels))
  for (j in seq_along(x$hierarchy$hierarchyCluster)) {
    for (i in x$hierarchy$hierarchyCluster[[j]]) {
      out[j, x$hierarchy$clusterMember[[i]]] <- x$pValues[i]
    }
  }
  colnames(out) <- x$hierarchy$labels
  rownames(out) <- x$hierarchy$hierarchyLevel
  out
} # p.hierarchy
