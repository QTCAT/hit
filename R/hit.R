
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
#' # summary(out)
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
  hier.names <- names(hierarchy[[1]])
  if (length(setdiff(hier.names, x.names)))
    stop("'hierarchy' includs variabels not in 'x'")
  if (identical(hier.names, x.names)) {
    additionalCovariates <- character(0L)
    x.nonTested <- NULL
  } else if (setequal(hier.names, x.names)) {
    hierarchy <- reorder(hierarchy, x.names)
    additionalCovariates <- character(0L)
    x.nonTested <- NULL
  } else {
    hierarchy <- reorder(hierarchy, x.names)
    additionalCovariates <- setdiff(x.names, names(hierarchy))
    x.nonTested <- match(additionalCovariates, x.names)
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
  max.allow.recursive <- as.integer(sqrt(mc.cores))
  pValues <- samp2.testing(1, 0, 0, x, y, allSamp1.ids, allActSet.ids, 
                            x.nonTested, hierarchy, max.p.esti, B, gamma, 
                            max.allow.recursive)
  # make output
  if (isTRUE(trace))
    cat("\nHIT has finished at:\n\t", as.character(Sys.time()), "\n")
  asi <- sort(unlist(allActSet.ids))
  sel.tab <- rep(0, p)
  sel.tab[unique(asi)] <- table(asi)/B
  # Output
  out <- list("pValues"=pValues,
              "selectFreq"=sel.tab,
              "hierarchy"=hierarchy,
              "additionalCovariates"=additionalCovariates)
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
#' @param cIndex index for cluster (mclapply index).
#' @param level hierarchy level counter for parallelism
#' @param upper.p p value upper huerarchy level of the clusters variables.
#' @param x design matrix, of dimension n x p.
#' @param y vector of quantitative response variable.
#' @param allSamp1.ids  list of subsampels.
#' @param allActSet.ids list of active sets.
#' @param x.nonTested  vector of indeces of non tested variabels.
#' @param hierarchy a hierarchy object
#' @param max.p.esti maximum alpha level. All p-values above this value are set 
#' to one. Small max.p.esti values reduce computing time.
#' @param B number of sample-splits.
#' @param gamma vector of gamma-values.
#' @param max.allow.recursive max. level of recursive parallelism
#' @keywords internal
samp2.testing <- function(cIndex, level, upper.p, x, y, allSamp1.ids, 
                          allActSet.ids, x.nonTested, hierarchy, max.p.esti, 
                          B, gamma, max.allow.recursive) {
  ## 2.3 Testing and multiplicity adjustment
  cluster <- hierarchy[[cIndex]]
  ## 2.3-1 Testing
  ## 2.3 Testing and multiplicity adjustment
  p.cluster <- sapply(1L:B, samp2.anova, 
                      x, y, cluster, x.nonTested, 
                      allSamp1.ids, allActSet.ids)
  ## 2.4-1 Aggregating 
  q.aggre <- sapply(gamma, 
                    function(i, x) { min(1, quantile(x/i, i)) }, 
                    x = p.cluster)
  p.aggre <- min(1, (1-log(min(gamma)))*min(q.aggre))
  ## 2.4-2 Hierarchical adjustment
  p.value <- max(p.aggre, upper.p)
  ## estimation at next lower level
  if (!is.null(cIndeces <- attr(cluster, "subset"))) {
    if (p.value < max.p.esti) {
      if (level <= max.allow.recursive) {
        mc.cores <- ifelse(max.allow.recursive == 1L, 1L, 2L)
        pValues <- mclapply(cIndeces, samp2.testing, level + 1L, p.value, x, y, 
                            allSamp1.ids, allActSet.ids, x.nonTested, 
                            hierarchy, max.p.esti, B, gamma, 
                            max.allow.recursive, mc.cores = mc.cores, 
                            mc.allow.recursive = TRUE)
      } else {
        pValues <- lapply(cIndeces, samp2.testing,  level + 1L, p.value, x, y, 
                          allSamp1.ids, allActSet.ids, x.nonTested, hierarchy, 
                          max.p.esti, B, gamma, max.allow.recursive)
      }
    } else {
      pOne <- function(cIndex) {
        if (!is.null(cIndeces <- attr(hierarchy[[cIndex]], "subset")))
          return(c(1, sapply(cIndeces, pOne)))
        return(1)
      }
      pValues <-  sapply(cIndeces, pOne)
    }
  } else {
    pValues <- c()
  }
  out <- c(p.value, unlist(pValues))
  out
}

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


#' @title Significant clusters at alpha threshold
#' @param object a hit object.
#' @param alpha alpha level.
#' @param ... further arguments passed to or from other methods (not used).
#' @method summary hit
#' @export 
summary.hit <- function(object, alpha=.05, ...) {
  make.pVal <- function(i) {
    p.value <- object$pValues[i]
    inx <- object$hierarchy[[i]]
    if (p.value <= alpha && all(is.na(p.cluster[inx]))) {
      p.cluster[inx] <<-  p.value
      id.cluster[inx] <<- counter
      h.cluster[inx] <<- attr(object$hierarchy[[i]], "height")
      counter <<- counter + 1L
    }
    return(NULL)
  } 
  stopifnot(is(object, "hit"))
  p.cluster <- rep(NA_real_, length(object$hierarchy[[1L]]))
  id.cluster <- rep(NA_integer_, length(object$hierarchy[[1L]]))
  h.cluster <- rep(NA_real_, length(object$hierarchy[[1L]]))
  counter <- 1L
  lapply(length(object$hierarchy):1L, make.pVal)
  non.na <- which(!is.na(id.cluster))
  out <- data.frame("clusters"=id.cluster[non.na], 
                    "pValues"=p.cluster[non.na],
                    "height"=h.cluster[non.na])
  rownames(out) <- names(object$hierarchy[[1L]])[non.na]
  out
} # summary.hit

#' @title Significant hierarchy
#' @description Significant hierarchy
#' @param x a hit object
#' @details makes a matrix of p-values for image(p.matrix(x))
#' @export
p.matrix <- function (x) {
  heig <- heightSets(x$hierarchy)
  allheig <- sapply(x$hierarchy, attr, "height")
  inx <- which(heig[1] == allheig)
  p.val <- rep(x$pValues[inx], sapply(allheig[inx], length))
  out <- list(rep(NA_real_, length(p.val)))
  out[[1]][unlist(x$hierarchy[inx])] <- p.val
  for (h in 2L:length(heig)) {
    out[[h]] <- out[[h-1]]
    inx <- which(heig[h] == allheig)
    p.val <- rep(x$pValues[inx], sapply(allheig[inx] , length))
    out[[h]][unlist(x$hierarchy[inx])] <- p.val
  }
  out <- do.call("rbind", out)
  colnames(out) <- names(x$hierarchy[[1]])
  rownames(out) <- heig
  out
} # p.matrix
