#' @title Hierarchical Inference Testing
#' 
#' @description Hierarchical inference testing for linear models with 
#' high-dimensional and/or correlated covariates by repeated sample splitting.
#' 
#' @param x Design matrix of dimension \code{n * p}, without intercept.
#' Variables not part of the dendrogram are added to the HO-model, see Details 
#' below.
#' @param y Quantitative response variable dimension \code{n}.
#' @param hierarchy Object of class \code{\link{as.hierarchy}}. Must include 
#' all variables of \code{x} which should be tested.
#' @param B Number of sample-splits.
#' @param p.samp1 Fraction of data used for the LASSO. The ANOVA uses 
#' \code{1 - p.samp1}.
#' @param sel.method Name of method for the selection of covariates via the 
#' LASSO, either "ASF" (default), the occurrence frequenze of covariates along 
#' the lambda path or "CV" n-fold cross-validation.
#' @param activeset.freq Frequency in which a covariate has to occur along the 
#' lambda path to be selected in the active set. Ignored if \code{sel.method} 
#' is "CV".
#' @param nfolds Number of folds (default is 10). Ignored if \code{sel.method} 
#' is "ASF". See \code{\link[glmnet]{cv.glmnet}} for more details. 
#' @param lambda.opt Criterion for optimum selection of cross validated lasso. 
#' Either "lambda.1se" (default) or "lambda.min". Ignored if \code{sel.method} 
#' is "ASF". See \code{\link[glmnet]{cv.glmnet}} for more details. 
#' @param gamma Vector of gamma-values.
#' @param max.p.esti Maximum alpha level. All p-values above this value are set 
#' to one. Small \code{max.p.esti} values reduce computing time.
#' @param mc.cores Number of cores for parallelising. Theoretical maximum is 
#' 'B'. For details see \code{\link[parallel]{mclapply}}.
#' @param trace If TRUE it prints current status of the program.
#' @param ... Additional arguments for \code{\link[glmnet]{cv.glmnet}}.
#' 
#' @details The H0-model contains variables, with are not tested, like 
#' experimental-design variables. These variables are not penalised in the 
#' LASSO model selection and are always include in the reduced ANOVA model.
#' 
#' @references Mandozzi, J. and Buehlmann, P. (2013). \emph{Hierarchical 
#' testing in the high-dimensional setting with correlated variables}. To 
#' appear in the Journal of the American Statistical Association. Preprint 
#' arXiv:1312.5556
#'   
#' @examples
#' set.seed(123)
#' n <- 80
#' p <- 90
#' # x with correlated columns
#' corMat <- toeplitz((p:1/p)^3)
#' corMatQ <- chol(corMat)
#' x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
#' colnames(x) <- paste0("x", 1:p)
#' # y
#' y <- x[, c(3, 5, 73)] %*% c(2, 5, 3) + rnorm(n)
#' # clustering
#' hc <- hclust(dist(t(x)))
#' # HIT with ASF
#' out <- hit(x, y, hc)
#' summary(out)
#' 
#' 
#' # HIT with CV
#' \dontrun{
#' out2 <- hit(x, y, hc, sel.method = "CV")
#' summary(out2)
#' }
#' 
#' @importFrom parallel mclapply
#' @importFrom stats reorder
#' @export 
hit <- function(x, y, hierarchy, B = 50, p.samp1 = 0.5, 
                sel.method = c("ASF", "CV"), activeset.freq = 1 / 3,
                nfolds = 10, lambda.opt = c("lambda.1se", "lambda.min"), 
                gamma = seq(0.05, 0.99, length.out = 100), max.p.esti = 1, 
                mc.cores = 1L, trace = FALSE, ...) {
  #   Mandozzi and Buehlmann (2015), 2 Description of method
  if (is.null(colnames(x)))
    stop("column names of 'x' are missing")
  n <- nrow(x)
  p <- ncol(x)
  if (class(hierarchy) == "dendrogram" || class(hierarchy) == "hclust")
    hierarchy <- as.hierarchy(hierarchy)
  stopifnot(class(hierarchy) == "hierarchy")
  ##### Checks order of variables
  x.names <- colnames(x)
  hier.names <- names(hierarchy)
  if (length(setdiff(hier.names, x.names)))
    stop("'hierarchy' includs variabels not in 'x'")
  if (identical(hier.names, x.names)) {
    additionalCovariates <- character(0L)
    x.nonTested <- integer(0L)
  } else if (setequal(hier.names, x.names)) {
    hierarchy <- reorder(hierarchy, x.names)
    additionalCovariates <- character(0L)
    x.nonTested <- integer(0L)
  } else {
    hierarchy <- reorder(hierarchy, x.names)
    additionalCovariates <- setdiff(x.names, names(hierarchy))
    x.nonTested <- match(additionalCovariates, x.names)
  }
  penalty.factor <- rep(1L, p)
  penalty.factor[x.nonTested] <- 0L
  ##### Subsample splits
  if (p.samp1 < .1 | p.samp1 > .9)
    stop("'p.samp1' must be between .1 and 0.9")
  n.samp1 <- as.integer(n * p.samp1)
  n.samp2 <- n - n.samp1
  allSamp1.ids <- replicate(B, sample.int(n, n.samp1), simplify = FALSE)
  ##  2.2 Screening
  if (trace)
    cat("LASSO has started at:\n\t", as.character(Sys.time()), "\n")
  allActSet.ids <- mclapply(allSamp1.ids, samp1.lasso,
                            x, y, n.samp2, penalty.factor, sel.method, 
                            activeset.freq, nfolds, lambda.opt, ...,
                            mc.cores = mc.cores, mc.cleanup = TRUE)
  ##  2.3 Testing and multiplicity adjustmen; and 
  ##  2.4 Aggregating and Hierarchical adjustment
  if (trace)
    cat("Significance testing has started at:\n\t", 
        as.character(Sys.time()), "\n")
  pValues <- samp2.sigHierarchy(1L, 0L, 0, x, y, allSamp1.ids, allActSet.ids, 
                                x.nonTested, hierarchy, 
                                max.p.esti, B, sort(gamma), 
                                mc.cores)
  ##### Results
  if (trace)
    cat("HIT has finished at:\n\t", as.character(Sys.time()), "\n")
  asi <- sort(unlist(allActSet.ids))
  sel.tab <- rep(0, p)
  sel.tab[unique(asi)] <- table(asi) / B
  out <- list(pValues = pValues,
              selectFreq = sel.tab,
              hierarchy = hierarchy,
              additionalCovariates = additionalCovariates)
  class(out) <- "hit"
  out
}


#' @title Variabel Screening
#' 
#' @description LASSO function of the HIT algorithem.
#' 
#' @param samp1 List of index for subsample (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param n.samp2 Number of individuals in samp2 which is the max. 
#' for non zero coefficients.
#' @param penalty.factor See glmnet.
#' @param sel.method Name of method for the selection of covariates via the 
#' LASSO, either "ASF" (default), the occurrence frequenze of covariates along 
#' the lambda path or "CV" n-fold cross-validation.
#' @param activeset.freq Frequency in which a covariate has to occur along the 
#' lambda path to be selected in the active set. Ignored if \code{sel.method} 
#' is "CV".
#' @param nfolds Number of folds (default is 10). Ignored if \code{sel.method} 
#' is "ASF". See \code{\link[glmnet]{cv.glmnet}} for more details. 
#' @param lambda.opt Criterion for optimum selection of cross validated lasso. 
#' Either "lambda.1se" (default) or "lambda.min". Ignored if \code{sel.method} 
#' is "ASF". See \code{\link[glmnet]{cv.glmnet}} for more details. 
#' @param ... Additional agruments
#' 
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom Matrix which
#' @importFrom stats coef
#' @keywords internal
samp1.lasso <- function(samp1, x, y, n.samp2, penalty.factor, sel.method, 
                        activeset.freq, nfolds, lambda.opt, ...) {
  ##  2.2 Screening
  if (match.arg(sel.method, c("ASF", "CV")) == "CV") {
    # n-fold cross-validation
    lambda.opt <- match.arg(lambda.opt, c("lambda.1se", "lambda.min"))
    lasso.fit <- cv.glmnet(x[samp1, ], y[samp1], 
                           penalty.factor = penalty.factor, 
                           nfolds = nfolds, dfmax = n.samp2 - 2L, ...)
    if (lambda.opt == "lambda.min")
      beta <- coef(lasso.fit, s = lasso.fit$lambda.min)[-1L]
    else 
      beta <- coef(lasso.fit, s = lasso.fit$lambda.1se)[-1L]
    actSet <- which(beta != 0 & penalty.factor == 1L)
  } else {
    # Select by the occurrence frequency along the lambda path
    lasso.fit <- glmnet(x[samp1, ], y[samp1], 
                        penalty.factor = penalty.factor, ...)
    inx <- which(lasso.fit$beta != 0, arr.ind = TRUE)[, 1L]
    inxcount <- table(inx)
    uniqinx <- sort(unique(inx))
    asThre <- ceiling(ncol(lasso.fit$beta) * activeset.freq)
    actSet <- uniqinx[inxcount > asThre]
    while (length(actSet) > n.samp2 - 2L) {
      asThre <- asThre + 1L
      actSet <- uniqinx[inxcount > asThre]
    }
  }
  actSet
}


#' @title Variabel Testing
#' 
#' @description ANOVA Testing, Multiplicity Adjustment, Aggregating and 
#' Hierarchical Adjustment
#' 
#' @param cIndex Index for cluster (mclapply index).
#' @param level Hierarchy level counter for parallelism
#' @param upper.p P value upper huerarchy level of the clusters variables.
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param allSamp1.ids  List of subsampels.
#' @param allActSet.ids List of active sets.
#' @param x.nonTested  Vector of indeces of non tested variabels.
#' @param hierarchy A hierarchy object
#' @param max.p.esti Maximum alpha level. All p-values above this value are set 
#' to one. Small max.p.esti values reduce computing time.
#' @param B Number of sample-splits.
#' @param gamma Vector of gamma-values.
#' @param cores Number of cores for parallelising.
#' 
#' @importFrom stats quantile
#' @keywords internal
samp2.sigHierarchy <- function(cIndex, level, upper.p, x, y, allSamp1.ids, 
                               allActSet.ids, x.nonTested, hierarchy, 
                               max.p.esti, B, gamma, cores) {
  ## 2.3 Testing and multiplicity adjustment
  cluster <- hierarchy[[cIndex]]
  p.cluster <- sapply(1L:B, samp2.sigNode, 
                      x, y, cluster, x.nonTested, 
                      allSamp1.ids, allActSet.ids)
  ##  2.4 Aggregating and Hierarchical adjustment
  ### 2.4-1 Aggregating 
  q.aggre <- sapply(gamma, 
                    function(i, x) { min(1, quantile(x / i, i)) }, 
                    x = p.cluster)
  p.aggre <- min(1, (1 - log(min(gamma))) * min(q.aggre))
  ### 2.4-2 Hierarchical adjustment
  p.value <- max(p.aggre, upper.p)
  ##### Estimation at next lower level and find a way to parallelize
  if (!is.null(cIndeces <- attr(cluster, "subset"))) {
    if (p.value < max.p.esti) {
      if (level == 0L && length(cIndeces) >= cores) {
        level <- as.integer(sqrt(cores))
        pValues <- mclapply(cIndeces, samp2.sigHierarchy, level + 1L, p.value, 
                            x, y, allSamp1.ids, allActSet.ids, x.nonTested, 
                            hierarchy, max.p.esti, B, gamma, cores, 
                            mc.cores = cores)
      } else if (level <= as.integer(sqrt(cores))) {
        mc.cores <- ifelse(as.integer(sqrt(cores)) == 1L, 1L, 2L)
        pValues <- mclapply(cIndeces, samp2.sigHierarchy, level + 1L, p.value, 
                            x, y, allSamp1.ids, allActSet.ids, x.nonTested, 
                            hierarchy, max.p.esti, B, gamma, cores, 
                            mc.cores = mc.cores, mc.allow.recursive = TRUE)
      } else {
        pValues <- lapply(cIndeces, samp2.sigHierarchy,  level + 1L, p.value, 
                          x, y, allSamp1.ids, allActSet.ids, x.nonTested, 
                          hierarchy, max.p.esti, B, gamma, cores)
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
#' 
#' @description ANOVA test (at single node) of the HIT algorithem.
#' 
#' @param k Index for subsample (mclapply index).
#' @param x Design matrix, of dimension n x p.
#' @param y Vector of quantitative response variable.
#' @param cluster Clusters to be tested.
#' @param x.nonTested Vector of indeces of non tested variabels.
#' @param allSamp1.ids  List of subsampels.
#' @param allActSet.ids List of active sets.
#' 
#' @keywords internal
samp2.sigNode <- function(k, x, y, cluster, x.nonTested, 
                          allSamp1.ids, allActSet.ids) {
  ## 2.3 Testing and multiplicity adjustment
  n <- nrow(x)
  actClust <- intersect(allActSet.ids[[k]], cluster)
  nonActClust <- setdiff(allActSet.ids[[k]], cluster)
  ### 2.3-1 Testing
  p.cluster <- 1
  if (l.actClust <- length(actClust)) {
    ##### ANOVA between active set and active set minus cluster
    samp2 <- (1L:n)[-allSamp1.ids[[k]]]
    y <- y[samp2]
    if (l.nonActClust <- length(nonActClust)) {
      if (l.nonTested <- length(x.nonTested)) {
        x  <- cbind(1L, x[samp2, c(x.nonTested, nonActClust, actClust)])
        assign <- rep(0L:3L, c(1L, l.nonTested, l.nonActClust, l.actClust))
        get.p <- 3L
      } else {
        x  <- cbind(1L, x[samp2, c(nonActClust, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonActClust, l.actClust))
        get.p <- 2L
      }
    } else {
      if (l.nonTested <- length(x.nonTested)) {
        x  <- cbind(1L, x[samp2, c(x.nonTested, actClust)])
        assign <- rep(0L:2L, c(1L, l.nonTested, l.actClust))
        get.p <- 2L
      } else {
        x  <- cbind(1L, x[samp2, actClust])
        assign <- rep(0L:1L, c(1L, l.actClust))
        get.p <- 1L
      }
    }
    p.test <- fast.anova(x, y, assign)[get.p]
    ### 2.3-2 Multiplicity adjustment
    p.cluster <- min(1, p.test * (l.actClust + l.nonActClust) / l.actClust)
  }
  p.cluster
}


#' @title Summary of HIT
#' 
#' @description Significant clusters at alpha threshold.
#' 
#' @param object A \code{\link{hit}} object.
#' @param alpha A alpha significance threshold.
#' @param max.height max. Height to consider.
#' @param ... Further arguments passed to or from other methods (not used).
#' 
#' @method summary hit
#' @export 
summary.hit <- function(object, alpha = 0.05, max.height, ...) {
  make.pVal <- function(i) {
    p.value <- object$pValues[i]
    inx <- object$hierarchy[[i]]
    if (p.value <= alpha && 
          (all(is.na(P.CLUSTER[inx])) || 
             (all(!is.na(P.CLUSTER[inx]) && 
                    P.CLUSTER[inx] <= alpha)))) {
      P.CLUSTER[inx] <<-  p.value
      ID.CLUSTER[inx] <<- COUNTER
      H.CLUSTER[inx] <<- attr(object$hierarchy[[i]], "height")
      COUNTER <<- COUNTER + 1L
    }
    return(NULL)
  } 
  stopifnot(inherits(object, "hit"))
  P.CLUSTER <- rep(NA_real_, length(object$hierarchy[[1L]]))
  ID.CLUSTER <- rep(NA_integer_, length(object$hierarchy[[1L]]))
  H.CLUSTER <- rep(NA_real_, length(object$hierarchy[[1L]]))
  COUNTER <- 1L
  lapply(length(object$hierarchy):1L, make.pVal)
  non.na <- which(!is.na(ID.CLUSTER))
  out <- data.frame(clusters = ID.CLUSTER[non.na], 
                    heights = H.CLUSTER[non.na],
                    pValues = P.CLUSTER[non.na])
  rownames(out) <- names(object$hierarchy[[1L]])[non.na]
  if (!missing(max.height)) 
    out <- out[out[, 2L] <= max.height, ]
  if (ll <- length(unique(out[, 1L])))
    out[, 1L] <- as.integer(factor(out[, 1L], labels = 1L:ll))
  out
}
