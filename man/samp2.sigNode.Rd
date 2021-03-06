% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hit.R
\name{samp2.sigNode}
\alias{samp2.sigNode}
\title{ANOVA testing and multiplicity adjustment}
\usage{
samp2.sigNode(k, x, y, cluster, family, test, x.notest, allSamp1.ids,
  allActSet.ids)
}
\arguments{
\item{k}{Index for subsample (mclapply index).}

\item{x}{Design matrix, of dimension n x p.}

\item{y}{Vector of quantitative response variable.}

\item{cluster}{Clusters to be tested.}

\item{family}{Distribution family of \code{y}.}

\item{test}{name of test.}

\item{x.notest}{Vector of indeces of non tested variabels.}

\item{allSamp1.ids}{List of subsampels.}

\item{allActSet.ids}{List of active sets.}
}
\description{
ANOVA test (at single node) of the HIT algorithem.
}
\keyword{internal}

