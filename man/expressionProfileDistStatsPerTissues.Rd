% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_funks.R
\name{expressionProfileDistStatsPerTissues}
\alias{expressionProfileDistStatsPerTissues}
\title{Extract different \code{stats} per tissue specific expression profile distances.}
\usage{
expressionProfileDistStatsPerTissues(expr.prof.dists.per.tissue, gene.classes,
  fam.name = NA, stats = c("max", "median", "mean", "maxMinusMin"))
}
\arguments{
\item{expr.prof.dists.per.tissue}{An instance of list with names identifying
tissues and values being instances of \code{base::dist}. See function
\code{expressionProfilesDists(..., per.tissue = TRUE)}.}

\item{gene.classes}{An instance of \code{list} where the names are the
different gene classes and the values character vectors intersecting with
the argument genes and thus segregating them into the respective
gene-classes.}

\item{fam.name}{The name of the group of genes the distances were
computed for. Will be the row-name of the returned numeric matrix. Default
is \code{NA}.}

\item{stats}{A character vector holding the names of the functions to apply
on the the tissue specific distances. Default is \code{stats = c('max',
'median', 'mean', 'maxMinusMin')}.}
}
\value{
An instance of \code{matrix} with numeric values holding the stats
inferred for each group of tissue specific distances as found in argument
\code{expr.prof.dists.per.tissue}. Returns \code{NULL} if
\code{all(is.na(expr.prof.dists.per.tissue))}.
}
\description{
Extract different \code{stats} per tissue specific expression profile distances.
}

