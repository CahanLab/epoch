% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_utils.R
\name{find_communities}
\alias{find_communities}
\title{Find communities in a static or dynamic network}
\usage{
find_communities(grn, use_weights = FALSE, weight_column = NULL)
}
\arguments{
\item{use_weights}{whether or not to use edge weights (for weighted graphs)}

\item{weight_column}{if using weights, name of the column containing edge weights}
}
\value{
community assignments of nodes in the dynamic network
}
\description{
Find communities in a static or dynamic network
}
