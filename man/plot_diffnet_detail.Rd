% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_networkComparison.R
\name{plot_diffnet_detail}
\alias{plot_diffnet_detail}
\title{Plot the dynamic differential network but colored by communities and optionally faded by betweenness}
\usage{
plot_diffnet_detail(
  grn,
  tfs,
  only_TFs = TRUE,
  order = NULL,
  weight_column = "zscore",
  compute_betweenness = TRUE
)
}
\arguments{
\item{grn}{the dynamic network}

\item{tfs}{TFs}

\item{only_TFs}{whether or not to only plot TFs and exclude non-regulators}

\item{order}{the order in which to plot epochs, or which epochs to plot}

\item{weight_column}{column name with edge weights}

\item{compute_betweenness}{whether or not to fade nodes by betweenness}
}
\value{

}
\description{
Plot the dynamic differential network but colored by communities and optionally faded by betweenness
}
