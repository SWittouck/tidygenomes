% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dist2pairs}
\alias{dist2pairs}
\title{Convert dist object to pair table}
\usage{
dist2pairs(dist)
}
\arguments{
\item{dist}{An object of the class "dist"}
}
\value{
A tibble
}
\description{
This functions converts an object of the class "dist" to a tidy table with
the columns object_1, object_2 and distance.
}
