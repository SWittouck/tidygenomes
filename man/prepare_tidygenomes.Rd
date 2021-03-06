% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importers.R
\name{prepare_tidygenomes}
\alias{prepare_tidygenomes}
\title{Prepare tidygenomes object}
\usage{
prepare_tidygenomes(
  genomes,
  pangenome = NULL,
  tree = NULL,
  phylogroups = NULL,
  root = NULL,
  genome_identifier = genome
)
}
\arguments{
\item{genomes}{A genome table with a column \code{genome}}

\item{tree}{An object of class phylo with tips corresponding to genomes}

\item{root}{Three tips that identify the root (see \link{root_tree.phylo})}
}
\value{
A tidygenomes object
}
\description{
This function prepares a tidygenomes object from a genome table, a
phylogenetic tree and a root location in the tree.
}
