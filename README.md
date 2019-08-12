# tidygenomes

## Overview

Tidygenomes is an R package for the combination, analysis and visualization of various types of genome-centered datasets.

Among other things, tidygenomes can:

* Combine various genome-centered datasets: genome tables, pangenomes, phylogenetic trees and genome pair tables (e.g. with genome-genome distances)
* Filter the data based on genome characteristics
* Add a set of phylogroups (defined as non-overlapping clades in the phylogeny)
* Compute "exclusivity" of phylogroups and related measures
* Compute orthogroup presence/absence patterns and map "signature patterns" (patterns universal in and unique to a clade) to a tree
* Compute genome-genome distances based on gene content
* Inflate species of a metapangenome (pangenome of multiple species-level pangenomes) and collapse species of a pangenome
* Create some visualizations: a heatmap of genome comparisons and an upset plot of presence/absence patterns

## Installation

Disclaimer: tidygenomes is software developed alongside some specific research projects, and the primary focus is currently on the research and not primarily on the re-usability of the package. Also, the package is still in the initial development phase.

You can install the latest, semi-stable version of tidygenomes by running the following code:

    # install.packages("devtools")
    devtools::install_github("swittouck/tidygenomes", ref = "v0.1.0")
