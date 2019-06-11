#' Filter the genomes of a tidygenomes object
#'
#' This function applies [dplyr::filter()] to the genomes table of a tidygenomes
#' object and subsequently filters all other components of the tidygenomes
#' object to retain only the requested genomes. It is also possible to filter on
#' variables of the phylogroups table, if present.
#'
#' @param tg A tidygenomes object
#' @param ... Filtering expression to pass on to [dplyr::filter()]
#' 
#' @return A tidygenomes object
#'
#' @examples
#' genomes <- data.frame(genome = c("genome 1", "genome 2"))
#' tg <- as_tidygenomes(genomes)
#' tg <- filter_genomes(genomes, genome == "genome 1")
#' 
#' @export
filter_genomes <- function(tg, ...) {
  
  if (! is.null(tg$patterns)) {
    stop("Tidygenomes object contains patterns")
  }
  
  genomes <-
    genomes_extended(tg) %>%
    filter(...) %>%
    pull(genome)
  
  tg$genomes <- filter(tg$genomes, genome %in% !! genomes)
  
  if (! is.null(tg$genes)) {
    tg$genes <- filter(tg$genes, genome %in% tg$genomes$genome)
  }
  
  if (! is.null(tg$orthogroups)) {
    tg$orthogroups <- 
      filter(tg$orthogroups, orthogroup %in% tg$genes$orthogroup)
  }
  
  if (! is.null(tg$pairs)) {
    tg$pairs <- 
      filter(
        tg$pairs, 
        genome_1 %in% tg$genomes$genome, genome_2 %in% tg$genomes$genome
      )
  }
  
  if (! is.null(tg$nodes)) {
    tips_to_keep <- 
      tg$genomes %>% 
      left_join(tg$nodes, by = "node") %>%
      pull(node)
    tg$tree <- ape::keep.tip(tg$tree, tips_to_keep)
    tg$nodes <- 
      filter(tg$nodes, node %in% c(tg$tree$tip.label, tg$tree$node.label))
  }
  
  if (! is.null(tg$species)) {
    
    tg$species <- filter(tg$species, species %in% tg$genomes$species)
    
  }
  
  tg
  
}

#' Filter the orthogroups of a tidygenomes object
#'
#' This function applies [dplyr::filter()] to the orthogroups table of a
#' tidygenomes object and subsequently filters other components of the
#' tidygenomes object to retain only the requested orthogroups.
#'
#' @param tg A tidygenomes object
#' @param ... Filtering expression to pass on to [dplyr::filter()]
#' 
#' @return A tidygenomes object
#' 
#' @export
filter_orthogroups <- function(tg, ...) {
  
  tg %>%
    modify_at("orthogroups", filter, ...) %>%
    modify_at("genes", filter, orthogroup %in% .$orthogroups$orthogroup)
  
}