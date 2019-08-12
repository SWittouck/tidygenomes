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
#' tidygenomes object to retain only the requested orthogroups. This function
#' will not remove genomes that have zero genes left after orthogroup filtering.
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

#' Update the names of the genomes
#'
#' This function updates the names of all genomes in all tables.
#' 
#' This function is deliberately not called "rename_genomes", because in the
#' tidyverse naming framework that would suggest that the function renames
#' variables in the genome table (which it doesn't do).
#'
#' @param tg A tidygenomes object
#' @param new_name An expression that evaluates to unique names within the
#'   genome table
#' 
#' @return A tidygenomes object
#' 
#' @export
update_genomes <- function(tg, new_name) {
  
  new_name <- rlang::enexpr(new_name)
  
  lut_genomes <-
    tg$genomes %>%
    mutate(new_name = !! new_name) %>%
    {structure(.$new_name, names = .$genome)}
  
  is_unique <- function(x) length(x) == length(unique(x))
  
  if (! is_unique(lut_genomes)) {
    stop("the new genome names are not unique")
  }
  
  tg$genomes <- tg$genomes %>% mutate(genome = !! new_name)
  
  if (! is.null(tg$genes)) {
    tg$genes <- tg$genes %>% mutate(genome = lut_genomes[genome] %>% unname())
  }
  
  if (! is.null(tg$pairs)) {
    
    tg$pairs <-
      tg$pairs %>%
      mutate(
        genome_1 = lut_genomes[genome_1] %>% unname(),
        genome_2 = lut_genomes[genome_2] %>% unname()
      )
    
    for (row in 1:nrow(tg$pairs)) {
      
      genome_1 <- tg$pairs[row, "genome_1"]
      genome_2 <- tg$pairs[row, "genome_2"]
      
      if (genome_1 < genome_2) {
        tg$pairs[row, "genome_1"] <- genome_2
        tg$pairs[row, "genome_2"] <- genome_1
      }
      
    }
    
  }
  
  if (! is.null(tg$nodes)) {
    
  }
  
  tg
  
}