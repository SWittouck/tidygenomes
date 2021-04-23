#' Return the tree of a tidygenomes object
#'
#' Internally, the tree of a tidygenomes object is stored with n-numbers as tip
#' labels. This function returns a tree with genome names as tips.
#'
#' @param tg A tidygenomes object
#'
#' @return A phylo object
#'
#' @export
get_tree <- function(tg) {
  
  # replace tipnode names with genome names
  lut_tips <- structure(tg$genomes$genome, names = tg$genomes$node)
  
  tree <- tg$tree 
  tree$tip.label <- lut_tips[tree$tip.label]
  names(tree$tip.label) <- NULL
  
  # return tree
  tree
  
}

#' Return a pangenome matrix
#'
#' This function returns the orthogroup content of the genomes in the form of a
#' matrix where the rows represent genomes, the columns represent orthogroups
#' and the cells contain the copy number or presence status of the orthogroups
#' in the genomes.
#'
#' @param tg A tidygenomes object
#' 
#' @return A matrix with the pangenome
#' 
#' @export
pangenome_matrix <- function(tg) {
  
  tg$genes %>%
    count(genome, orthogroup) %>%
    spread(key = orthogroup, value = n, fill = 0) %>%
    `class<-`("data.frame") %>%
    `rownames<-`(.$genome) %>%
    select(- genome) %>%
    as.matrix()
  
}

#' Return a genome table
#'
#' This function returns the genome table of a tidygenomes object, possibly
#' extended with extra genome-related data such as phylogroups.
#'
#' @param tg A tidygenomes object
#' @param extend Whether to add genome-related metadata
#' 
#' @return A genome table
#' 
#' @export
genomes <- function(tg, extend = F) {
  
  if(extend & tibble::has_name(tg, "nodes")) {
    tg$genomes <- tg$genomes %>% left_join(tg$nodes, by = "node")
  }
  
  if(extend & tibble::has_name(tg, "phylogroups")) {
    tg$genomes <- tg$genomes %>% left_join(tg$phylogroups, by = "phylogroup")
  }
  
  tg$genomes
  
}

#' Return a pair table
#'
#' This function returns the pair table of a tidygenomes object, possibly
#' completed with inverse pairs and possibly extended with extra
#' genome-related data such as phylogroups.
#'
#' @param tg A tidygenomes object
#' @param extend Whether to add genome-related metadata
#' @param complete Whether to add inverse pairs
#' @param add_reflections Whether to add reflexive pairs
#' 
#' @return A pair table
#' 
#' @export
pairs <- function(tg, extend = F, complete = F, add_reflections = F) {
  
  if (complete) {
    pairs <- complete_pairs(tg$pairs, genome_1, genome_2)
  } else {
    pairs <- tg$pairs
  }
  
  if (add_reflections) {
    reflections <- 
      tibble(genome_1 = tg$genomes$genome, genome_2 = tg$genomes$genome)
    pairs <- pairs %>% bind_rows(reflections)
  }
  
  if (extend) {
    
    genomes <- genomes(tg, extend = T) 
    
    pairs <-
      pairs %>%
      left_join(genomes %>% rename(genome_1 = genome), by = "genome_1") %>%
      left_join(
        genomes %>% rename(genome_2 = genome), by = "genome_2", 
        suffix = c("_1", "_2")
      )
    
  }
  
  pairs
  
}

orthogroups <- function(tg) {
  
  tg$orthogroups
  
}
