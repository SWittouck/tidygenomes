#' Call `ggtree()` and add metadata
#'
#' This function calls [ggtree::ggtree()] on the tree component of a tidygenomes object and
#' adds all available genome and node metadata.
#'
#' @param tg A tidygenomes object
#' @param ... Extra arguments to pass on to `ggtree()`
#' 
#' @return A ggplot object
#' 
#' @importFrom ggtree %<+%
#' @export
ggtree_augmented <- function(tg, ...) {
  
  if (! is.null(tg$phylogroups)) {
    
    tg <- tg %>% modify_at("nodes", left_join, tg$phylogroups, by = "phylogroup")
    
  }
  
  if (! is.null(tg$patterns)) {
    
    tg <- tg %>% modify_at("nodes", left_join, tg$patterns, by = "node")
    
  }
  
  nodes <- 
    tg$nodes %>%
    left_join(tg$genomes, by = "node") %>%
    select(label = node, everything())
  
  ggtree::ggtree(tg$tree, ...) %<+%
    nodes 
  
}