expand_clade <- function(tree, tip_labels) {
  
  tip_labels <- 
    tip_labels %>%
    unique() %>%
    purrr::keep(~ . %in% tree$tip.label) 
  
  if (length(tip_labels) <= 1) return(tip_labels)
  
  tip_labels %>%
    phytools::findMRCA(tree = tree) %>%
    phytools::getDescendants(tree = tree) %>%
    {tree$tip.label[.]} %>%
    purrr::discard(is.na)
  
}

expand_clades <- function(clades, genome_1 = genome_1, genome_2 = genome_2, tree) {
  
  genome_1 <- rlang::enexpr(genome_1)
  genome_2 <- rlang::enexpr(genome_2)
  
  clades %>%
    mutate(genomes = map2(
      !! genome_1, !! genome_2, 
      ~ expand_clade(!! tree, c(.x, .y))
    )) %>%
    unnest(genomes) %>%
    rename(genome = genomes) %>%
    select(- !! genome_1, - !! genome_2)
  
}

# function to convert tibble to matrix
# (use as_tibble to do the opposite)
as_matrix <- function(tibble, rownames = NULL) {
  
  if (! is.null(rownames)) {
    
    tibble %>%
      `class<-`("data.frame") %>%
      `rownames<-`(.[[rownames]]) %>%
      select(- one_of(rownames)) %>%
      as.matrix()
    
  } else {
    
    as.matrix(tibble)
    
  }
  
}