expand_clade <- function(tree, tip_labels) {
  
  tip_labels <- 
    tip_labels %>%
    unique() %>%
    keep(~ . %in% tree$tip.label) 
  
  if (length(tip_labels) <= 1) return(tip_labels)
  
  mrca <- phytools::findMRCA(tree, tip = tree$tip.label %in% tip_labels %>% which()) 
  descendants <- phytools::getDescendants(tree, node = mrca) 
  tree$tip.label[descendants] %>% 
    na.omit() %>%
    `attr<-`("na.action", NULL)
  
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