mrca <- function(tips, tree) {
  
  tips <- 
    tips %>%
    unique() %>%
    {.[. %in% c(tree$tip.label, tree$node.label)]}
  
  if (length(tips) == 0) stop("tips not found in tree")
  
  if (length(tips) == 1) return(tips)
  
  tips %>%
    phytools::findMRCA(tree = tree) %>%
    {c(tree$tip.label, tree$node.label)[.]}
  
}

# descendantes include node itself! 
descendants <- function(node, tree) {
  
  node %>%
    {which(c(tree$tip.label, tree$node.label) == .)} %>%
    phytools::getDescendants(tree = tree) %>%
    {c(tree$tip.label, tree$node.label)[.]} %>%
    {c(., node)} %>%
    unique()
  
}

ancestor_if_complete <- function(tip_labels, tree) {
  
  if (length(tip_labels) == 1) return(as.character(NA))
  
  ancestor <- tip_labels %>% phytools::findMRCA(tree = tree)
  
  descendant_labels <-
    ancestor %>%
    phytools::getDescendants(tree = tree) %>%
    {tree$tip.label[.]} %>%
    purrr::discard(is.na)
  
  ancestor_label <-
    ancestor %>%
    {. - ape::Ntip(tree)} %>%
    {tree$node.label[.]}
  
  if (setequal(tip_labels, descendant_labels)) {
    return(ancestor_label)
  } else {
    return(as.character(NA))
  }
  
}

genomes_extended <- function(tg) {
  
  if(tibble::has_name(tg, "nodes")) {
    tg$genomes <- tg$genomes %>% left_join(tg$nodes, by = "node")
  }
  
  tg$genomes
  
}

complete_pairs <- function(pairs) {
  
  pairs_2 <-
    pairs %>%
    rename(genome_1_new = genome_2, genome_2_new = genome_1) %>%
    rename(genome_1 = genome_1_new, genome_2 = genome_2_new)
  
  bind_rows(pairs, pairs_2)
  
}