#' Add a branch to the root of the tree
#' 
#' This function adds a branch to the root of the tree, adding one extra node in
#' the process.
#' 
#' @param tree A phylogeny of class `phylo`, with edge lengths
#' @param branch_length The length of the root branch
#' 
#' @return An object of class `phylo`
#' 
#' @example 
#' tree <- ape::read.tree(text = "(a:1, b:1);")
#' tree <- add_rootbranch(tree)
#' plot(tree)
#' 
#' @export
add_rootbranch <- function(tree, branch_length = 1) {
  
  n_tips <- length(tree$tip.label)
  tree$edge[tree$edge > n_tips] <- tree$edge[tree$edge > n_tips] + 1
  tree$edge <- rbind(c(n_tips + 1, n_tips + 2), tree$edge)
  tree$edge.length <- c(branch_length, tree$edge.length)
  tree$Nnode <- tree$Nnode + 1
  
  tree
  
}

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
  
  if(tibble::has_name(tg, "phylogroups")) {
    tg$genomes <- tg$genomes %>% left_join(tg$phylogroups, by = "phylogroup")
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

tipnodes_ladderized <- function(tree) {
  
  tree <- ape::ladderize(tree)
  
  tree$edge[, 2] %>%
  {.[. <= length(tree$tip.label)]} %>%
  {tree$tip.label[.]}
  
}

# function to convert various data types to tibble
as_tibble <- function(data) {
  
  if("dist" %in% class(data)) {
    
    n <- attr(data, "Size")
    
    expand.grid(i = 1:n, j = 1:n) %>%
      filter(i < j) %>%
      mutate(object_1 = labels(data)[i]) %>%
      mutate(object_2 = labels(data)[j]) %>%
      mutate(distance = data[n * (i - 1) - i * (i - 1) / 2 + j - i]) %>%
      select(- i, - j)
    
  }
  
}