ancestor_if_complete <- function(tree, tip_labels) {
  
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

add_tree <- function(tog, tree) {
  
  tree <- ape::keep.tip(tree, tog$genomes$genome)
  
  node_labels <- str_c("n", 1:tree$Nnode)
  
  nodes <- tibble(support = tree$node.label, node = node_labels)
  tree$node.label <- node_labels
  
  c(tog, list(tree = tree, nodes = nodes))
  
}

map_patterns <- function(tog) {
  
  if(is.null(tog$patterns)) stop("no patterns found")
  
  nodes_patterns <-
    tog$components %>%
    group_by(pattern) %>%
    summarize(node = ancestor_if_complete(tog$tree, genome))
  
  tog %>%
    modify_at("nodes", left_join, nodes_patterns)
  
}

# assumes that phylogroups are defined in terms of species names
add_phylogroups <- function(tog, phylogroups) {
  
  genomes <- tog$genomes %>% select(genome, species)
  
  tog$phylogroups <- 
    phylogroups %>%
    left_join(
      genomes %>% rename(species_type = species, genome_type = genome)
    ) %>%
    left_join(
      genomes %>% rename(species_peripheral = species, genome_peripheral = genome)
    ) %>%
    select(phylogroup, genome_type, genome_peripheral)
  
  genomes_phylogroups <- 
    tog$phylogroups %>%
    expand_clades(genome_type, genome_peripheral, tog$tree)
  
  tog %>%
    modify_at("genomes", left_join, genomes_phylogroups)
  
}

map_phylogroups <- function(tog) {
  
  if(is.null(tog$phylogroups)) stop("no phylogroups found")
  
  phylogroups_patterns <-
    tog$genomes %>%
    group_by(phylogroup) %>%
    summarize(node = ancestor_if_complete(tog$tree, genome))
  
  tog %>%
    modify_at("phylogroups", left_join, phylogroups_patterns)
  
}

# assumes tree will be plotted in ladderized fashion
add_phylogroup_color <- function(tog, n = 12) {
  
  n_phylogroups <- nrow(tog$phylogroups)
  
  tog$tree <- ape::ladderize(tog$tree)
  
  genomes_ordered <-
    tog$tree$edge[, 2] %>%
    {.[. <= length(tog$tree$tip.label)]} %>%
    {tog$tree$tip.label[.]}
  
  tog$phylogroups <-
    tog$phylogroups %>%
    mutate(genome_type_fct = factor(genome_type, levels = !! genomes_ordered)) %>%
    arrange(genome_type_fct) %>%
    mutate(phylogroup_color = 1:UQ(n) %>% as.character() %>% rep_len(!! n_phylogroups)) %>%
    select(- genome_type_fct)
  
  tog
  
} 

ggtree_augmented <- function(tog, ...) {
  
  if (! is.null(tog$phylogroups)) {
    
    tog <- tog %>% modify_at("genomes", left_join, tog$phylogroups)
    
  }
  
  # mask variable that gives ggtree trouble
  tog$genomes$node <- NULL
  
  nodes <- 
    tog$nodes %>%
    left_join(tog$patterns) %>%
    select(label = node, everything())
  
  ggtree(tog$tree, ...) %<+%
    nodes %<+%
    tog$genomes
  
}

add_n_genomes_in_phylogroup <- function(tog) {
  
  if (! is.null(tog$phylogropus)) stop("no phylogroups found")
  
  phylogroups_n_genomes <- 
    tog$genomes %>% 
    count(phylogroup) %>%
    rename(n_genomes = n)
  
  tog %>%
    modify_at("phylogroups", left_join, phylogroups_n_genomes) 
  
}