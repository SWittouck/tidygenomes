
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

add_n_genomes_in_phylogroup <- function(tog) {
  
  if (! is.null(tog$phylogropus)) stop("no phylogroups found")
  
  phylogroups_n_genomes <- 
    tog$genomes %>% 
    count(phylogroup) %>%
    rename(n_genomes = n)
  
  tog %>%
    modify_at("phylogroups", left_join, phylogroups_n_genomes) 
  
}