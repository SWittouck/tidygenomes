
add_n_genomes_in_phylogroup <- function(tog) {
  
  if (! is.null(tog$phylogropus)) stop("no phylogroups found")
  
  phylogroups_n_genomes <- 
    tog$genomes %>% 
    count(phylogroup) %>%
    rename(n_genomes = n)
  
  tog %>%
    modify_at("phylogroups", left_join, phylogroups_n_genomes) 
  
}