# Returns tibble with for every genome the number of orthogroups that are absent
# in that genome but present in all other genomes of its species.
genomes_n_um_orthogroups <- function(tog) {
  
  # initiate empty character vector
  genomes_n_um_orthogroups <- tibble()
  
  # make list of species
  species_list <- unique(tog$genomes$clade) 
  
  # loop over species
  for (species in species_list) {
    
    matrix <- orthogroup_content_matrix_clade(tog, species)
    matrix <- matrix == 0
    matrix <- matrix[, apply(matrix, 2, function(x) sum(x) == 1)]
    
    genomes_n_um_orthogroups_new <- tibble(
      genome = rownames(matrix),
      n_um_orthogroups = apply(matrix, 1, sum)
    )
    
    genomes_n_um_orthogroups <- genomes_n_um_orthogroups %>%
      bind_rows(genomes_n_um_orthogroups_new)
    
  }
  
  genomes_n_um_orthogroups
  
}

# Returns tibble with for every genome the number of core orthogroups of its
# species if that genome were to be left out.
genomes_n_loo_core_orthogroups <- function(tog) {
  
  species_n_core_orthogroups <- species_n_core_orthogroups(tog)
  genomes_n_um_orthogroups <- genomes_n_um_orthogroups(tog)
  
  tog$genomes %>%
    select(genome, clade) %>%
    left_join(species_n_core_orthogroups) %>%
    left_join(genomes_n_um_orthogroups) %>%
    mutate(n_loo_core_orthogroups = n_core_orthogroups + n_um_orthogroups) %>%
    select(genome, n_loo_core_orthogroups)
  
}

# Returns a tible with for each species the mean and sd of the number of
# orthogroups per genome
# Expects "n_orthogroups_in_genome" to be present in genomes tibble. 
species_av_orthogroups <- function(tog) {
  
  tog$genomes %>%
    select(clade, genome, n_orthogroups_in_genome) %>%
    group_by(clade) %>%
    summarize(
      av_orthogroups = mean(n_orthogroups_in_genome), 
      sd_orthogroups = sd(n_orthogroups_in_genome)
    )
  
}
