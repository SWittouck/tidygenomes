# These functions are for computing information from a tog object and returing
# it. The function names follow the logic of just specifying what they return,
# as seems to be the custom in R (e.g. function class() or length()).

# Naming of tibbles that might be returned: 
# <what_do_the_rows_represent>_<information_in_extra_variables>

# synonym to avoid devtools::load_all() freaking out
metaorthogroup_content_matrix <- pangenome_matrix

# Returns the orthogroup content of a all genomes of a species in the form of a
# matrix.
orthogroup_content_matrix_species <- function(tog, species) {
  
  tog$genes %>%
    left_join(tog$genomes) %>%
    filter(clade == !! species) %>%
    left_join(tog$orthogroups) %>%
    count(genome, orthogroup_name) %>%
    spread(key = orthogroup_name, value = n, fill = 0) %>%
    `class<-`("data.frame") %>%
    `rownames<-`(.$genome) %>%
    select(- genome) %>%
    as.matrix()
  
}

# Returns tibble with core orthogroup count per species
species_n_core_orthogroups <- function(tog, exceptions = 0) {
  
  # tibble with number of genomes for each species
  species_n_genomes <- tog$genomes %>%
    count(clade) %>%
    rename(n_genomes_in_species = n)
  
  tog$genes %>%
    left_join(tog$genomes) %>%
    count(clade, orthogroup, genome) %>%
    rename(n_copies = n) %>%
    left_join(species_n_genomes) %>%
    filter(n_copies != 0) %>%
    group_by(clade) %>%
    count(orthogroup, clade, n_genomes_in_species) %>%
    filter(n >= n_genomes_in_species - !! exceptions) %>%
    summarize(n_core_orthogroups = n())
  
}

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

# Returns a character vector with orthogroup names in order of hierarchical
# clustering on presence in genomes.
orthogroups_ordered <- function(tog) {
  
  # initiate empty character vector
  orthogroups_ordered <- character()
  
  # make list of species
  species_list <- unique(tog$genomes$clade) 
  
  # loop over species
  for (species in species_list) {
    
    # hierarchical clustering of orthogroups on presence in genomes
    clust_out <- orthogroup_content_matrix_clade(tog, species = species) %>%
      t() %>%
      vegan::vegdist(method = "bray", binary = T) %>%
      hclust(method = "average") 
    
    # extract ordered orthogroups from clustering output
    orthogroups_ordered_new <- clust_out$labels[clust_out$order]
    
    orthogroups_ordered <- c(orthogroups_ordered, orthogroups_ordered_new)
    
  }
  
  orthogroups_ordered
  
}

# Returns a character vector with genome names in order of hierarchical
# clustering on orthogroup or metaorthogroup content. 
genomes_ordered <- function(tog, method = "metaorthogroups") {
  
  if (method == "metaorthogroups") return(genomes_ordered_metaorthogroups(tog))
  if (method == "orthogroups") return(genomes_ordered_orthogroups(tog))
  
  stop("supply a valid method: either 'orthogroups' or 'metaorthogroups'")
  
}

# See function genomes_ordered()
genomes_ordered_orthogroups <- function(tog) {
  
  # initiate empty character vector
  genomes_ordered <- character()
  
  # make list of species
  species_list <- unique(tog$genomes$clade)
  
  # loop over species
  for (species in species_list) {
    
    # hierarchical clustering of genomes on orthogroup content
    clust_out <- orthogroup_content_matrix_clade(tog, species = species) %>%
      vegan::vegdist(method = "bray", binary = T) %>%
      hclust(method = "average") 
    
    genomes_ordered_new <- clust_out$labels[clust_out$order]
    
    genomes_ordered <- c(genomes_ordered, genomes_ordered_new)
    
  }
  
  genomes_ordered
  
}

# See function genomes_ordered()
genomes_ordered_metaorthogroups <- function(tog) {
  
  # hierarchical clustering of genomes on metaorthogroup content
  clust_out <- metaorthogroup_content_matrix(tog) %>%
    vegan::vegdist(method = "bray", binary = T) %>%
    hclust(method = "average") 
  
  genomes_ordered <- clust_out$labels[clust_out$order]
  
  genomes_ordered
  
}
