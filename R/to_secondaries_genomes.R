# Returns tog object, with added to the genomes table the first and second
# components of a PCoA on orthogroup content.
add_pcoa <- function(tog) {
  
  gene_count_matrix <- get_gene_count_matrix(tog)
  dist_matrix <- vegan::vegdist(gene_count_matrix, method = "jaccard", binary = T)
  
  # perform PCoA
  pcoa <- cmdscale(dist_matrix, k = 2, eig = T, list = T)
  pcoa_variances <- pcoa$eig/sum(pcoa$eig)
  pcoa_dimensions <- pcoa$points %>%
    as_tibble() %>%
    mutate(genome = !! rownames(pcoa$points)) %>%
    rename(pcoa_1 = V1, pcoa_2 = V2)
  
  # add PCoA dimensions to sample table
  tog$genomes <- tog$genomes %>%
    left_join(pcoa_dimensions)
  
  # add PCoA variances to ta object
  tog$pcoa_variances <- pcoa_variances
  
  tog
  
}