#' Add presence/absence patterns to a tidygenomes object
#'
#' This function determines unique orthogroup presence/absence patterns in
#' genomes and adds this information to a tidygenomes object.
#' 
#' The following will be added to the tidygenomes object:
#' 
#' * A variable `pattern` in the orthogroups table, with for each orthogroup its
#' presence/absence pattern.
#' * A table `components`, listing the genomes defining each presence/absence
#' pattern.
#' * A table `patterns` with for each pattern, the number of orthogroups that
#' follow it.
#'
#' @param tg A tidygenomes object
#' 
#' @return A tidygenomes object
#'
#' @examples
#' genes <- data.frame(
#'   gene = c("gene A", "gene B", "gene C"),
#'   genome = c("genome 1", "genome 1", "genome 2"),
#'   orthogroup = c("og alpha", "og beta", "og alpha")
#' )
#' tg <- as_tidygenomes(genes)
#' tg <- add_patterns(tg)
#' 
#' @export
add_patterns <- function(tg) {
  
  if (is.null(tg$genes)) stop("No gene table present")
  if (is.null(tg$orthogroups)) stop("No orthogroup table present")
  
  patterns_raw <-
    tg$genes %>%
    distinct(genome, orthogroup) %>%
    mutate(present = TRUE) %>%
    spread(key = genome, value = present, fill = FALSE) %>%
    nest(orthogroup) %>%
    mutate(pattern = str_c("p", 1:n()))
  
  orthogroups_patterns <-
    patterns_raw %>%
    unnest() %>%
    select(orthogroup, pattern)
  
  tg$components <-
    patterns_raw %>%
    select(- data) %>%
    gather(key = "genome", value = "present", - pattern) %>%
    filter(present) %>%
    select(- present)
  
  tg$patterns <- 
    patterns_raw %>%
    mutate(frequency = map_int(data, ~ length(.$orthogroup))) %>%
    select(pattern, frequency)
  
  tg %>%
    modify_at("orthogroups", left_join, orthogroups_patterns, by = "orthogroup")
  
}

#' Map presence/absence patterns to the tree
#'
#' This function determines for each unique presence/absence pattern whether it
#' is a signature for a clade in the tree.
#' 
#' A variable `node` will be added to the pattern table. 
#'
#' @param tg A tidygenomes object
#' 
#' @return A tidygenomes object
#' 
#' @export
map_patterns <- function(tg) {
  
  if(is.null(tg$patterns)) stop("No patterns present")
  
  patterns_nodes <-
    tg$components %>%
    left_join(tg$genomes %>% select(genome, node), by = "genome") %>%
    group_by(pattern) %>%
    summarize(node = ancestor_if_complete(node, tg$tree))
  
  tg %>%
    modify_at("patterns", left_join, patterns_nodes, by = "pattern")
  
}

#' Add various orthogroup measures to the orthogroup table
#'
#' This function adds the variables `og_genes` and `og_genomes` to the
#' orthogroups table of a tidygenomes object.
#'
#' @param tg A tidygenomes object
#' 
#' @return A tidygenomes object
#' 
#' @export
add_orthogroup_measures <- function(tg) {
  
  if (is.null(tg$genes)) stop("No gene table present")
  if (is.null(tg$orthogroups)) stop("No orthogroup table present")
  
  orthogroups_measures <-
    tg$genes %>%
    group_by(orthogroup) %>%
    summarize(og_genes = n(), og_genomes = length(unique(genome)))
  
  tg %>%
    modify_at("orthogroups", left_join, orthogroups_measures, by = "orthogroup")
  
}

#' Add various phylogroup measures to the phylogroup table
#'
#' This function adds the variable `pg_genomes` to the orthogroups table of a
#' tidygenomes object.
#'
#' @param tg A tidygenomes object
#' 
#' @return A tidygenomes object
#' 
#' @export
add_phylogroup_measures <- function(tg) {
  
  if (is.null(tg$phylogroups)) stop("No phylogroups present")
  
  phylogroups <-
    tg$genomes %>%
    left_join(tg$nodes) %>%
    count(phylogroup) %>%
    rename(pg_genomes = n)
  
  tg %>%
    modify_at("phylogroups", left_join, phylogroups)
  
}


#' Add exclusivity of phylogroups
#'
#' This function calculates the minimum similarity within phylogroups and the
#' maximum similarity between phylogroups, assesses exclusivity and adds these
#' measures to the phylogroup table.
#'
#' @param tg A tidygenomes object
#' @param similarity An expression defining a genome similarity measure
#' 
#' @return A tidygenomes object
#' 
#' @export
add_exclusivity <- function(tg, similarity) {
  
  if (is.null(tg$phylogroups)) stop("No phylogroups present")
  
  similarity <- rlang::enexpr(similarity)
  
  genomes <-
    tg$genomes %>%
    left_join(tg$nodes, by = "node") %>%
    left_join(tg$phylogroups, by = "phylogroup") %>%
    select(genome, phylogroup)
  
  genome_pairs_full <- 
    complete_pairs(tg$pairs) %>%
    left_join(
      genomes %>% select(genome_1 = genome, phylogroup_1 = phylogroup),
      by = "genome_1"
    ) %>%
    left_join(
      genomes %>% select(genome_2 = genome, phylogroup_2 = phylogroup),
      by = "genome_2"
    ) %>%
    mutate(within = phylogroup_1 == phylogroup_2) %>%
    mutate(similarity = !! similarity)
  
  genomes_membership <-
    genome_pairs_full %>%
    rename(genome = genome_1, phylogroup = phylogroup_1) %>%
    group_by(genome, phylogroup) %>%
    group_map(function(df, groups) {
      if (sum(df$within) == 0) {
        df_within <- tibble(min_similarity_within = 1, furthest_within = groups$genome)
      } else {
        df_within <- df %>%
          filter(within) %>%
          filter(similarity == min(similarity)) %>%
          slice(1) %>%
          select(min_similarity_within = similarity, furthest_within = genome_2)
      }
      df_between <- df %>%
        filter(! within) %>%
        filter(similarity == max(similarity)) %>%
        slice(1) %>%
        select(max_similarity_between = similarity, closest_between = genome_2)
      bind_cols(df_within, df_between)
    }) %>%
    ungroup() %>%
    mutate(
      consensus_phylogroup_member = 
        min_similarity_within > max_similarity_between
    )
  
  phylogroups_exclusivity <-
    genomes_membership %>%
    group_by(phylogroup) %>%
    summarize(
      exclusive = all(consensus_phylogroup_member),
      min_similarity_within = min(min_similarity_within),
      max_similarity_between = max(max_similarity_between)
    )
  
  genomes_membership <-
    genomes_membership %>% 
    select(
      genome, 
      min_similarity_within, furthest_within,
      max_similarity_between, closest_between,
      consensus_phylogroup_member
    )
  
  tg$genomes <- left_join(tg$genomes, genomes_membership, by = "genome")
  tg$phylogroups <- 
    left_join(tg$phylogroups, phylogroups_exclusivity, by = "phylogroup")
  
  tg
  
}