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
    left_join(tg$nodes, by = "node") %>%
    count(phylogroup) %>%
    rename(pg_genomes = n)
  
  tg %>%
    modify_at("phylogroups", left_join, phylogroups, by = "phylogroup")
  
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

#' Add color variable for phylogroups
#'
#' This function adds a variable phylogroup_color to the phylogroup table that
#' can be used to plot phylogroup colors in cases where there are too many
#' phylogroups for the color scale you want to use.
#' 
#' The variable phylogroup_color assigns a number to each phylogroup, in order
#' of the ladderized tree. If there are more phylogroups than numbers, the
#' numbers will be recycled. The goal is that two phylogroups that are directly
#' next to each other in the tree labels will never have the same color.
#'
#' @param tg A tidygenomes object
#' @param n The number of colors
#' 
#' @return A tidygenomes object
#' 
#' @export
add_phylogroup_color <- function(tg, n = 12) {
  
  if (is.null(tg$phylogroups)) stop("No phylogroups present")
  
  n_phylogroups <- nrow(tg$phylogroups)
  
  tipnodes_ladderized <- tipnodes_ladderized(tg$tree)
  
  phylogroups_ordered <-
    tibble(node = tipnodes_ladderized) %>%
    left_join(tg$genomes, by = "node") %>%
    filter(is_phylogroup_type) %>%
    left_join(tg$nodes, by = "node") %>%
    pull(phylogroup)

  tg$phylogroups <-
    tg$phylogroups %>%
    mutate(
      phylogroup_fct = factor(phylogroup, levels = !! phylogroups_ordered)
    ) %>%
    arrange(phylogroup_fct) %>%
    mutate(
      phylogroup_color = 
        1:UQ(n) %>% as.character() %>% rep_len(!! n_phylogroups)
    ) %>%
    select(- phylogroup_fct)
  
  tg
  
} 

#' Add gene content distance
#'
#' This function adds a variable gcd (gene content distance) to the pair table.
#' As the name suggests, these are distances between genomes based on shared
#' gene content. Various types of distances can be computed; the function uses
#' the vegan package for this.
#'
#' @param tg A tidygenomes object
#' @param method See [vegan::vegdist()]
#' @param binary See [vegan::vegdist()]
#' 
#' @return A tidygenomes object
#' 
#' @export
add_gcd <- function(tg, method, binary) {
  
  pairs_gcd <-
    tg %>%
    pangenome_matrix() %>%
    vegan::vegdist(method = method, binary = binary) %>%
    as_tibble() %>%
    rename(genome_2 = object_1, genome_1 = object_2, gcd = distance)
  
  tg %>%
    modify_at("pairs", left_join, pairs_gcd, by = c("genome_1", "genome_2"))
  
}

#' Collapse species pangenomes
#'
#' This function collapses each species to a single genome, and each orthogroup
#' within a species to a single gene.
#' 
#' A variable "species" must be present in the genome table for this function to
#' work.
#' 
#' A variable "status" will be added to the gene table, indicating for each gene
#' whether the orthogroup that gene belongs to is core or accessory within the
#' species.
#'
#' @param tg A tidygenomes object
#' @param core_threshold The minimum percentage of genomes within a species
#'   where a gene needs to be present to be considered a core gene
#' 
#' @return A tidygenomes object
#' 
#' @export
collapse_species <- function(tg, core_threshold = 0.9) {
  
  if (! is.null(tg$tree) | ! is.null(tg$pairs)) {
    stop("Remove the tree and/or genome pairs before collapsing species")
  }
  
  if (is.null(tg$genomes$species)) {
    stop("Add a variable 'species' to the genome table")
  }
  
  genomes <-
    tg$genomes %>%
    add_count(species) %>%
    rename(sp_genomes = n)
  
  tg$genes <-
    tg$genes %>%
    left_join(genomes, by = "genome") %>%
    distinct(orthogroup, genome, species, sp_genomes) %>%
    group_by(orthogroup, species) %>%
    summarize(status = if_else(
      n() >= sp_genomes[1] * !! core_threshold,
      "core", "accessory"
    )) %>%
    ungroup() %>%
    rename(genome = species) %>%
    mutate(gene = as.character(1:n()))
  
  tg$genomes <-
    genomes %>%
    select(genome = species, sp_genomes) %>%
    distinct(genome, sp_genomes)
  
  tg
  
}