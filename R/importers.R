#' Create a new tidygenomes object
#'
#' This function creates a new tidygenomes object from a genome table, pangenome
#' table, genome tree or genome pair table.
#' 
#' `data` should be one of the following:
#' 
#' * a **genome table**: contains a variable `genome`, with unique values
#' * a **pangenome table**: contains the variables `gene`, `genome` and
#' `orthogroup`
#' * a **genome tree**: should be of class `phylo`
#' * a **genome pair table**: contains the variables `genome_1` and
#' `genome_2`, with only unique combinations
#'
#' @param data A data type that can be converted into a tidygenomes object
#' 
#' @return A tidygenomes object: this is a list consisting of a genome table and
#'   optionally extra data objects related to the genomes in the genome table
#'
#' @examples
#' genomes <- data.frame(genome = c("genome 1", "genome 2"))
#' tg <- as_tidygenomes(genomes)
#' 
#' @export
as_tidygenomes <- function(data) {
  
  is_genome_table <- function(data) {
    if (! is.data.frame(data)) return(FALSE)
    if (! "genome" %in% names(data)) return(FALSE)
    if (! length(unique(data$genome)) == length(data$genome)) {
      return(FALSE)
    }
    TRUE
  }
  
  is_pangenome_table <- function(data) {
    is.data.frame(data) & 
      all(c("gene", "genome", "orthogroup") %in% names(data))
  }
  
  is_genome_tree <- function(data) {
    "phylo" %in% class(data)
  }
  
  is_genome_pair_table <- function(data) {
    is.data.frame(data) &
      all(c("genome_1", "genome_2") %in% names(data))
  }
  
  if ("tidygenomes" %in% class(data)) {
    
    return(data)
    
  } else if (is_genome_table(data)) {
    
    message("Creating tidygenomes object from genome table")
    tg <- list(genomes = data)
    
  } else if (is_pangenome_table(data)) {
    
    message("Creating tidygenomes object from pangenome table")
    tg <- list(
      genomes = tibble(genome = unique(data$genome)),
      genes = data,
      orthogroups = tibble(orthogroup = unique(data$orthogroup))
    )
    
  } else if (is_genome_tree(data)) {
    
    message("Creating tidygenomes object from genome tree")
    Nnode <- data$Nnode
    Ntip <- length(data$tip.label)
    node_label <- data$node.label
    tip_label <- data$tip.label
    node_label_new <- paste0("n", 1:Nnode)
    tip_label_new <- paste0("n", (Nnode + 1):(Nnode + Ntip))
    data$node.label <- node_label_new
    data$tip.label <- tip_label_new
    nodes <- tibble(
      node = c(node_label_new, tip_label_new),
      node_label = c(node_label, rep(NA, Ntip))
    )
    genomes <- tibble(genome = tip_label, node = tip_label_new)
    tg <- list(genomes = genomes, nodes = nodes, tree = data)
    
  } else if (is_genome_pair_table(data)) {
    
    message("Creating tidygenomes object from genome pair table")
    genomes <- 
      c(data$genome_1, data$genome_2) %>%
      unique() %>%
      {tibble(genome = .)}
    tg <- list(genomes = genomes, pairs = data)
    
  } else {
    
    stop("Data is not recognized")
    
  }
  
  class(tg) <- "tidygenomes"
  
  tg
  
}

#' Add a tidygenomes object to a tidygenomes object
#'
#' This function adds a second tidygenomes object, or data from which a
#' tidygenomes object can be created, to a tidygenomes object.
#' 
#' See [as_tidygenomes()] for a list of what `data` can be. 
#'
#' @param tg A tidygenomes object
#' @param data A tidygenomes object or data type that can be converted into a
#'   tidygenomes object
#' 
#' @return A tidygenomes object
#'
#' @examples
#' genomes <- data.frame(genome = c("genome 1", "genome 2"))
#' genes <- data.frame(
#'   gene = c("gene A", "gene B", "gene C"),
#'   genome = c("genome 1", "genome 1", "genome 2"),
#'   orthogroup = c("og alpha", "og beta", "og alpha")
#' )
#' tg <- as_tidygenomes(genomes)
#' tg <- add_tidygenomes(genes)
#' 
#' @export
add_tidygenomes <- function(tg, data) {
  
  data <- as_tidygenomes(data)
  
  common_components <- intersect(names(tg), names(data))
  if (length(common_components) != 1) {
    stop("Tg and data contain overlapping components")
  }
  
  common_genomevars <- intersect(names(tg$genomes), names(data$genomes))
  if (length(common_genomevars) != 1) {
    stop("Tg and data contain overlapping genome variables")
  }
  
  if (! setequal(tg$genomes$genome, data$genomes$genome)) {
    genomes_to_keep <- intersect(tg$genomes$genome, data$genomes$genome)
    message(
      "Tg and data contain ", nrow(tg$genomes), " and ", nrow(data$genomes), 
      " genomes respectively; the overlap of ", length(genomes_to_keep), 
      " genomes was retained"
    )
    tg <- filter_genomes(tg, genome %in% !! genomes_to_keep)
    data <- filter_genomes(data, genome %in% !! genomes_to_keep)
  }
  
  tg$genomes <- tg$genomes %>% left_join(data$genomes, by = "genome")
  data$genomes <- NULL
  
  c(tg, data)
  
}

#' Add genome metadata to tidygenomes object
#'
#' This function joins a table containing genome data to the genome table of a
#' tidygenomes object.
#'
#' @param tg A tidygenomes object
#' @param genome_metadata A table having at least one variable in common with
#'   the genomes table of the tidygenomes object
#' @param by The variable to join by
#' 
#' @return A tidygenomes object
#' 
#' @export
add_genome_metadata <- function(tg, genome_metadata, by = "genome") {
  
  tg %>%
    modify_at("genomes", left_join, genome_metadata, by = by)
  
}

#' Add phylogroups to a tidygenomes object
#'
#' Phylogroups are defined as clades in a phylogenetic tree. This function will
#' add user-defined phylogroups to a tidygenomes object containing a tree.
#' 
#' The phylogroups table should consist of one row per phylogroup and should
#' contain the following variables:
#' 
#' * phylogroup: the name of the phylogroup
#' * genome_type: the type genome of the phylogroup
#' * genome_peripheral: a peripheral genome of the phylogroup
#' 
#' A phylogroup is then defined as the most recent common ancestor of the type
#' and peripheral genomes in the tree, along with all its descendants.
#' 
#' The following will be added to the tidygenomes object:
#' 
#' * A variable `phylogroup` in the nodes table.
#' * A variable `is_phylogroup_ancestor` in the nodes table. 
#' * A variable `is_phylogroup_type` in the genomes table. 
#'
#' @param tg A tidygenomes object
#' @param phylogroups A table specifying phylogroups
#' @param genome_identifier The name of a variable from the genomes table that
#'   is used to identify genomes in the phylogroups table
#' 
#' @return A tidygenomes object
#' 
#' @export
add_phylogroups <- function(tg, phylogroups, genome_identifier = genome) {
  
  genome_identifier <- rlang::enexpr(genome_identifier)
  
  lut_nodes <-
    tg$genomes %>%
    mutate(genome_identifier = !! genome_identifier) %>%
    {structure(.$node, names = (.$genome_identifier))}
  
  phylogroups <- 
    phylogroups %>%
    mutate(genome_type = lut_nodes[genome_type]) %>%
    mutate(genome_peripheral = lut_nodes[genome_peripheral]) %>%
    mutate(anc_node = map2_chr(
      genome_type, genome_peripheral, 
      ~ mrca(c(.x, .y), tree = tg$tree)
    ))
  
  nodes_phylogroups <-
    phylogroups %>%
    mutate(node = map(anc_node, ~ descendants(., tree = tg$tree))) %>%
    unnest(node) %>%
    select(node, phylogroup)
  
  is_unique <- function(x) length(unique(x)) == length(x)
  
  if (! is_unique(nodes_phylogroups$node)) {
    stop("Phylogroups overlap")
  } 
  
  tg$nodes <- 
    tg$nodes %>% 
    left_join(nodes_phylogroups, by = "node") %>%
    mutate(is_phylogroup_ancestor = node %in% !! phylogroups$anc_node) %>%
    mutate(is_phylogroup_type = node %in% !! phylogroups$genome_type)
  
  tg$phylogroups <-
    phylogroups %>%
    select(- genome_type, - genome_peripheral, - anc_node)
  
  tg
  
}