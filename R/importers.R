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
  
  c(tg, data) %>%
    `class<-`("tidygenomes")
  
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
#' The phylogroup table should consist of one row per phylogroup and should
#' contain the following variables:
#' 
#' * phylogroup: the name of the phylogroup
#' * genome_type: the type genome of the phylogroup
#' * genome_peripheral: a peripheral genome of the phylogroup (optional)
#' 
#' If peripheral genomes are supplied, a phylogroup is defined as all
#' descendants of the most recent common ancestor of the type and peripheral
#' genomes. If peripheral genomes are not supplied, a phylogroup is defined as
#' all descendants of the oldest ancestor of the type genome that is not also an
#' ancestor of any other type genomes. 
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
  
  # convert genome identifiers to node names
  lut_nodes <-
    tg$genomes %>%
    mutate(genome_identifier = !! genome_identifier) %>%
    {structure(.$node, names = (.$genome_identifier))}
  phylogroups <- 
    phylogroups %>% mutate(genome_type = unname(lut_nodes[genome_type])) 
  
  # determine the phylogroup of all nodes and whether they are ancestral nodes
  # of phylogroups
  if ("genome_peripheral" %in% names(phylogroups)) {
    
    # situation 1: phylogroups defined by type and peripheral genomes
    
    phylogroups <- 
      phylogroups %>% 
      mutate(genome_peripheral = unname(lut_nodes[genome_peripheral]))
    nodes_phylogroups <-
      phylogroups %>%
      mutate(anc_node = map2_chr(
        genome_type, genome_peripheral, 
        ~ mrca(c(.x, .y), tree = tg$tree)
      )) %>%
      mutate(node = map(anc_node, ~ descendants(., tree = tg$tree))) %>%
      unnest(node) %>%
      mutate(is_phylogroup_ancestor = node == anc_node) %>%
      select(node, phylogroup, is_phylogroup_ancestor)
    
    phylogroups <- phylogroups %>% select(- genome_peripheral)
    
  } else {
    
    # situation 2: phylogroups defined by types only 
    
    # put known phylogroups in numeric vector, where numbers represent node
    # numbers
    nodes_phylogroups <- 
      tg$nodes %>%
      left_join(phylogroups, by = c("node" = "genome_type")) %>%
      right_join(
        tibble(node = c(tg$tree$tip.label, tg$tree$node.label)),
        by = "node"
      ) %>%
      pull(phylogroup)
    
    # initialize empty vector with node numbers of phylogroup ancestors
    ancestral_nodes <- numeric()
    
    # children-first traversal: infer phylogroups from children
    tree <- tg$tree %>% ape::reorder.phylo("postorder") 
    parents <- unique(tree$edge[, 1])
    for (parent in parents) {
      children <- tree$edge[tree$edge[, 1] == parent, 2]
      children_phylogroup <- 
        nodes_phylogroups[children] %>%
        {.[! is.na(.)]}
      if (length(children_phylogroup) == 1) {
        nodes_phylogroups[parent] <- children_phylogroup
      } else if (length(children_phylogroup) == 2) {
        nodes_phylogroups[parent] <- "no phylogroup"
        ancestral_nodes <- c(ancestral_nodes, children)
      }
    }
    
    # remove ancestral nodes with the label "no phylogroup"
    ancestral_nodes <- 
      ancestral_nodes[nodes_phylogroups[ancestral_nodes] != "no phylogroup"]
    
    # parents-first traversal: infer phylogroups from parents
    root <- tree$edge[nrow(tree$edge), 1]
    nodes_phylogroups[root] <- "no phylogroup"
    for (child_row in nrow(tree$edge):1) {
      child <- tree$edge[child_row, 2]
      if (is.na(nodes_phylogroups[child])) {
        parent <- tree$edge[child_row, 1]
        nodes_phylogroups[child] <- nodes_phylogroups[parent]
      }
    }
    
    # make node table with their phylogroups
    nodes_phylogroups <-
      tibble(
        node_number = 1:length(nodes_phylogroups), 
        phylogroup = nodes_phylogroups
      ) %>%
      mutate(
        node = c(tree$tip.label, tree$node.label)[node_number],
        is_phylogroup_ancestor = node_number %in% ancestral_nodes
      ) %>%
      select(- node_number)
    
  }
  
  # stop if nodes are not unique (phylogroups overlap)
  is_unique <- function(x) length(unique(x)) == length(x)
  if (! is_unique(nodes_phylogroups$node)) {
    stop("Phylogroups overlap")
  } 
  
  # add phylogroups to node table
  tg$nodes <- 
    tg$nodes %>% 
    left_join(nodes_phylogroups, by = "node") %>%
    replace_na(list(phylogroup = "no phylogroup"))
  
  # add to genome table whether genomes are type of a phylogroup
  tg$genomes <- 
    tg$genomes %>%
    mutate(is_phylogroup_type = node %in% !! phylogroups$genome_type)
  
  # add the phylogroup table
  tg$phylogroups <- phylogroups %>% select(- genome_type)
  
  tg
  
}

#' Remove phylogroups from a tidygenomes object
#'
#' This function undoes the effects of the function [add_phylogroups()]
#'
#' @param tg A tidygenomes object
#' 
#' @return A tidygenomes object
#' 
#' @export
remove_phylogroups <- function(tg) {
  
  tg$phylogroups <- NULL
  tg$nodes$phylogroup <- NULL
  tg$nodes$is_phylogroup_ancestor <- NULL
  tg$genomes$is_phylogroup_type <- NULL
  
  tg
  
}

#' Inflate one of the pangenomes of a metapangenome
#'
#' This function inflates one of the pangenomes of what you could call a
#' "metapangenome": a pangenome where one or more individual genomes
#' represent(s) the entire pangenome of one species. The orthogroups of this
#' species are represented as individual genes in the metapangenome. This
#' function will replace the single species to be inflated by the individual
#' genomes of which is consists, and its "genes" with the complete species-level
#' orthogroups.
#'
#' The genes in the metapangenome that belong to the genome (species) to be
#' inflated should correspond to orthogroups in the species pangenome.
#'
#' @param tg_meta A tidygenomes object containing the metapangenome
#' @param tg_species A tidgenomes object containing the species pangenome
#' @param species The genome in the metapangenome that represents the species to
#'   inflate
#'
#' @return A tidygenomes object
#'
#' @export
inflate_pangenome <- function(tg_meta, tg_species, species) {
  
  if (! is.null(tg_meta$tree)) {
    stop("Remove the tree from the metapangenome tidygenomes object")
  }
  
  if (! is.null(tg_meta$phylogroups)) {
    stop("Remove the phylogroups and tree from the metapangenome ", 
         "tidygenomes object")
  }
  
  if (! is.null(tg_meta$pairs)) {
    stop("Remove the genome pairs from the metapangenome tidygenomes object")
  }
  
  if (! all(tg_species$orthogroups$orthogroup %in% tg_meta$genes$gene)) {
    stop("Make sure that all species orthogroup are present as genes in the ",
         "metapangenome gene table")
  }
  
  # split metapangenome gene table in genes to inflate and genes not to inflate
  tg_meta$genes <-
    tg_meta$genes %>%
    mutate(split = if_else(genome == !! species, "inflate", "notinflate")) %>%
    split(f = .$split)
  
  # use the metapangenomes genes to inflate to give new orthogroup names to the
  # species genes
  tg_meta$genes$inflate <- 
    tg_meta$genes$inflate %>%
    select(orthogroup, orthogroup_species = gene)
  tg_species$genes <-
    tg_species$genes %>%
    rename(orthogroup_species = orthogroup) %>%
    left_join(tg_meta$genes$inflate, by = "orthogroup_species") %>%
    select(- orthogroup_species)
  
  # add the updated species genes to the metapangenome genes not to inflate
  tg_meta$genes$inflate <- tg_species$genes 
  tg_meta$genes$notinflate$split <- NULL
  tg_meta$genes <- bind_rows(tg_meta$genes$inflate, tg_meta$genes$notinflate)
  
  # inflate the genome table 
  tg_meta$genomes <-
    tg_meta$genomes %>%
    filter(genome != !! species) %>%
    bind_rows(tg_species$genomes)
  
  tg_meta
  
}