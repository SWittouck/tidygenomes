#' Perform postorder tree traversal
#'
#' This functions returns the node numbers of a tree in order of postorder tree
#' traversal (tips first).
#'
#' @param tree An object of the class phylo
#'
#' @return A numeric vector with node numbers
#'
#' @export
nodes_postorder <- 
  function(tree, node = setdiff(tree$edge[, 1], tree$edge[, 2]), nodes = c()) {
    # note on the implementation:
    # - it is recursive 
    # - the node argument starts on the root node 
    # - the nodes argument starts as an empty vector
    children <- tree$edge %>% {.[.[, 1] == node, 2]}
    for (child in children) {
      nodes <- nodes_postorder(tree, child, nodes) 
    }
    nodes <- c(nodes, node)
    nodes
  }

#' Perform preorder tree traversal
#'
#' This functions returns the node numbers of a tree in order of preorder tree
#' traversal (root first).
#'
#' @param tree An object of the class phylo
#'
#' @return A numeric vector with node numbers
#'
#' @export
nodes_preorder <- 
  function(tree, node = setdiff(tree$edge[, 1], tree$edge[, 2]), nodes = c()) {
    # note on the implementation:
    # - it is recursive 
    # - the node argument starts on the root node 
    # - the nodes argument starts as an empty vector
    nodes <- c(nodes, node)
    children <- tree$edge %>% {.[.[, 1] == node, 2]}
    for (child in children) {
      nodes <- nodes_preorder(tree, child, nodes) 
    }
    nodes
  }

#' Complement the pairs of a pair table
#'
#' For each unique pair of objects (a, b), the function adds a row for the pair
#' (b, a) with the same information.
#' 
#' @param pairs Data frame where each row represents a pair of objects
#' @param object_1 Variable defining the first object (unquoted)
#' @param object_2 Variable defining the second object (unquoted)
#' 
#' @return A complemented table
#' 
#' @export
complete_pairs <- function(pairs, object_1, object_2) {
  
  object_1 <- rlang::enexpr(object_1)
  object_2 <- rlang::enexpr(object_2)
  
  pairs_2 <-
    pairs %>%
    rename(object_1_new = !! object_2, object_2_new = !! object_1) %>%
    rename(!! object_1 := object_1_new, !! object_2 := object_2_new)
  
  bind_rows(pairs, pairs_2)
  
}

#' Convert dist object to pair table
#'
#' This functions converts an object of the class "dist" to a tidy table with
#' the columns object_1, object_2 and distance.
#' 
#' @param dist An object of the class "dist"
#' 
#' @return A tibble
#' 
#' @export
dist2pairs <- function(dist) {
  
  if (! "dist" %in% class(dist)) stop("data is not a dist object")
  
  n <- attr(dist, "Size")
  
  expand.grid(i = 1:n, j = 1:n) %>%
    filter(i < j) %>%
    mutate(object_1 = labels(dist)[i]) %>%
    mutate(object_2 = labels(dist)[j]) %>%
    mutate(distance = dist[n * (i - 1) - i * (i - 1) / 2 + j - i]) %>%
    select(- i, - j)
  
}

#' Convert pair table to matrix
#'
#' This functions converts a data frame where each row is a unique pairwise
#' comparison to a (dis)similarity matrix.
#' 
#' @param pairs A table with pairwise comparisons
#' @param object_1 Name of the variable specifying the first object of the pair
#' @param object_2 Name of the variable specifying the second object of the pair
#' @param measure Name of the (dis)similarity variable
#' @param diag Something to fill the diagonal with
#' 
#' @return A matrix
#' 
#' @export
pairs2matrix <- 
  function(pairs, object_1 = object_1, object_2 = object_2, measure, diag) {
    
    object_1 <- rlang::enquo(object_1)
    object_2 <- rlang::enquo(object_2)
    measure <- rlang::enquo(measure)
    
    pairs %>%
      transmute(
        object_1 = {{object_1}}, object_2 = {{object_2}}, measure = {{measure}}
      ) %>% 
      complete_pairs(object_1, object_2) %>%
      pivot_wider(names_from = object_2, values_from = measure) %>%
      `class<-`("data.frame") %>%
      `rownames<-`(.$object_1) %>%
      select(- object_1) %>%
      as.matrix() %>%
      {.[rownames(.), rownames(.)]} %>%
      `diag<-`(diag)
    
  }

#' Fix the pair order of a genome pair table
#'
#' For each row in a genome pair table (table with at least the columns genome_1
#' and genome_2), this function assesses whether genome_1 < genome_2. If not, it
#' swaps the values for genome_1 and genome_2, as well as for other variables
#' whose names end in "_1" and "_2".
#' 
#' @param pairs A table with at least the variables genome_1 and genome_2
#' 
#' @return A tibble
#' 
#' @export
fix_pair_order <- function(pairs) {
  
  namechange1 <- c("_1" = "_1prev", "_2" = "_2prev")
  namechange2 <- c("_1prev" = "_2", "_2prev" = "_1")
  bind_rows(
    pairs[pairs$genome_1 <= pairs$genome_2, ],
    pairs[pairs$genome_1 > pairs$genome_2, ] %>%
      rename_with(~ stringr::str_replace_all(., namechange1)) %>%
      rename_with(~ stringr::str_replace_all(., namechange2))
  )
  
}

#' Root tree given three tips
#'
#' This function roots a phylogenetic tree given three tip labels.
#'
#' Tips a, b and c define exactly one internal node in the unrooted tree. The
#' tree will be rooted on the branch leading from this node to tip a.
#' 
#' The branch labels (e.g. support values) are handled correctly, i.e. moved to
#' the correct nodes. When the tree is rooted on a non-tip branch, its branch
#' label will be copied to both child notes of the root node.
#' 
#' @param tree An object of class phylo
#' @param tips Three tip labels 
#' 
#' @return An object of class phylo
#' 
#' @export
root_tree.phylo <- function(tree, tips) {
  
  if (! "phylo" %in% class(tree)) {
    stop("tree should be of class phylo")
  }
  
  # root on node defined by tips
  root_new <- 
    ape::mrca(tree) %>%
    {.[tips, tips]} %>%
    {.[. > length(tree$tip.label)]} %>%
    {.[. == max(.)]} %>%
    {.[1]}
  tree <- tree %>% ape::root.phylo(node = root_new, edgelabel = T)
  
  # resolve root node such that first tip is (part of) outgroup
  outgroup <-
    ape::mrca(tree) %>%
    {.[tips[1], ]} %>%
    {.[. != length(tree$tip.label) + 1]} %>%
    names()
  tree <- 
    tree %>% 
    ape::root.phylo(outgroup = outgroup, edgelabel = T, resolve.root = T)
  
  # divide root branch length
  if ("edge.length" %in% names(tree)) {
    n_tips <- length(tree$tip.label) 
    l <- sum(tree$edge.length[tree$edge[, 1] == n_tips + 1])
    tree$edge.length[tree$edge[, 1] == n_tips + 1] <- l / 2
  }
  
  # copy the branch support of one rootchild to the other
  # except when one of them is a tip! 
  # (phangorn midpoint rooting also copies the support value of the rootchilds)
  n_tips <- length(tree$tip.label)
  children <- tree %>% {.$edge[.$edge[, 1] == n_tips + 1, 2]}
  if (all(children > n_tips)) {
    label <- c(tree$tip.label, tree$node.label)[children] %>% keep(~ . != "")
    tree$node.label[children - n_tips] <- label
  }
  
  # return tree
  tree
  
}

#' Root tree given three genomes
#'
#' This applies [root_tree.phylo] to the tree component of a tidygenomes object.
#' 
#' @param tg An tidygenomes object
#' @param root Three tips that identify the root (see [root_tree.phylo])
#' @param genome_identifier Variable of the genome table that corresponds to the
#'   given genomes
#' 
#' @return A tidygenomes object
#' 
#' @export
root_tree <- function(tg, root, genome_identifier = genome) {
  
  if (! "tree" %in% names(tg)) {
    stop("Tg should contain a tree")
  }
  
  genome_identifier <- rlang::enexpr(genome_identifier)
  
  genome_to_node <-
    tg$genomes %>%
    mutate(genome_identifier = !! genome_identifier) %>%
    {structure(.$node, names = .$genome_identifier)}
  
  tips <- genome_to_node[root]
  
  if (! length(tips) == 3) {
    stop("Not all nodes were found")
  }
  
  tg$tree <- tg$tree %>% root_tree.phylo(tips) 
  
  tg
  
}

#' Add a branch to the root of the tree
#' 
#' This function adds a branch to the root of the tree, adding one extra node in
#' the process.
#' 
#' @param tree A phylogeny of class `phylo`, with edge lengths
#' @param branch_length The length of the root branch
#' 
#' @return An object of class `phylo`
#' 
#' @example 
#' tree <- ape::read.tree(text = "(a:1, b:1);")
#' tree <- add_rootbranch(tree)
#' plot(tree)
#' 
#' @export
add_rootbranch <- function(tree, branch_length = 1) {
  
  n_tips <- length(tree$tip.label)
  tree$edge[tree$edge > n_tips] <- tree$edge[tree$edge > n_tips] + 1
  tree$edge <- rbind(c(n_tips + 1, n_tips + 2), tree$edge)
  tree$edge.length <- c(branch_length, tree$edge.length)
  tree$Nnode <- tree$Nnode + 1
  
  tree
  
}

mrca <- function(tips, tree) {
  
  tips <- 
    tips %>%
    unique() %>%
    {.[. %in% c(tree$tip.label, tree$node.label)]}
  
  if (length(tips) == 0) stop("tips not found in tree")
  
  if (length(tips) == 1) return(tips)
  
  tips %>%
    phytools::findMRCA(tree = tree) %>%
    {c(tree$tip.label, tree$node.label)[.]}
  
}

# descendantes include node itself! 
descendants <- function(node, tree) {
  
  node %>%
    {which(c(tree$tip.label, tree$node.label) == .)} %>%
    phytools::getDescendants(tree = tree) %>%
    {c(tree$tip.label, tree$node.label)[.]} %>%
    {c(., node)} %>%
    unique()
  
}

ancestor_if_complete <- function(tip_labels, tree) {
  
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

genomes_extended <- function(tg) {
  
  if(tibble::has_name(tg, "nodes")) {
    tg$genomes <- tg$genomes %>% left_join(tg$nodes, by = "node")
  }
  
  if(tibble::has_name(tg, "phylogroups")) {
    tg$genomes <- tg$genomes %>% left_join(tg$phylogroups, by = "phylogroup")
  }
  
  tg$genomes
  
}

tipnodes_ladderized <- function(tree) {
  
  tree <- ape::ladderize(tree)
  
  tree$edge[, 2] %>%
  {.[. <= length(tree$tip.label)]} %>%
  {tree$tip.label[.]}
  
}