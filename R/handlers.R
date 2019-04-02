filter_orthogroups <- function(tog, ...) {
  
  tog %>%
    modify_at("orthogroups", filter, ...) %>%
    modify_at("genes", filter, orthogroup %in% .$orthogroups$orthogroup)
  
}

filter_patterns <- function(tp, ...) {
  
  tp %>%
    modify_at("patterns", filter, ...) %>%
    modify_at("components", filter, pattern %in% .$patterns$pattern)
  
}

mutate_genomes <- function(tog, ...) {
  
  modify_at(tog, "genomes", mutate, ...)
  
}

filter_genomes <- function(to, ...) {
  
  to %>%
    purrr::modify_at("genomes", filter, ...) %>%
    purrr::modify_at("genes", filter, genome %in% .$genomes$genome) %>%
    purrr::modify_at("orthogroups", filter, orthogroup %in% .$genes$orthogroup)
  
}