#' Return a pangenome matrix
#'
#' This function returns the orthogroup content of the genomes in the form of a
#' matrix where the rows represent genomes, the columns represent orthogroups
#' and the cells contain the copy number or presence status of the orthogroups
#' in the genomes.
#'
#' @param tg A tidygenomes object
#' 
#' @return A matrix with the pangenome
#' 
#' @export
pangenome_matrix <- function(tg) {
  
  tg$genes %>%
    count(genome, orthogroup) %>%
    spread(key = orthogroup, value = n, fill = 0) %>%
    `class<-`("data.frame") %>%
    `rownames<-`(.$genome) %>%
    select(- genome) %>%
    as.matrix()
  
}