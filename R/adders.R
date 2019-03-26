# These functions add external data to the tog object. 

create_tidyorthogroups <- function(orthofinder_dir = NULL, clade = NULL) {
  
  if (is.null(orthofinder_dir)) {
    
    tog <- list(
      genes = tibble(),
      genomes = tibble(),
      orthogroups = tibble() 
    )
    
    class(tog) <- "tidyorthogroups"
    
  } else {
    
    ! is.null(clade) || stop("please supply a clade name")
    
    tog <- create_tidyorthogroups()
    tog <- add_orthogroups_clade(tog, orthofinder_dir, clade)

  }
  
  tog
  
}

add_genome_tibble <- function(tog, genome_tibble) {
  
  modify_at(tog, "genomes", left_join, genome_tibble)
  
}

add_orthogroup_tibble <- function(tog, orthogroup_tibble) {
  
  modify_at(tog, "orthogroups", left_join, orthogroup_tibble)
  
}

add_n_genomes_in_orthogroup <- function(tog) {
  
  orthogroups_n_genomes <- 
    tog$genes %>%
    distinct(genome, orthogroup) %>%
    count(orthogroup) %>%
    rename(n_genomes = n)
  
  add_orthogroup_tibble(tog, orthogroups_n_genomes)
  
}

add_orthogroups_species <- function(tog, orthofinder_dir, clade) {
  
  genes_new <- read_orthogroups(orthofinder_dir)
  
  genomes_new <- tibble(
    genome = genes_new$genome %>% unique(),
    clade = !! clade
  ) 
  
  orthogroups_new <- tibble(
    orthogroup_name = genes_new$orthogroup %>% unique(),
    orthogroup = 1:length(orthogroup_name) + nrow(tog$orthogroups),
    clade = !! clade
  )
  
  # change orthogroup names by indices defined in table orthogroups_new
  genes_new <- genes_new %>%
    left_join(select(orthogroups_new, orthogroup_index = orthogroup, orthogroup = orthogroup_name)) %>%
    select(gene, genome, orthogroup = orthogroup_index)
  
  tog$genes <- bind_rows(tog$genes, genes_new)
  tog$genomes <- bind_rows(tog$genomes, genomes_new)
  tog$orthogroups <- bind_rows(tog$orthogroups, orthogroups_new)
  
  tog
  
}

add_metaorthogroups <- function(tog, orthofinder_dir) {
  
  orthogroups <- read_orthogroups(orthofinder_dir) %>%
    select(orthogroup_name = gene, metaorthogroup = orthogroup)

  tog$orthogroups <- tog$orthogroups %>%
    left_join(orthogroups)
  
  tog
  
}

add_genome_data <- function(tog, genomes) {
  
  tog$genomes <- left_join(tog$genomes, genomes)
  
  tog
  
}

left_join2 <- function(df1, df2) {
  
  df1 <- keep(df1, ~ ! every(., is.na))
  left_join(df1, df2)
  
}

add_bedfiles <- function(tog, bedfiles) {

  tog$genes <- split(tog$genes, tog$genes$genome) 
  
  for (bedfile in bedfiles) {
    
    bedfile_genes <- readr::read_tsv(file = bedfile, col_names = F) %>%
      `names<-`(c("contig", "start", "end", "gene", "score", "strand")) %>%
      mutate_at("gene", str_c, "_1", sep = "")
    
    genome <- str_match(bedfile, "([^/]+)\\.bed")[, 2]

    tog$genes <- modify_at(tog$genes, genome, left_join2, bedfile_genes) 
    
  }

  tog$genes <- bind_rows(tog$genes)
  
  tog
  
}

add_gff_files <- function(tog, gff_files) {
  
  tog$genes <- split(tog$genes, tog$genes$genome) 
  
  for (gff_file in gff_files) {
    
    gff_file_genes <- readr::read_tsv(file = gff_file, col_names = F, comment = "#") %>%
      `names<-`(c("contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
      mutate(gene = str_match(attributes, "ID=[0-9]+_([^;]*);")[,2]) %>%
      mutate(gene = str_c(contig, gene, sep = "_"))
    
    genome <- str_match(gff_file, "([^/]+)\\.gff")[, 2]
    
    tog$genes[[genome]] <- left_join(tog$genes[[genome]], gff_file_genes)
    
  }
  
  tog$genes <- bind_rows(tog$genes)
  
  tog
  
}
