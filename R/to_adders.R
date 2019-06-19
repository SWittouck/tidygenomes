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
