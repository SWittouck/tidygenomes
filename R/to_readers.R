# Parser for hmmer domtables. 
# Returns tidy table where each row is a hit. 
# Remark: the implementation of this function should avoid using the function
# read_table2(), since it can't deal with comment lines after the data (probably
# a bug in readr)
read_hmmer_domtbl <- function(fin) {
  
  read_lines(fin) %>%
    tibble(line = .) %>%
    filter(! str_detect(line, "^#")) %>%
    separate(line, into = str_c("X", 1:23), sep = "[ ]+", convert = T) %>%
    `names<-`(c(
      "target", "target_accession", "target_length", 
      "query", "query_accession", "query_length", 
      "e_value_full", "score_full", "bias_full",
      "i", "i_max", "c_e_value", "i_e_value", "score", "bias", 
      "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
      "acc", "target_description"
    ))
  
}
