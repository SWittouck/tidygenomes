# if one metric given: cladevar on x-axis and metric1 on y-axis
# if two metrics given: metric1 on x-axis, metric1 on y-axis and facet wrap on cladevar
# to do: make the "put sp. last" functionality more general than Lactobacillus
# to do: make the color flexible (e.g. enable coloring of different variable than assembly level)
plot_genomes_per_species <- function(genomes, metric1, metric2, cladevar = clade_name, col = NULL, n_min = 5) {
  
  metric1_quosure <- enquo(metric1)
  metric1_formula <- substitute(metric1)
  
  metric2_quosure <- enquo(metric2)
  metric2_formula <- substitute(metric2)
  
  cladevar_quosure <- enquo(cladevar)
  cladevar_formula <- substitute(cladevar)
  
  col_quosure <- enquo(col)
  col_formula <- substitute(col)
  
  if (missing(metric2)) {
    
    clade_name_ordered <- genomes %>%
      group_by(!! cladevar_quosure) %>%
      summarize(med_genes = median(!! metric1_quosure)) %>%
      arrange(- med_genes) %>%
      pull(!! cladevar_quosure) %>%
      setdiff("Lactobacillus sp.") %>%
      union("Lactobacillus sp.")
    
    genomes <- genomes %>%
      mutate(!! cladevar_formula := factor(!! cladevar_quosure, levels = !! clade_name_ordered))
    
    genomes %>%
      add_count(!! cladevar_quosure) %>%
      filter(n >= !! n_min) %>%
      ggplot(aes_(x = cladevar_formula, y = metric1_formula, group = cladevar_formula, col = col_formula)) +
      geom_boxplot(col = "grey", outlier.alpha = 0) +
      geom_jitter(width = 0.2, height = 0) +
      scale_color_brewer(palette = "Paired") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0.1, 0.1, 0.1, 1.2), "cm")
      ) 
    
  } else {
    
    ggplot(genomes, aes_(x = metric1_formula, y = metric2_formula, col = col_formula)) +
      geom_point() +
      scale_color_brewer(palette = "Paired") +
      facet_wrap(cladevar_formula, scales = "free")
    
  }
  
}

upset_plot <- function(tog, genome_name, genome_col, genome_bold, n = 50) {
    
  genome_name <- rlang::enexpr(genome_name) 
  genome_col <- rlang::enexpr(genome_col)
  genome_bold <- rlang::enexpr(genome_bold)
  
  tidypatterns <- 
    tidypatterns(tog) %>%
    modify_at("patterns", arrange, desc(n)) %>%
    modify_at("patterns", mutate, rank = 1:n()) %>%
    filter_patterns(rank <= !! n) %>%
    modify_at("patterns", mutate, pattern_ordered = factor(pattern, levels = pattern))
  
  theme_upset <- 
    theme(
      axis.text.x = element_blank(), 
      axis.title = element_blank(),
      axis.ticks = element_blank(), 
      panel.background = element_rect(fill = "transparent"), 
      plot.background = element_rect(fill = "transparent")
  )
  
  tog$genomes <-
    tog$genomes %>%
    arrange(!! genome_name) %>%
    group_by(!! genome_col) %>%
    nest() %>%
    mutate(genome_col = rep_len(as.character(1:12), nrow(.))) %>%
    unnest() %>%
    ungroup()
  
  tog$genomes <- rename(tog$genomes, genome_bold = !! genome_bold)
  
  plot_main <- 
    tidypatterns$components %>%
    left_join(tidypatterns$patterns) %>%
    left_join(tog$genomes) %>%
    ggplot(aes(x = pattern_ordered, y = !! genome_name, group = pattern_ordered)) +
    geom_line(col = "grey") +
    geom_point(size = 3, aes(col = genome_col)) +
    scale_color_brewer(palette = "Paired", guide = "none") +
    theme_bw() +
    theme_upset +
    theme(axis.text.y = element_text(face = if_else(tog$genomes$genome_bold, "bold", "plain")))

  plot_marg <- 
    tidypatterns$patterns %>% 
    ggplot(aes(x = pattern_ordered, y = n)) +
    geom_histogram(stat = "identity") +
    # geom_text(aes(label = n), vjust= - 0.25, size = 2) + 
    ylim(c(0, 1.1 * max(tidypatterns$patterns$n))) +
    theme_bw() +
    theme_upset +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
  
  egg::ggarrange(
    plot_marg, plot_main,
    nrow = 2, ncol = 1, heights = c(1, nrow(lgc$genomes) / 50)
  )
  
}
