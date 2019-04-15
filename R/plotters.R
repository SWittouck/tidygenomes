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

upset_plot <- function(
  tog, genome_name, genome_col, genome_bold, 
  color_scale = scale_color_brewer(palette = "Paired", guide = "none"), 
  n = 50) {
  
  genome_name <- rlang::enexpr(genome_name) 
  genome_col <- rlang::enexpr(genome_col)
  genome_bold <- rlang::enexpr(genome_bold)
  
  if (is.null(tog$patterns)) stop("no patterns found")
  
  tog$patterns <-
    tog$patterns %>%
    arrange(desc(n_orthogroups)) %>%
    slice(1:UQ(n)) %>%
    mutate(pattern_fct = factor(pattern, levels = pattern)) 
  
  tog$genomes <- 
    tog$genomes %>%
    mutate(genome_fct = factor(
      genome, levels = !! genomes_ladderized(tog$tree)
    )) %>%
    arrange(desc(genome_fct)) %>%
    mutate(genome_name_fct = factor(
      !! genome_name, levels = !! genome_name
    )) %>%
    rename(genome_bold = !! genome_bold)
  
  theme_upset <- 
    theme(
      axis.text.x = element_blank(), 
      axis.title = element_blank(),
      axis.ticks = element_blank(), 
      panel.background = element_rect(fill = "transparent"), 
      plot.background = element_rect(fill = "transparent")
    )
  
  plot_main <- 
    tog$components %>%
    right_join(tog$patterns) %>%
    left_join(tog$genomes) %>%
    ggplot(aes(x = pattern_fct, y = genome_name_fct, group = pattern_fct)) +
    geom_line(col = "grey") +
    geom_point(size = 3, aes(col = !! genome_col)) +
    color_scale +
    theme_bw() +
    theme_upset +
    theme(axis.text.y = element_text(face = if_else(tog$genomes$genome_bold, "bold", "plain")))
  
  plot_marg <- 
    tog$patterns %>% 
    ggplot(aes(x = pattern_fct, y = n_orthogroups)) +
    geom_histogram(stat = "identity") +
    # geom_text(aes(label = n), vjust= - 0.25, size = 2) + 
    ylim(c(0, 1.1 * max(tog$patterns$n_orthogroups))) +
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
