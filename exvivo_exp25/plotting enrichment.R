library(dplyr)
library(ggplot2)
library(gridExtra)

## enrichment_df needs a pathway column, adj_p_val column, and enrichment column (NES from GSEA). Enrichment is expected to be a continuous variable.
## top is the number of pathways (up regulated and down regulated) to take.

plot_enrichment_result <- function(enrichment_result, enrichment_label, top = 10,
                                   adj_pval_cutoff = 0.05, enrich_width = 9, sig_width = 4,
                                   enrichment_title = "Enrichment"){
   
   #' Used to make reversed logarithmic scales
   #' @import scales
   reverselog_trans <- function(base = exp(1)) {
      library(scales)
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                        breaks = scales::log_breaks(base = base),
                        domain = c(1e-100, Inf))
   }
   
   
   plot_df <- enrichment_result %>%
      dplyr::select(pathway, adj_p_val, enrichment) %>%
      filter(adj_p_val < adj_pval_cutoff) %>%
      mutate(sign = case_when(enrichment > 0 ~ "Up",
                              TRUE ~ "Down")) %>%
      group_by(sign) %>%
      arrange(adj_p_val) %>%
      slice(1:top) %>%
      ungroup()
   
   p_enrichment <- ggplot(plot_df, aes(x = enrichment, y = reorder(pathway, enrichment))) +
      geom_bar(stat = 'identity', aes(fill = sign)) +
      scale_fill_manual(values = c("Down" = "dodgerblue3", "Up" = "firebrick2", "Not significant" = "black")) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
            legend.position = "none",
            axis.title.x = element_text(size=16),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 12),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(x = enrichment_label) +
      ggtitle(enrichment_title)
   
   
   p_sig <- ggplot(plot_df, aes(x = adj_p_val, y = reorder(pathway, enrichment))) +
      geom_bar(stat='identity') +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
            legend.position="none",
            axis.title.x = element_text(size=16),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.ticks.y = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.y = element_blank()) +
      scale_x_continuous(trans = reverselog_trans(10)) +
      labs(x = "Adjusted p-value") +
      ggtitle("Significance")
   
   arrange_matrix <- as.matrix(c(rep(1, enrich_width), rep(2, sig_width))) %>% t()
   p_both <- grid.arrange(p_enrichment, p_sig, layout_matrix = arrange_matrix)
}