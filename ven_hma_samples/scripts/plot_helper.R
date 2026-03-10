### Function to annotate t_test pvalues in a boxplot (y = value, x = category) comparing two values in the x-axis.
t_test_ann <- function(q, comparison, pval_cutoff = 1.1, y_start = NULL, test = "t_test"){
   x_var = rlang::quo_get_expr(q$mapping$x)
   y_var = rlang::quo_get_expr(q$mapping$y)
   test_df = q$data %>% dplyr::select(x_var, y_var) %>%
      mutate(y = !!y_var, x = !!x_var)
   if (!(is.factor(test_df$x))){
      test_df <- test_df %>%
         mutate(x = factor(x))
   }
   boxplot_df = layer_data(q, 1)
   factor_levels = sapply(comparison, function(x){which(levels(test_df$x) == x)})
   factor_levels_inbetween = min(factor_levels):max(factor_levels)
   x_start = min(factor_levels) - 0.05
   x_end = max(factor_levels) + 0.05
   seg_y = 0.04 * (max(boxplot_df$ymax) - min(boxplot_df$ymin))
   if (is.null(y_start)){
      y_start = max(boxplot_df[factor_levels_inbetween, "upper"]) + seg_y/2
   }
   y_end = y_start + seg_y
   
   test_df = test_df %>% filter(!!x_var %in% comparison)
   t_test = t.test(y ~ x, data = test_df,
                   alternative = "two.sided", var.equal = TRUE)
   w_test = wilcox.test(y ~ x, data = test_df,
                        alternative = "two.sided", var.equal = TRUE)
   if (test == "t_test"){
      stat_test = t_test
   } else if (test == "wilcox_test"){
      stat_test = w_test
   }
   if (stat_test[[3]] < pval_cutoff){
      q = q + annotate("segment", x = x_start, xend = x_end, y = y_end, yend = y_end, linewidth = 0.8) +
         annotate("segment", x = x_start, xend = x_start, y = y_start, yend = y_end, linewidth = 0.8) + 
         annotate("segment", x = x_end, xend = x_end, y = y_start, yend = y_end, linewidth = 0.8) +
         annotate("text", x = (x_start + x_end)/2, y = y_end + seg_y, label = formatC(stat_test[[3]], digits = 2),
                  size = 4.3)
      
   }
   return(q)
}

# q1 = t_test_ann(q, c("Refractory", "Response_no_relapse"))
# q2 = t_test_ann(q1, c("Relapse", "Response_no_relapse"))
# q2

## Helper function to compute t_test between the combined group 'Repsonse_no_relapse + Relapse' vs 'Refractory'
stat_test_ <- function(plot_df){
   test_df = plot_df %>%
      mutate(x_ = case_when(x %in% c("Response_no_relapse", "Relapse") ~ "combined_response_group",
                            TRUE ~ x)) %>%
      filter(x_ != "Paired_relapse_sample")
   t_test = t.test(y ~ x_, data = test_df,
                   alternative = "two.sided", var.equal = TRUE)
   w_test = wilcox.test(y ~ x_, data = test_df,
                        alternative = "two.sided", var.equal = TRUE)
   message_ = paste("t_test_pval =", formatC(t_test[[3]], digits = 2), "\n") %>%
      paste("wilcox_test_pval =", formatC(w_test[[3]], digits = 2))
   print(message_)
}

## Plot CD14 and CD34
## Example of how to use the t_test_ann function
# plot_df <- exprs(m_exp28[c("CD14", "CD34"), ]) %>% as.data.frame() %>% tibble::rownames_to_column("feature") %>%
#    tidyr::pivot_longer(-feature, names_to = "sample_name", values_to = "abundance") %>% 
#    inner_join(p_data, by = "sample_name")
# 
# q = ggplot(plot_df %>% filter(feature == "CD14"), 
#            aes(x = subcohort, y = abundance, fill = group)) + geom_boxplot(width = 0.75) + 
#    ylab("Ven signal") + scale_fill_manual(values = subcohort_colors) +
#    ggtitle("CD14") +
#    theme(text = element_text(size = 15))
# q
# 
# q = ggplot(plot_df %>% filter(feature == "CD34"), 
#            aes(x = subcohort, y = abundance, fill = group)) + geom_boxplot(width = 0.75) + 
#    ylab("Ven signal") + scale_fill_manual(values = subcohort_colors) +
#    ggtitle("CD34") +
#    theme(text = element_text(size = 15))
# q
# 
# stat_test_(plot_df %>% mutate(x = subcohort, y = abundance))
# ## Example of how to use the t_test_ann function
# q1 = t_test_ann(q, c("Refractory", "Response_no_relapse"))
# q2 = t_test_ann(q1, c("Relapse", "Response_no_relapse"), y_start = 3.5)
# q2 = t_test_ann(q2, c("Refractory", "Relapse"), y_start = 4)
# q2
# q1 = t_test_ann(q, c("Refractory", "Response_no_relapse"), test = "wilcox_test")
# q2 = t_test_ann(q1, c("Relapse", "Response_no_relapse"), test = "wilcox_test")
# q2 = t_test_ann(q2, c("Refractory", "Relapse"), test = "wilcox_test", y_start = -0.12)
# q2