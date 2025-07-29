#' Plot feature-level summary using sim_imbalanced_groups summary output
#'
#' @param summary_list The `summary` list from `sim_imbalanced_groups`.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_sim_interaction_effect_demographics <- function(
    summary_list, 
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    fill_name = NULL
) {

  # Information on feature summary
  feature_summary <- summary_list$interaction[, c("ancestry", "total_DEGs", "total_non_DEGs")]
  
  # Reshape
  long_df <- reshape(
    feature_summary,
    direction = "long",
    varying = c("total_DEGs", "total_non_DEGs"),
    v.names = "count",
    times = c("DEGs", "Non-DEGs"),
    timevar = "status",
    idvar = "ancestry"
  )
  
  # Calculate percentage within ancestry
  total_counts <- tapply(long_df$count, long_df$ancestry, sum)
  long_df$pct <- round(100 * long_df$count / total_counts[long_df$ancestry], 1)
  
  # Convert Status to factor for stacking order
  long_df$status <- factor(long_df$status, levels = c("DEGs", "Non-DEGs"))
  
  # Plot
  p <- ggplot(
    data = long_df, 
    mapping = aes(
        x = ancestry, 
        y = count, 
        fill = status
      )
    ) +
    geom_bar(
      stat = "identity", 
      position = "dodge", 
      color = "black"
    ) +
    geom_text(
      mapping = aes(
        label = paste0(count, " (", pct, "%)")
      ),
      position = position_dodge(
        width = 1,
      ),
      vjust = -0.5,
      size = 2,
      color = "black"
    ) +
    scale_fill_manual(
      values = setNames(
        c(
          "gold", 
          "gray80"
          ),
        c(
          "DEGs", 
          "Non-DEGs"
        )
      ),
      breaks = c(
        "DEGs",
        "Non-DEGs"
      )
    ) +
    labs(
      title = title,
      x = x_label,
      y = ifelse(is.null(y_label), "Count", y_label),
      fill = fill_name
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()
  
  return(p)
}
