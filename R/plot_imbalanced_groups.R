#' Plot Imbalanced Groups
#'
#' Create a bar plot comparing two groups with counts by category.
#'
#' @param MX A data frame for the first group.
#' @param MY A data frame for the second group.
#' @param x_var Name of the variable for the x-axis (string).
#' @param fill_var Name of the variable for the fill color (string).
#' @param title Plot title (optional).
#' @param x_label Label for the x-axis (optional).
#' @param y_label Label for the y-axis (optional).
#'
#' @return A ggplot object.
#' @export
plot_imbalanced_groups <- function(
  MX,
  MY,
  x_var,
  fill_var,
  title = NULL,
  x_label = NULL,
  y_label = NULL
){

  M <- rbind(MX, MY)

  p <- ggplot(
    data = M,
    mapping = aes(
      x = .data[[x_var]],
      fill = .data[[fill_var]]
    )
  ) +
  geom_bar(
    position = "dodge",
    color = "black",
    linewidth = 0.1
  ) +
  labs(
    title = title,
    x = ifelse(is.null(x_label), x_var, x_label),
    y = ifelse(is.null(y_label), "count", y_label),
    fill = fill_var
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}