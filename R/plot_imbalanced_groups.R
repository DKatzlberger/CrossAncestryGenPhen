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
#' @param y_label Label for the fill (optional).
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
  y_label = NULL,
  fill_label = NULL
){

  ## --- Input data structure check ---
  assert_input(
    MX = MX, 
    MY = MY,
    g_col = fill_var, 
    a_col = x_var
  )

  ## --- Ancetsry levels ----
  a_1 <- unique(MX[[x_var]]); a_2 <- unique(MY[[x_var]])
  a_levels <- c(a_1, a_2)

  ## --- Bind frames ---
  M <- rbind(MX, MY)
  M[[x_var]] <- factor(M[[x_var]], levels = a_levels)

  ## --- Bar plot ---
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
    y = ifelse(is.null(y_label), "Count", y_label),
    fill = ifelse(is.null(fill_label), fill_var, fill_label),
  ) +
  theme_CrossAncestryGenPhen()

  ## --- Return ----
  return(p)
}