#' Plot stratified ancestry sets.
#'
#' Create a bar plot comparing stratified ancestry sets.
#'
#' @param MX A data frame of train meta.
#' @param MY A data frame of test meta.
#' @param MR A data frame of inference meta.
#' @param x_var Name of the variable for the x-axis (string).
#' @param fill_var Name of the variable for the fill color (string).
#' @param title Plot title (optional).
#' @param x_label Label for the x-axis (optional).
#' @param y_label Label for the y-axis (optional).
#'
#' @return A ggplot object.
#' @export
plot_stratified_sets <- function(
  MX, 
  MY, 
  MR, 
  x_var,
  fill_var,
  title = NULL,
  x_label = NULL,
  y_label = NULL
) {

  ## --- Input data structure check ---
  assert_input(
    MX = MX, 
    MY = MY,
    MR = MR,
    g_col = fill_var, 
    a_col = x_var
  )

  ## --- Prepare x_var label ---
  MR[[x_var]] <- paste0(MR[[x_var]], "\n(Reference, R)")
  MX[[x_var]] <- paste0(MX[[x_var]], "\n(Subset, X)")
  MY[[x_var]] <- paste0(MY[[x_var]], "\n(Inference, Y)")

  # Ensure factor levels
  M <- rbind(MX, MY, MR)
  M[[x_var]] <- factor(
    M[[x_var]], 
    levels = c(
      unique(MR[[x_var]]),  
      unique(MX[[x_var]]),  
      unique(MY[[x_var]]) 
    )
  )

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
    fill = fill_var
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}

