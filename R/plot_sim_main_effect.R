#' Scatter Plot of Realized log2 Fold-Change X vs Y
#'
#' Create a scatter plot comparing realized or observed log2 fold-change
#' effects for two ancestries (populations). This helps visualize how the 
#' main effects in X and Y relate to each other, highlighting features 
#' that are DEGs in one or both populations.
#'
#' Points are colored by DEG status (Both, X-only, Y-only, None) and shaped 
#' by presence of an interaction effect. Intended effect sizes and a 
#' biological threshold band can be shown for reference.
#'
#' @param fX Data frame for ancestry X features (must contain `true_log2FC` or `observed_log2FC`).
#' @param fY Data frame for ancestry Y features.
#' @param x_var Column name in `fX` to use for X-axis effect (e.g., "true_log2FC").
#' @param y_var Column name in `fY` to use for Y-axis effect (e.g., "true_log2FC").
#' @param intended_main_effect_X Optional numeric intended log2FC for X ancestry.
#' @param intended_main_effect_Y Optional numeric intended log2FC for Y ancestry.
#' @param effect_thr Optional numeric threshold for interaction effect band (e.g., 1).
#' @param features Optional character vector of feature names to label.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param caption Plot caption (optional).
#' @param point_size Numeric size of scatter points (default: 1).
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline geom_line geom_ribbon
#' @importFrom ggplot2 scale_color_manual scale_shape_manual scale_fill_manual scale_linetype_manual
#' @importFrom ggplot2 labs theme
#' @importFrom ggrepel geom_text_repel
plot_sim_main_effect <- function(
  fX,
  fY,
  x_var,
  y_var,
  intended_main_effect_X = NULL,
  intended_main_effect_Y = NULL,
  effect_thr = NULL,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  caption = NULL,
  point_size = 1
) {

  # Unique intended FCs for X and Y
  X_intended <- intended_main_effect_X
  Y_intended <- intended_main_effect_Y

  X_a <- unique(fX$ancestry)
  Y_a <- unique(fY$ancestry)

  # Combined effect data frame
  combined_df <- data.frame(
    X_effect = fX[[x_var]],
    Y_effect = fY[[y_var]],
    X_true_effect = fX$true_log2FC,
    Y_true_effect = fY$true_log2FC,
    is_DE_X = as.logical(fX$is_DE),
    is_DE_Y = as.logical(fY$is_DE),
    feature = fX$feature,
    stringsAsFactors = FALSE
  )

  # Annotate features if provided
  combined_df$annotate <- if (!is.null(features)) {
    ifelse(
      combined_df$feature %in% features, 
      combined_df$feature, 
      NA
    )
  } else {
    NA
  }

  # Fill by true DEGs
  X_status <- paste(X_a, "only")
  Y_status <- paste(Y_a, "only")

  combined_df$DEG_status <- with(
    combined_df, ifelse(
      is_DE_X & is_DE_Y, "Both",
      ifelse(is_DE_X, X_status,
        ifelse(is_DE_Y, Y_status, "None")
      )
    )
  )

  combined_df$DEG_status <- factor(
    combined_df$DEG_status,
    levels = c("Both", X_status, Y_status, "None")
  )

  # Order of point layers
  plot_order <- ifelse(
    combined_df$DEG_status == "Both", 3,
    ifelse(combined_df$DEG_status == "None", 1, 2)
  )
  combined_df <- combined_df[order(plot_order), ]

  # Shape by interactions 
  true_interaction <- combined_df$Y_true_effect - combined_df$X_true_effect
  is_DE <- as.numeric(abs(true_interaction) > 0)
  combined_df$Interaction_status <- factor(
    is_DE,
    levels = c(0, 1),
    labels = c("No interaction", "Interaction")
  )

  # Lines for intended only
  vline_type <- paste(X_a, "injected effect")
  vline_df <- data.frame(x = c(0, 0), y = c(0, 0), type = vline_type)

  hline_type <- paste(Y_a, "injected effect")
  hline_df <- data.frame(x = c(0, 0), y = c(0, 0), type = hline_type)

  # Abline for 0 interactions
  slope <- 1
  intercept <- 0
  # Span the full x range you want
  x_min <- min(combined_df$X_effect, na.rm = TRUE)
  x_max <- max(combined_df$X_effect, na.rm = TRUE)
  x_seq <- seq(x_min, x_max, length.out = 200)
  y_seq <- slope * x_seq + intercept

  line_df <- data.frame(
    x = x_seq, 
    y = y_seq
  )

  # Ribbon for threshold
  if (!is.null(effect_thr)) {
    line_df$ymin <- y_seq - effect_thr
    line_df$ymax <- y_seq + effect_thr
    line_df$Band <- factor(paste("±", effect_thr, "log2FC"))
  }

  # Base plot
  p <- ggplot(
    data = combined_df,
    mapping =  aes(
      x = X_effect, 
      y = Y_effect
      )
    ) +
    geom_point(
      mapping = aes(
        color = DEG_status,
        shape = Interaction_status
      ),
      size = point_size
    ) +
    scale_color_manual(
      values = setNames(
        c(
          "purple", 
          "#1f77b4", 
          "#ff7f0e", 
          "gray80"
        ),
        c(
          "Both", 
          X_status, 
          Y_status, 
          "None"
        )
      ),
      breaks = c(
        "Both", 
        X_status, 
        Y_status, 
        "None"
      )
    ) +
    scale_shape_manual(
      values = c(
        "No interaction" = 1, 
        "Interaction" = 17
      ),
      drop = FALSE
    )


    # Ribbon if effect thr is defined (biological interpretation)
    if (!is.null(effect_thr)) {
      p <- p + geom_line(
        data = line_df,
        mapping = aes(
          x = x, 
          y = y,
        ),
        color = "grey50",
        linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_ribbon(
        data = line_df,
        mapping = aes(
          x = x, 
          ymin = ymin, 
          ymax = ymax,
          fill = Band,
          group = Band
        ),
        alpha = 0.2,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(
        values = setNames(
          "grey50", 
          unique(line_df$Band)
        )
      )
    } 

    # Add intended effect size X
    if (!is.null(intended_main_effect_X)){
      p <- p + geom_vline(
        xintercept = c(
          -X_intended,
          X_intended
        ),
        linetype = "dashed",
        linewidth = 0.3,
        color = "red",
        show.legend = FALSE
      ) +
      geom_line(
        data = vline_df,
        aes(
          x = x, 
          y = y, 
          linetype = type
        ),
        color = "red"
      )
    }

    if (!is.null(intended_main_effect_Y)){
      p <- p + geom_hline(
        yintercept = c(
          -Y_intended,
          Y_intended
        ),
        linetype = "dashed",
        linewidth = 0.3,
        color = "red",
        show.legend = FALSE
      ) +
      geom_line(
        data = hline_df,
        aes(
          x = x, 
          y = y, 
          linetype = type
        ),
        color = "red"
      )
    }

    # Only add scale_linetype_manual if either was set
    linetype_vals <- c()
    linetype_names <- c()

    if (!is.null(intended_main_effect_X)) {
      linetype_vals <- c(
        linetype_vals, 
        "dashed"
      )
      linetype_names <- c(
        linetype_names, 
        vline_type
      )
    }
    if (!is.null(intended_main_effect_Y)) {
      linetype_vals <- c(
        linetype_vals, 
        "dashed"
      )
      linetype_names <- c(
        linetype_names, 
        hline_type
      )
    }
    if (length(linetype_vals) > 0) {
      p <- p + scale_linetype_manual(
        values = setNames(
          linetype_vals, 
          linetype_names
        )
      )
    }

    # Final styling of the plot
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), paste(X_a, x_var), x_label),
      y = ifelse(is.null(y_label), paste(Y_a, y_var), y_label),
      caption = caption,
      linetype = "Intended effect size",
      fill = "Interaction threshold",
      color = "Main effect",
      shape = "Interaction effect"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  # Add labels if features are provided
  if (!is.null(features)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(combined_df, !is.na(annotate)),
      aes(
        x = X_effect,
        y = Y_effect,
        label = annotate
      ),
      size = 2,
      min.segment.length = 0,
      segment.size = 0.3,
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5,
      bg.color = "white",
      bg.r = 0.2,
      inherit.aes = FALSE
    )
  }

  return(p)
}
