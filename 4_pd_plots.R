library(pdp)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)

# Set file paths
root <- "C:/Users/Mike/OneDrive/MAYV_NHP_prediction/"
model_run <- Sys.Date()
run_folder <- paste0(root, "/Analysis/NHP_model_run_", model_run)

# Import variables from previous script
red_aggr <- read.csv(paste0(run_folder, "/var_importance_reduced.csv"))

# Select top 12 variables by mean relative importance
top_vars <- red_aggr |>
  arrange(desc(mean_rel.inf)) |>
  slice_head(n = 12) 

# Cleaning up names and adding category for plotting later
names <- c("Latitude Range", "Max. Longevity (days)", "Max. Temperature (C)",
           "Min. Pop. Density", "Ecoregions (n)", "Neonate Mass (g)",  "Mass Female (kg)",   
           "Mean Precipitation (in)", "Home Range (km2)", "Social Group Size (n)",  "Human Footprint",        
           "Weaning Age (days)")
type <- c("biogeo", "biogeo", "biogeo", "human", "biogeo", "life", 
          "life", "biogeo", "life", "life", "human", "life")
top_vars <- data.frame(cbind(top_vars, names = names, type = type))

# Import model object from previous script
models <- final_reduced$models

# Function to compute PDPs for one model and one variable
get_pdp <- function(model_obj, variable, train_df) {
  pdp::partial(
    object   = model_obj$fit,
    pred.var = variable,
    train    = train_df,
    prob     = TRUE,
    n.trees  = model_obj$best_iter,
    type     = "classification",  # <-- explicitly set this
    plot     = FALSE
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(variable = variable)
}

# Loop over variables and models to get PDPs
all_pdps <- purrr::map_dfr(top_vars, function(v) {
  purrr::map_dfr(models, ~get_pdp(.x, v, train_df = prim_data_selected))
})

# Summarize across bootstraps: mean + 95% CI
# Split the data frame by variable
pdp_list <- split(all_pdps, all_pdps$variable)

# For each variable, summarize by the values of that variable's column
pdp_summary_list <- map2(pdp_list, names(pdp_list), function(df, var_name) {
  df %>%
    group_by(val = .data[[var_name]]) %>%
    summarise(
      mean_yhat = mean(yhat, na.rm = TRUE),
      lower     = quantile(yhat, 0.025, na.rm = TRUE),
      upper     = quantile(yhat, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = var_name)
})

# Combine all back into one data frame
pdp_summary <- bind_rows(pdp_summary_list)

# Create list of plots with error handling
# Create list of plots with error handling
pdp_plots <- lapply(top_vars$var, function(v) {
  tryCatch({
    pdp_subset <- filter(pdp_summary, variable == v)
    
    hist_data <- prim_data_selected[[v]]
    hist_obj <- hist(hist_data, plot = FALSE)
    hist_df <- data.frame(
      mids = hist_obj$mids,
      counts = hist_obj$counts
    )
    
    # Fixed PDP y-axis range: 0 to 0.33
    pdp_min <- 0
    pdp_max <- 0.33
    hist_max <- max(hist_df$counts, na.rm = TRUE)
    scale_factor <- hist_max / (pdp_max - pdp_min)
    
    pdp_scaled <- pdp_subset %>%
      mutate(
        scaled_mean = mean_yhat * scale_factor,
        scaled_lower = lower * scale_factor,
        scaled_upper = upper * scale_factor
      )
    
    # Set color based on variable type
    type <- top_vars$type[top_vars$var == v]
    color <- if (type == "human") "#F8766D"
    else if (type == "biogeo") "#619CFF"
    else "#00BA38"
    name <- top_vars$name[top_vars$var == v]
      
    ggplot() +
      geom_col(data = hist_df,
               aes(x = mids, y = counts),
               fill = "gray80", color = "black", alpha = 0.5,
               width = diff(hist_obj$breaks)[1]) +
      geom_ribbon(data = pdp_scaled,
                  aes(x = val, ymin = scaled_lower, ymax = scaled_upper),
                  fill = color, alpha = 0.3) +
      geom_line(data = pdp_scaled,
                aes(x = val, y = scaled_mean),
                color = color, size = 1.2) +
      scale_y_continuous(
        name = NULL,
        sec.axis = sec_axis(~ . / scale_factor, name = NULL, breaks = seq(0, 0.5, 0.1))
      ) +
      labs(x = name, title = NULL) +
      theme_bw() +
      theme(
        #axis.title = element_blank(),
        axis.text = element_text(size = 8),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "none"
      )
  }, error = function(e) {
    message(paste("Skipping", v, "due to error:", e$message))
    NULL
  })
})

# Remove NULLs
pdp_plots <- Filter(Negate(is.null), pdp_plots)

# Plot in grid
library(patchwork)
grid_plot <- wrap_plots(pdp_plots, ncol = 4)
ggsave(filename = paste0("PDP.png"), plot = grid_plot, width = 10, height = 6)

