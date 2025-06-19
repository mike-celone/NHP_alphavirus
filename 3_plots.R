

# ── 0. Setup ─────────────────────────────────────────────────
library(terra)
library(dplyr)
library(ggplot2)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Set file paths
root <- "C:/Users/Mike/OneDrive/MAYV_NHP_prediction/"
model_run <- Sys.Date()
run_folder <- paste0(root, "/Analysis/NHP_model_run_", model_run)
raster_path <- paste0(root, "/Rasters/Mammals_primates/")

# Use worldclim as a template raster to match extents
template <- rast(paste0(root, "/Rasters/wc10/bio1.bil"))

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Import the dataset with labels
prim_data <- read.csv(paste0(root, "/Analysis/prim_data_final.csv"))

# Select virus-positive species
positive_species <- prim_data %>%
  filter(vir_pos == 1) %>%
  pull(species)

species_files <- paste0(gsub(" ", "_", positive_species), ".tif")

# Import the predicted probability scores
probs <- read.csv(paste0(run_folder, "/predictions_reduced.csv"))

# Average scores across bootstrap runs
probs <- probs %>% group_by(species) %>% summarize(pred = mean(prediction))

# ── 1. Setup ─────────────────────────────────────────────────

# Load and align rasters
load_and_align <- function(species_name) {
  path <- file.path(raster_path, species_name)
  if (!file.exists(path)) return(NULL)
  r <- rast(path)
  r <- crop(r, template)
  resample(r, template, method = "near")
}

rasters <- lapply(species_files, load_and_align)
r_stack <- rast(rasters)
overlap <- sum(r_stack, na.rm = TRUE)
overlap_df <- as.data.frame(overlap, xy = TRUE, na.rm = TRUE)
colnames(overlap_df) <- c("x", "y", "count")
max_value <- max(overlap_df$count)
# Plot virus-positive map
g <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "gray40", size = 0.2) +
  geom_tile(data = overlap_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis(name = "Number of known\nhost species", option = "H", breaks = seq(0, max_value, by = 2)) +
  coord_sf(xlim = c(-150, 150), ylim = c(-40, 40), expand = FALSE) +
  ggthemes::theme_map(base_size = 14) +
  theme(legend.text = element_text(size = 11),
        #legend.position = c(.05, .1)
        )

# Save to PNG
ggsave("known_host_plot.png", g, width = 15, height = 12, dpi = 300)


# ── 2. Create plots using NHP area of habitat rasters ─────────

thresholds <- c(0.1, 0.3, 0.5)
all_plots <- list()

for (thresh in thresholds) {
  cat("Processing threshold:", thresh, "\n")
  
  # Filter species
  selected_species <- probs %>%
    filter(pred >= thresh) %>%
    pull(species)
  
  species_files <- paste0(gsub(" ", "_", selected_species), ".tif")
  
  # Load rasters
  rasters <- lapply(species_files, load_and_align)
  
  if (length(rasters) == 0) next  # skip if none
  
  r_stack <- rast(rasters)
  overlap <- sum(r_stack, na.rm = TRUE)
  overlap_df <- as.data.frame(overlap, xy = TRUE, na.rm = TRUE)
  colnames(overlap_df) <- c("x", "y", "count")
  
  # Plot filename
  filename <- paste0("overlap_threshold_", thresh, ".png")
  
  png(filename, width = 8, height = 6, units = "in", res = 300)
  
  # Plot
  p <- ggplot() +
    geom_sf(data = world, fill = "grey90", color = "gray40", size = 0.2) +
    geom_tile(data = overlap_df, aes(x = x, y = y, fill = count)) +
    scale_fill_viridis(name = paste0("NHP overlap\n(p ≥ ", thresh, ")"),
                       option = "H") +
    coord_sf(xlim = c(-150, 150), ylim = c(-40, 40), expand = FALSE) +
    ggthemes::theme_map(base_size = 14) +
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 10))
  
  all_plots[[as.character(thresh)]] <- p
}

# Combine vertically if any plots exist
if (length(all_plots) > 0) {
  final_plot <- wrap_plots(all_plots, ncol = 1)
  
  # Save to one PNG
  ggsave("combined_overlap_plots.png", final_plot, width = 8, height = 6 * length(all_plots), dpi = 300)
}
