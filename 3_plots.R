

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
prim_data <- 

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

# Plot virus-positive map
ggplot() +
  geom_sf(data = world, fill = "grey90", color = "gray40", size = 0.2) +
  geom_tile(data = overlap_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis(name = "NHP overlap\n(virus +)", option = "H") +
  coord_sf(xlim = c(-150, 150), ylim = c(-40, 40), expand = FALSE) +
  ggthemes::theme_map(base_size = 14) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

# ── 2. Create plots using NHP area of habitat rasters ─────────

thresholds <- c(0.1, 0.2, 0.3, 0.5)

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
  
  print(p)
  dev.off()
}

# ── 2. Create partial dependence plots ────────────────────────
