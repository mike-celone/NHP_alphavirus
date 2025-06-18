# Michael Celone
# Date: June 2025
# Project: Predicting undiscovered non-human primate hosts of Semliki Forest complex Alphaviruses
#
# Code for creating several climate and biogeographical variables
# And importing/merging/cleaning the data for BRT analysis

# Data sources include: 
#    (1) COMBINE: a coalesced mammal database of intrinsic and extrinsic traits
#        (https://doi.org/10.1002/ecy.3344)
#
#    (2) Ecological Traits of the World's Primates
#        (https://doi.org/10.1038/s41597-019-0059-9)
#
#    (3) Raster files from several sources, including WorldClim, IUCN, and SEDAC               
##
library(dplyr)
library(tidyr)
library(tidyverse)
library(sp)
library(sf)
library(rgdal)
library(raster)
library(terra)
library(exactextractr)
library(RISmed)
options(scipen = 999) 

#####################################################
# Step 1: Import COMBINE and Ecological Traits data #
#####################################################
setwd("C:/Users/Mike/OneDrive/MAYV_NHP_prediction/data/")

# COMBINE trait data
combine <- read.csv("COMBINE_trait_data_imputed.csv") %>%
  
  # Remove the extinct species
  filter(iucn2020_binomial != "Not recognised") %>%
  
  # Remove a few columns
  dplyr::select(-c(order, genus, species)) %>%
  
  # Fixing up some of the species vs sub-species discrepencies
  mutate(iucn2020_binomial = case_when(
    iucn2020_binomial == "Lagothrix lagothricha" & phylacine_binomial == "Lagothrix cana" ~ "Lagothrix cana",
    iucn2020_binomial == "Lagothrix lagothricha" & phylacine_binomial == "Lagothrix lugens" ~ "Lagothrix lugens",
    iucn2020_binomial == "Lagothrix lagothricha" & phylacine_binomial == "Lagothrix poeppigii" ~ "Lagothrix poeppigii",
    iucn2020_binomial == "Sapajus apella" & phylacine_binomial == "Sapajus macrocephalus" ~ "Sapajus macrocephalus",
    TRUE ~ iucn2020_binomial))

# Bring in Ecological Traits of the World's Primates data
mass <- read.csv("ETWP_BodyMass.csv") 
diel <- read.csv("ETWP_DielActivity.csv")
habitat <- read.csv("ETWP_Habitat.csv")
iucn <- read.csv("ETWP_IUCN_Poptrend_Realm.csv")
locomotion <- read.csv("ETWP_Locomotion.csv")
trophic <- read.csv("ETWP_TrophicGuild.csv")
etwp <- list(diel, mass, habitat, iucn, locomotion, trophic)

# Fixing up a few issues with ETWP species
etwp_new <- list()
for (i in 1:length(etwp)){
  etwp_new[[i]] <- etwp[[i]] %>%
    mutate(
      
      #Remove the underscore
      Species = str_replace(Species, "_", " "),
      Species = case_when(
        
        #If binomial is NA fill in with the ITIS binomial
        is.na(Species) ~ Species..ITIS.,
        
        #Fix up a few binomials to match the IUCN spelling
        Species == "Galagoides demidovii" ~ "Galagoides demidoff",
        Species == "Lepilemur ahmansonorum" ~ "Lepilemur ahmansonori",
        Species == "Lepilemur hubbardi" ~ "Lepilemur hubbardorum",
        Species == "Lepilemur sahamalazensis" ~ "Lepilemur sahamalaza",
        Species == "Piliocolobus waldronae" ~ "Piliocolobus waldroni",
        TRUE ~ as.character(Species)))
}
# Note that some species have conflicting entries for Trophic Guild and Locomotion and these were changed to NA
# When there are multiple entries for mass per species I took the average

### One-hot encoding/cleaning some of the variables

# IUCN variables
iucn_new <- data.frame(etwp_new[[4]]) %>%
  mutate(valuea = 1, valueb = 1) %>%
  spread(IUCN, valuea, fill = 0) %>%
  spread(Pop_T, valueb, fill = 0) %>%
  dplyr::select(Family, Species, 10:19) %>%
  # renaming some variables for easier processing
  rename(threat_critical=CR, threat_endangered=EN, threat_vulnerable=VU, 
         threat_near_threatened=NT, threat_least_concern=LC, threat_data_deficient=DD,
         threat_not_evaluated=NE, ppoulation_increasing=I, population_stable=S,
         population_decreasing=D)

# ETWP: Habitat
habitat_new <- data.frame(etwp_new[[3]]) %>%
  dplyr::select(Family, Species, 6:12)

# ETWP: Trophic
trophic <- data.frame(etwp_new[[6]]) %>%
  dplyr::select(Family, Species, TrophicGuild) %>%
  distinct()  

# Remove species with conflicting data -- setting these to NA
trophic$dup <- duplicated(trophic$Species)
dups <- trophic %>% filter(dup==T) %>% dplyr::select(Species)
trophic <- trophic[ ! trophic$Species %in% dups$Species, ]
trophic_new <- trophic %>%
  mutate(value = 1) %>%
  spread(TrophicGuild, value, fill = 0) %>%
  dplyr::select(1:2, 4:9)

# ETWP: Locomotion
locomotion <- data.frame(etwp_new[[5]]) %>%
  dplyr::select(Family, Species, Locomotion) %>%
  distinct()

# Remove duplicate species with conflicting data -- setting these to NA
locomotion$dup <- duplicated(locomotion$Species)
dups_loc <- locomotion %>% filter(dup==T) %>% dplyr::select(Species)
locomotion <- locomotion[ ! locomotion$Species %in% dups_loc$Species, ]

locomotion_new <- locomotion %>%
  mutate(value = 1) %>%
  spread(Locomotion, value, fill = 0) %>%
  dplyr::select(Family, Species, locomotion_arboreal=AR, locomotion_both=BOTH, locomotion_terrestrial=`T`) 

# ETWP: Mass
mass_new <-data.frame(etwp_new[[2]]) %>%
  group_by(Species) %>%
  dplyr::select(Family, Species, BodyMassMale_kg, BodyMassFemale_kg) %>% 
  summarize(massMale_kg = mean(BodyMassMale_kg, na.rm =TRUE, na.action = na.pass),
            massFemale_kg = mean(BodyMassFemale_kg, na.rm = TRUE, na.action = na.pass))

# Join the ETWP datasets together
join_etwp <- mass_new %>%
  full_join(habitat_new) %>% 
  full_join(iucn_new) %>% 
  full_join(locomotion_new) %>% 
  full_join(trophic_new)  

# Merge ETWP and COMBINE
combined_traits <- right_join(join_etwp, combine, by=c("Species"="iucn2020_binomial", "Family"="family")) %>%
  mutate(biogeographical_realm = case_when(
    
    # Create a new category for "multiple realms" 
    grepl(",", biogeographical_realm) ~ "realm_multiple",
    TRUE ~ as.character(biogeographical_realm)),
    
    # Recode activity cycle
    activity_cycle = case_when(
      activity_cycle == 1 ~ "activity_cycle_nocturnal",
      activity_cycle == 3 ~ "activity_cycle_diurnal",
      activity_cycle == 2 ~ "activity_cycle_multiple"),
    
    # One-hot encoding
    valuea = 1, valueb = 1, valuec = 1, valued = 1)  %>% 
  spread(Family, valuea,  fill = 0) %>%
  spread(biogeographical_realm, valueb, fill = 0) %>%
  spread(foraging_stratum, valuec, fill = 0) %>%
  spread(activity_cycle, valued, fill = 0) %>%
  
  # Rename some variables
  rename(foraging_stratum_arboreal=Ar, 
         foraging_stratum_ground=G, 
         foraging_stratum_scansorial=S) %>%
  dplyr::select(-`<NA>`)

# Calculating derived variables from COMBINE data
combined_traits <- combined_traits %>%
  mutate(
    body_size_ratio = adult_mass_g/neonate_mass_g,
    rel_age_fb_days = age_first_reproduction_d/max_longevity_d,
    rel_sexual_maturity_days = maturity_d/max_longevity_d,
    postnatal_growth_rate = weaning_mass_g/neonate_mass_g,
    # Formula for production from Hamilton et al., 2010
    mass_specific_production = (neonate_mass_g/adult_mass_g)*litter_size_n*litters_per_year_n,
    home_range_scaled_km2_g = home_range_km2/adult_mass_g)

names(combined_traits) <- tolower(names(combined_traits))

# Add the WoS hits
wos <- read.csv("wos_hits.csv")
combined_traits <- left_join(combined_traits, wos, by = "species")

#####################################################
# Step 2: Extract variables from species range maps #
#####################################################

# Set working directory
setwd("C:/Users/Mike/OneDrive/MAYV_NHP_prediction/Rasters/")

# List the species range rasters
files <- list.files("Mammals_primates/")

# Download the WorldClim rasters if necessary
#worldclim_rasters <- getData("worldclim",var="bio",res=10)
raster_list <- list.files(path = paste0(getwd(), "/wc10/"), pattern = ".bil", full.names = T)
worldclim_rasters <- stack(raster_list)
bio1 <- rast(worldclim_rasters[["bio1"]])
bio5 <- rast(worldclim_rasters[["bio5"]])
bio6 <- rast(worldclim_rasters[["bio6"]])
bio12 <- rast(worldclim_rasters[["bio12"]])
bio13 <- rast(worldclim_rasters[["bio13"]])
bio14 <- rast(worldclim_rasters[["bio14"]])

# Import additional rasters
popdens <- rast("gpw_v4_population_density_rev11_2020_15_min.asc")
species_rich <- rast("species_richness.tif")
hfp <- rast("hfp_resampled.tif")

# Bring in WWF ecoregions
ecor <- st_read("tnc_terr_ecoregions.shp")
ecor_rast <- rasterize(ecor, bio1, field = "ECO_CODE") # Rasterize the polygons

# Choose a reference raster — for example, bio1
ref_rast <- bio1

# Resample others to match the reference (do not need to resample the other WorldClim rasters)
popdens_res <- resample(popdens, ref_rast)
species_rich_res <- resample(species_rich, ref_rast)
hfp_res <- resample(hfp, ref_rast)
ecor_res <- resample(ecor_rast, ref_rast)

# Stack them
env_stack <- c(ref_rast, bio5, bio6, bio12, bio13, bio14,
               popdens_res, species_rich_res, hfp_res, ecor_res)

names(env_stack) <- c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14", "popdens", "rich", "hfp", "ecor")

# Loop through each species raster and extract the relevant information

# Create output directory if it doesn't exist
dir.create("results_tmp_16Jun", showWarnings = FALSE)

# Placeholder for results
results_list <- list()

# Loop through all files
for (file_path in files) {
  message("▶ Processing: ", file_path)
  
  result_df <- tryCatch({
    # Load raster
    r <- rast(file.path("Mammals_primates", file_path))
    
    # Crop and mask
    crop_ras <- crop(env_stack, r)
    r_resampled <- resample(r, crop_ras, method = "near")
    masked <- mask(crop_ras, r_resampled, maskvalue = NA)
    
    # Latitude stats
    non_na_cells <- which(!is.na(values(r)))
    coords <- xyFromCell(r, non_na_cells)
    latitudes <- coords[, 2]
    latitude <- data.frame(ymin = min(latitudes), ymax = max(latitudes)) %>%
      mutate(range_lat = ymax - ymin, centroid = (ymin + ymax) / 2)
    
    # Climate values
    mean_temp  <- global(masked[["bio1"]], "mean", na.rm = TRUE)[1,1]
    max_temp   <- global(masked[["bio5"]], "max",  na.rm = TRUE)[1,1]
    min_temp   <- global(masked[["bio6"]], "min",  na.rm = TRUE)[1,1]
    range_temp <- global(masked[["bio5"]] - masked[["bio6"]], "mean", na.rm = TRUE)[1,1]
    
    mean_prec  <- global(masked[["bio12"]], "mean", na.rm = TRUE)[1,1]
    max_prec   <- global(masked[["bio13"]], "max",  na.rm = TRUE)[1,1]
    min_prec   <- global(masked[["bio14"]], "min",  na.rm = TRUE)[1,1]
    range_prec <- global(masked[["bio13"]] - masked[["bio14"]], "mean", na.rm = TRUE)[1,1]
    
    # Population and richness
    mean_popdens <- global(masked[["popdens"]], "mean", na.rm = TRUE)[1,1]
    max_popdens  <- global(masked[["popdens"]], "max",  na.rm = TRUE)[1,1]
    min_popdens  <- global(masked[["popdens"]], "min",  na.rm = TRUE)[1,1]
    species_avg  <- global(masked[["rich"]], "mean",   na.rm = TRUE)[1,1]
    hfp_avg      <- global(masked[["hfp"]], "mean",    na.rm = TRUE)[1,1]
    
    # Unique ecoregions
    unique_ecos <- unique(values(masked[["ecor"]]))
    unique_ecos <- unique_ecos[!is.na(unique_ecos)]
    ecor_n <- length(unique_ecos)
    
    # Assemble result
    binomial <- gsub("_", " ", sub("\\.tif$", "", file_path))
    
    data.frame(
      binomial = binomial,
      ecor_n = ecor_n,
      hfp_avg = hfp_avg,
      species_avg = species_avg,
      max_popdens = max_popdens,
      min_popdens = min_popdens,
      mean_popdens = mean_popdens,
      mean_temp = mean_temp,
      max_temp = max_temp,
      min_temp = min_temp,
      range_temp = range_temp,
      mean_prec = mean_prec,
      max_prec = max_prec,
      min_prec = min_prec,
      range_prec = range_prec,
      ymin = latitude$ymin,
      ymax = latitude$ymax,
      range_lat = latitude$range_lat,
      centroid = latitude$centroid
    )
    
  }, error = function(e) {
    message("❌ Skipped ", file_path, ": ", e$message)
    NULL
  })
  
  if (!is.null(result_df)) {
    out_file <- file.path("results_tmp", paste0(tools::file_path_sans_ext(file_path), ".csv"))
    write_csv(result_df, out_file)
    message("✅ Saved: ", out_file)
  }
}

# List all batch files
batch_files <- list.files("results_tmp", pattern = ".csv$", full.names = TRUE)

# Combine them
all_results <- bind_rows(lapply(batch_files, read_csv))

# Save final full result
write_csv(all_results, "species_environment_summary.csv")

# Join the extracted variables with the trait dataset 
combined_traits <- combined_traits %>%
  distinct(species, .keep_all = TRUE)

final_df <- right_join(combined_traits, all_results, by=c("species"="binomial")) %>%
  dplyr::select(-phylacine_binomial)

##############################################
# Step 3: Final steps to clean up trait data #
##############################################

#Remove variables with >80% NA (i.e., 300 or more NA)
na_rem <- nrow(final_df)*.8
prim_data_reduced <- final_df[, colSums(is.na(final_df)) < na_rem]

#Check the variance of the variables
no_var <- caret::nearZeroVar(
  prim_data_reduced[2:ncol(prim_data_reduced)],
  freqCut = 90/10,
  uniqueCut = 10,
  saveMetrics = T,
  names = F)

# Get variable names (i.e. rownames) where zeroVar is TRUE
zero_var_cols <- rownames(no_var)[no_var$zeroVar == TRUE]

#Remove the variables with little variance
prim_data_final <- prim_data_reduced %>%
  dplyr::select(-c(zero_var_cols,
                   
                   # Irrelevant variables
                   hibernation_torpor, 
                   fossoriality,
                   
                   # Redundant
                   dphy_invertebrate,
                   dphy_vertebrate,
                   dphy_plant,
                   trophic_level,
                   island_endemicity))

# Replace some NaN with NA
prim_data_final[prim_data_final == "NaN"] <- NA

# Re-order a few cols
prim_data_final <- prim_data_final %>%
  dplyr::select(species, vir_pos, everything())

# Variable coverage
missing_vars <- prim_data_final %>%
  dplyr::select(everything()) %>%  
  summarise_all(
    funs(sum(!is.na(.))/nrow(prim_data_final)))

missing_species <- rowSums(is.na(prim_data_final))

# Checking skewness
skew <- sapply(c(3:4, 31:50, 60:62, 94:118), function(x) moments::skewness(na.omit(prim_data_final[, x])))
skew <- colnames(prim_data_final)[c(3:4, 31:50, 60:62, 94:118)][which(skew > 2)]

# Log transform the skewed vars
for (i in skew) {
  col_vals <- prim_data_final[[i]]
  
  # Check if any values are ≤ 0
  if (any(col_vals <= 0, na.rm = TRUE)) {
    # Shift values: add 1 to avoid log(0), or shift entire vector if needed
    prim_data_final[[i]] <- log(col_vals + 1)
  } else {
    # Safe to log directly
    prim_data_final[[i]] <- log(col_vals)
  }
}
sapply(prim_data_final[skew], function(x) sum(is.infinite(x) | is.nan(x)))
write.csv(prim_data_final, "prim_data_final.csv")
