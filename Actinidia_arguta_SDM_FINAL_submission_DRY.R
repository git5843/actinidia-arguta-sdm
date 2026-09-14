# Actinidia arguta climatic suitability and early spring frost risk in Europe.
# One-script workflow: model fitting, projection, figures, tables, rasters.

# ---- 1. SESSION SETUP ----

rm(list = ls())
graphics.off()

set.seed(123)
options(stringsAsFactors = FALSE)
terra::terraOptions(progress = 1)

# ---- 2. PACKAGES ----

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "sf", "terra", "geodata", "usdm",
  "biomod2", "ggplot2", "ggspatial", "ggnewscale", "grid", "scales",
  "viridis", "rnaturalearth", "rnaturalearthdata", "biscale", "cowplot",
  "patchwork", "pROC", "spThin", "maxnet", "writexl", "reshape2"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing. Install them before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ---- 3. USER SETTINGS ----

# Falls back to the project path if occurrence.txt isn't in the working directory.
occurrence_file <- "occurrence.txt"
if (!file.exists(occurrence_file)) {
  occurrence_file <- "C:/Program Files/RStudio/Minikiwiproject/occurrence.txt"
}
if (!file.exists(occurrence_file)) {
  stop(
    "Occurrence file not found. Put occurrence.txt in the working directory or ",
    "update 'occurrence_file' in Section 3."
  )
}

output_root <- "FINAL_EXPORTS"

dir.create(output_root, showWarnings = FALSE)
dir.create(file.path(output_root, "MAIN"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "SUPPLEMENTARY"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "TABLES"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "RASTERS"), recursive = TRUE, showWarnings = FALSE)

EU_XLIM <- c(-10, 40)
EU_YLIM <- c(34, 70)
ASIA_XLIM <- c(100, 160)
ASIA_YLIM <- c(20, 60)

extent_asia <- terra::ext(ASIA_XLIM[1], ASIA_XLIM[2], ASIA_YLIM[1], ASIA_YLIM[2])
extent_eu   <- terra::ext(EU_XLIM[1], EU_XLIM[2], 30, EU_YLIM[2])

FIG_WIDTH  <- 10
FIG_HEIGHT <- 7
FIG_DPI    <- 600

# ---- 4. GLOBAL STYLE AND EXPORT FUNCTIONS ----

theme_map_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.30, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0),
      plot.caption = element_text(size = 8, hjust = 0),
      plot.margin = margin(8, 8, 8, 8)
    )
}

north_arrow_layer <- function(size = 0.75) {
  annotation_north_arrow(
    location = "tl",
    pad_x = unit(0.35, "cm"),
    pad_y = unit(0.35, "cm"),
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(size, "cm"),
    width = unit(size, "cm")
  )
}

# Shared styling for the Europe and East Asia map extents.
apply_map_style <- function(p, xlim, ylim, base_size = 11) {
  p +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    north_arrow_layer() +
    theme_map_publication(base_size)
}

save_publication_figure <- function(plot, filename, folder = "MAIN",
                                    width = FIG_WIDTH, height = FIG_HEIGHT,
                                    dpi = FIG_DPI) {
  ggsave(
    filename = file.path(output_root, folder, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

save_csv_table <- function(x, filename) {
  readr::write_csv(x, file.path(output_root, "TABLES", filename))
}

save_xlsx_table <- function(x, filename) {
  writexl::write_xlsx(x, path = file.path(output_root, "TABLES", filename))
}

add_europe_background <- function(fill = "grey96", border = "grey70", linewidth = 0.08) {
  list(
    geom_sf(data = europe, fill = fill, color = NA),
    geom_sf(data = europe, fill = NA, color = border, linewidth = linewidth)
  )
}

# ---- 5. OCCURRENCE DATA CLEANING AND SPATIAL FILTERING ----

# 5.1. Load global country polygons and GBIF occurrences.
world_sf   <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_vect <- terra::vect(world_sf)

gbif_data <- readr::read_tsv(occurrence_file, show_col_types = FALSE)
message("Raw GBIF records: ", nrow(gbif_data))

# 5.2. Remove records with missing or duplicated coordinates.
points_unique <- gbif_data %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

message("After removing missing/duplicate coordinates: ", nrow(points_unique))

points_sf <- st_as_sf(
  points_unique,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

# 5.3. Retain terrestrial records only.
world_union <- st_union(world_sf)
land_index <- lengths(st_intersects(points_sf, world_union)) > 0
points_land <- points_sf[land_index, ]
points_land <- points_land[!duplicated(points_land$geometry), ]

message("After land filtering: ", nrow(points_land))

# 5.4. Restrict occurrence data to the East Asian calibration region.
study_region_sf <- st_crop(
  world_sf,
  xmin = ASIA_XLIM[1], xmax = ASIA_XLIM[2],
  ymin = ASIA_YLIM[1], ymax = ASIA_YLIM[2]
)

points_region <- points_land[
  st_within(points_land, st_union(study_region_sf), sparse = FALSE),
]

message("After restricting to East Asian calibration area: ", nrow(points_region))

# 5.5. Apply 5-km spatial thinning.
thin_data <- data.frame(
  species = "MiniKiwi",
  lat = st_coordinates(points_region)[, 2],
  long = st_coordinates(points_region)[, 1]
)

thinned <- spThin::thin(
  loc.data = thin_data,
  lat.col = "lat",
  long.col = "long",
  spec.col = "species",
  thin.par = 5,
  reps = 10,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  verbose = TRUE
)

pres_final <- as.data.frame(thinned[[1]])
pres_final <- pres_final[, c("Longitude", "Latitude")]
colnames(pres_final) <- c("x", "y")

message("After 5-km spatial thinning: ", nrow(pres_final))

pres_vect <- terra::vect(pres_final, geom = c("x", "y"), crs = "EPSG:4326")

# ---- 6. ENVIRONMENTAL PREDICTORS AND VARIABLE SELECTION ----

# 6.1. Download or load WorldClim predictors.
bio_global  <- geodata::worldclim_global(var = "bio",  res = 2.5, path = "data_env/")
tmin_global <- geodata::worldclim_global(var = "tmin", res = 2.5, path = "data_env/")
tavg_global <- geodata::worldclim_global(var = "tavg", res = 2.5, path = "data_env/")

# 6.2. Prepare East Asian predictor layers.
bio_asia <- crop(bio_global, extent_asia)
names(bio_asia) <- paste0("BIO", 1:19)

spring_min <- app(crop(tmin_global[[3:4]], extent_asia), fun = min) / 10
names(spring_min) <- "spring_min"

chilling_proxy <- app(crop(tavg_global[[c(1, 2, 3, 11, 12)]], extent_asia), fun = mean) / 10
names(chilling_proxy) <- "chilling_proxy"

# 6.3. Screen candidate predictors for multicollinearity.
candidate_stack <- c(bio_asia, spring_min, chilling_proxy)
candidate_vals <- terra::extract(candidate_stack, pres_vect)[, -1]
candidate_vals <- as.data.frame(candidate_vals)
candidate_vals <- candidate_vals[complete.cases(candidate_vals), ]

cor_mat_candidates <- cor(candidate_vals, use = "complete.obs", method = "pearson")

vif_res <- usdm::vifstep(candidate_vals, th = 10)
selected_vars_vif <- vif_res@results$Variables

# 6.4. Define the final predictor set used in modelling.
final_keep <- c("spring_min", "BIO4", "BIO5", "BIO15", "BIO18")
bio_selected <- candidate_stack[[final_keep]]

final_vals <- terra::extract(bio_selected, pres_vect)[, -1]
final_vals <- as.data.frame(final_vals)
final_vals <- final_vals[complete.cases(final_vals), ]

vif_final <- usdm::vifstep(final_vals, th = 10)
cor_mat_final <- cor(final_vals, use = "complete.obs", method = "pearson")

message("Final predictor stack used in modelling:")
print(names(bio_selected))
message("Final correlation matrix among retained predictors:")
print(round(cor_mat_final, 2))

# ---- 7. BACKGROUND SAMPLING AND BIOMOD2 FORMATTING ----

# 7.1. Generate terrestrial background points in East Asia.
asia_land_vect <- crop(world_vect, extent_asia)
asia_land_mask <- rasterize(asia_land_vect, bio_selected[[1]], field = 1, background = NA)

set.seed(123)
bg_points <- terra::spatSample(
  asia_land_mask,
  size = 10000,
  method = "random",
  xy = TRUE,
  na.rm = TRUE
)

bg_points <- as.data.frame(bg_points)[, 1:2]
colnames(bg_points) <- c("x", "y")

resp_xy <- rbind(pres_final, bg_points)
resp_var <- c(rep(1, nrow(pres_final)), rep(0, nrow(bg_points)))

# 7.2. Remove records with missing environmental values.
env_vals <- terra::extract(
  bio_selected,
  terra::vect(resp_xy, geom = c("x", "y"), crs = "EPSG:4326")
)[, -1]

env_vals <- as.data.frame(env_vals)
keep_complete <- complete.cases(env_vals)

resp_xy  <- resp_xy[keep_complete, ]
resp_var <- resp_var[keep_complete]

message("Presences after environmental completeness check: ", sum(resp_var == 1))
message("Background after environmental completeness check: ", sum(resp_var == 0))

# 7.3. Format the response and predictor data for BIOMOD2.
myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp_var,
  expl.var = bio_selected,
  resp.xy = resp_xy,
  resp.name = "MiniKiwi",
  PA.nb.rep = 0,
  filter.raster = TRUE
)

# ---- 8. BIOMOD2 MODEL FITTING, ENSEMBLE MODELLING, AND EUROPEAN PROJECTION ----

# 8.1. Fit individual algorithms using environmental cross-validation.
myBiomodModelOut_env <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = c("GLM", "GAM", "GBM", "RF", "MAXNET"),
  CV.strategy = "env",
  CV.k = 5,
  CV.balance = "presences",
  metric.eval = c("TSS", "ROC"),
  var.import = 10,
  seed.val = 123
)

# 8.2. Extract model evaluation metrics.
evals_env <- get_evaluations(myBiomodModelOut_env)

tss_env <- evals_env %>%
  filter(metric.eval == "TSS") %>%
  select(algo, run, calibration, validation) %>%
  as_tibble()

tss_summary_env <- tss_env %>%
  group_by(algo) %>%
  summarise(
    Mean_TSS = mean(validation, na.rm = TRUE),
    SD_TSS = sd(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_TSS))

roc_env <- evals_env %>%
  filter(metric.eval == "AUCroc") %>%
  select(algo, run, calibration, validation) %>%
  as_tibble()

roc_summary_env <- roc_env %>%
  group_by(algo) %>%
  summarise(
    Mean_ROC = mean(validation, na.rm = TRUE),
    SD_ROC = sd(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_ROC))

# 8.3. Build ensemble models.
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut_env,
  models.chosen = "all",
  em.by = "all",
  em.algo = c("EMmean", "EMwmean"),
  metric.select = "TSS",
  metric.select.thresh = 0.6,
  metric.eval = c("TSS", "ROC"),
  var.import = 3
)

# 8.4. Prepare European projection predictors.
bio_eu_base <- crop(bio_global, extent_eu)
names(bio_eu_base) <- paste0("BIO", 1:19)

spring_eu <- app(crop(tmin_global[[3:4]], extent_eu), fun = min) / 10
names(spring_eu) <- "spring_min"

env_eu_all <- c(bio_eu_base, spring_eu)
env_eu <- env_eu_all[[names(bio_selected)]]

if (!identical(names(bio_selected), names(env_eu))) {
  stop("Training and projection predictor names do not match.")
}

# 8.5. Prepare March minimum temperature as a visualization-only frost layer.
tmin_march_eu <- crop(tmin_global[[3]], extent_eu) / 10
frost_df <- as.data.frame(tmin_march_eu, xy = TRUE, na.rm = FALSE)
colnames(frost_df) <- c("x", "y", "MarchTmin")

# 8.6. Project individual models to Europe.
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut_env,
  new.env = env_eu,
  proj.name = "Europe",
  models.chosen = "all",
  metric.binary = "TSS",
  compress = "xz",
  build.clamping.mask = TRUE
)

# 8.7. Project ensemble models to Europe.
myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  metric.binary = "TSS"
)

all_ensemble_predictions <- get_predictions(myBiomodEMProj)
message("Available ensemble layers:")
print(names(all_ensemble_predictions))

suit_map_main <- all_ensemble_predictions[["MiniKiwi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo"]]
suit_map_emmean <- all_ensemble_predictions[["MiniKiwi_EMmeanByTSS_mergedData_mergedRun_mergedAlgo"]]

suit_df <- as.data.frame(suit_map_main, xy = TRUE, na.rm = TRUE)
colnames(suit_df) <- c("x", "y", "Suitability")
if (max(suit_df$Suitability, na.rm = TRUE) > 1) {
  suit_df$Suitability <- suit_df$Suitability / 1000
}

suit_df_emmean <- as.data.frame(suit_map_emmean, xy = TRUE, na.rm = TRUE)
colnames(suit_df_emmean) <- c("x", "y", "Suitability")
if (max(suit_df_emmean$Suitability, na.rm = TRUE) > 1) {
  suit_df_emmean$Suitability <- suit_df_emmean$Suitability / 1000
}

# 8.8. Load Europe polygons for all map outputs.
europe <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe <- europe[europe$continent == "Europe", ]

# ---- 9. DERIVED MAP LAYERS ----

# Frost-risk classes for a given high-risk threshold (°C).
classify_frost_risk <- function(df, threshold) {
  hi_lab  <- sprintf("High risk (≤ %s°C)", threshold)
  mod_lab <- sprintf("Moderate risk (%s to 0°C)", threshold)
  low_lab <- "Low warning (0 to 1°C)"

  df %>%
    filter(!is.na(MarchTmin)) %>%
    mutate(
      Risk = case_when(
        MarchTmin <= threshold ~ hi_lab,
        MarchTmin > threshold & MarchTmin <= 0 ~ mod_lab,
        MarchTmin > 0 & MarchTmin <= 1 ~ low_lab,
        TRUE ~ "No risk"
      )
    ) %>%
    filter(Risk != "No risk") %>%
    mutate(Risk = factor(Risk, levels = c(hi_lab, mod_lab, low_lab)))
}

frost_10_df <- classify_frost_risk(frost_df, -1)
frost_15_df <- classify_frost_risk(frost_df, -1.5)

# Uncertainty layer from individual projected models.
model_stack <- get_predictions(myBiomodProj)
if (terra::global(model_stack, "max", na.rm = TRUE)[1, 1] > 1) {
  model_stack <- model_stack / 1000
}

sd_map <- app(model_stack, fun = sd, na.rm = TRUE)
names(sd_map) <- "SD_uncertainty"

sd_eur <- mask(crop(sd_map, terra::vect(europe)), terra::vect(europe))
sd_df <- as.data.frame(sd_eur, xy = TRUE, na.rm = TRUE)
colnames(sd_df) <- c("x", "y", "SD")

# Save main raster outputs.
raster_exports <- list(
  suitability_EMwmean.tif = suit_map_main,
  suitability_EMmean.tif = suit_map_emmean,
  uncertainty_SD.tif = sd_eur,
  march_minimum_temperature.tif = tmin_march_eu
)
for (nm in names(raster_exports)) {
  terra::writeRaster(
    raster_exports[[nm]],
    filename = file.path(output_root, "RASTERS", nm),
    overwrite = TRUE
  )
}

# ---- 10. MAIN FIGURES ----

# 10.1. Figure 1: Study area and cleaned/thinned occurrences in East Asia.
points_region_plot <- points_region
points_thinned_plot <- st_as_sf(pres_final, coords = c("x", "y"), crs = 4326)

p_fig01 <- ggplot() +
  geom_sf(data = study_region_sf, fill = "grey96", color = "grey70", linewidth = 0.20) +
  geom_sf(data = points_region_plot, color = "grey40", size = 0.9, alpha = 0.55) +
  geom_sf(data = points_thinned_plot, color = "#990000", size = 1.0, alpha = 0.78) +
  labs(
    title = expression(bold("Occurrence records of ") * bolditalic("Actinidia arguta") * bold(" in East Asia")),
    subtitle = "Cleaned records before thinning (grey) and 5-km thinned records used for modelling (red)",
    caption = "Occurrence records were filtered for missing coordinates, duplicate coordinates, land overlap, and calibration-area extent.",
    x = "Longitude",
    y = "Latitude"
  )

p_fig01 <- apply_map_style(p_fig01, ASIA_XLIM, ASIA_YLIM)
save_publication_figure(p_fig01, "Figure_01_EastAsia_occurrences.png", "MAIN", width = 8, height = 6)

# 10.2-10.3. Figures 2 & 3: predictor correlation heatmaps (candidate set and final set).

make_corr_heatmap <- function(cor_mat, title, subtitle, label_threshold,
                               tile_color, tile_linewidth,
                               low, mid, high, text_size, base_size,
                               title_size, subtitle_size, axis_text_size,
                               legend_title_size, legend_text_size) {
  cm <- cor_mat
  cm[upper.tri(cm)] <- NA
  cm_long <- reshape2::melt(cm, na.rm = TRUE)
  colnames(cm_long) <- c("Var1", "Var2", "Correlation")

  cm_long$Label <- if (is.null(label_threshold)) {
    sprintf("%.2f", cm_long$Correlation)
  } else {
    ifelse(
      abs(cm_long$Correlation) >= label_threshold & cm_long$Correlation != 1,
      sprintf("%.2f", cm_long$Correlation), ""
    )
  }

  ggplot(cm_long, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = tile_color, linewidth = tile_linewidth) +
    geom_text(aes(label = Label), size = text_size, color = "black") +
    scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0,
                         limits = c(-1, 1), name = "Pearson r") +
    coord_fixed() +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      plot.title = element_text(face = "bold", size = title_size),
      plot.subtitle = element_text(size = subtitle_size),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size)
    )
}

p_fig02 <- make_corr_heatmap(
  cor_mat_candidates,
  title = "Correlation heatmap of candidate predictors",
  subtitle = "Bioclimatic variables and derived thermal predictors before final selection",
  label_threshold = 0.75,
  tile_color = "grey92", tile_linewidth = 0.20,
  low = "#3b4cc0", mid = "white", high = "#d73027",
  text_size = 2.6, base_size = 9,
  title_size = 13, subtitle_size = 9, axis_text_size = 7,
  legend_title_size = 9, legend_text_size = 7
)
save_publication_figure(p_fig02, "Figure_02_Candidate_predictor_correlations.png", "MAIN", width = 7.5, height = 7)

p_fig03 <- make_corr_heatmap(
  cor_mat_final,
  title = "Correlation matrix of the final predictor set",
  subtitle = "Reduced multicollinearity after VIF filtering and ecological selection",
  label_threshold = NULL,
  tile_color = "white", tile_linewidth = 0.30,
  low = "#5e7bb6", mid = "#f7f7f7", high = "#c75b5b",
  text_size = 3.9, base_size = 11,
  title_size = 16, subtitle_size = 12, axis_text_size = 10,
  legend_title_size = 11, legend_text_size = 10
)
save_publication_figure(p_fig03, "Figure_03_Final_predictor_correlations.png", "MAIN", width = 7.5, height = 6.5)

# 10.4. Figure 4: ROC curve for GLM environmental cross-validation runs.

roc_plot_created <- FALSE

tryCatch({
  preds <- get_predictions(myBiomodModelOut_env)
  obs <- myBiomodData@data.species
  calib_lines <- as.data.frame(get_calib_lines(myBiomodModelOut_env))

  # Exclude "allRun": fitted on the full dataset, no independent validation subset.
  glm_dat <- preds %>%
    filter(
      algo == "GLM",
      run != "allRun"
    ) %>%
    mutate(obs = obs[points]) %>%
    filter(
      !is.na(obs),
      !is.na(pred)
    )

  run_ids <- unique(glm_dat$run)

  roc_list <- lapply(run_ids, function(run_id) {
    df_run <- glm_dat %>%
      filter(run == run_id)

    calib_col <- paste0("_allData_", run_id)
    if (!calib_col %in% names(calib_lines)) {
      stop("Could not find calibration column for ", run_id)
    }

    # FALSE = validation observations in BIOMOD2 calibration table
    valid_points <- which(calib_lines[[calib_col]] == FALSE)
    df_valid <- df_run %>%
      filter(points %in% valid_points)

    if (nrow(df_valid) == 0) {
      stop("No validation observations found for ", run_id)
    }
    if (dplyr::n_distinct(df_valid$obs) < 2) {
      stop("Validation subset does not contain both classes for ", run_id)
    }

    pROC::roc(
      response = df_valid$obs,
      predictor = df_valid$pred,
      quiet = TRUE
    )
  })

  names(roc_list) <- run_ids

  # Common FPR grid for averaging ROC curves
  fpr_grid <- seq(0, 1, length.out = 300)

  interp_mat <- sapply(roc_list, function(r) {
    tmp <- data.frame(
      fpr = 1 - r$specificities,
      tpr = r$sensitivities
    )
    tmp <- tmp[order(tmp$fpr), ]
    tmp <- tmp[!duplicated(tmp$fpr), ]

    approx(
      x = tmp$fpr,
      y = tmp$tpr,
      xout = fpr_grid,
      rule = 2
    )$y
  })

  roc_summary_df <- data.frame(
    fpr = fpr_grid,
    mean_tpr = rowMeans(interp_mat, na.rm = TRUE),
    lower_tpr = apply(interp_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    upper_tpr = apply(interp_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )

  # Use BIOMOD2 validation AUC values so Figure 4 agrees with Table 2
  glm_auc_summary <- roc_summary_env %>%
    filter(algo == "GLM")

  auc_label <- paste0(
    "Mean validation AUC = ",
    sprintf("%.3f", glm_auc_summary$Mean_ROC[1]),
    " ± ",
    sprintf("%.3f", glm_auc_summary$SD_ROC[1])
  )

  p_fig04 <- ggplot(roc_summary_df, aes(x = fpr, y = mean_tpr)) +
    geom_ribbon(aes(ymin = lower_tpr, ymax = upper_tpr), fill = "#3182bd", alpha = 0.25) +
    geom_line(color = "#b30000", linewidth = 1.15) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
    annotate("text", x = 0.64, y = 0.14, label = auc_label, size = 4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      title = "ROC curve for the GLM model",
      subtitle = "Environmental cross-validation runs",
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )

  save_publication_figure(p_fig04, "Figure_04_ROC_curve_GLM.png", "MAIN", width = 7.5, height = 6.5)

  roc_plot_created <- TRUE
  message("Figure 4 ROC plot created successfully using GLM environmental validation runs.")
}, error = function(e) {
  message("Figure 4 ROC plot failed: ", conditionMessage(e))
})

if (!roc_plot_created) {
  message("Figure 4 was not exported. Check the error message above.")
}

# 10.5. Figure 5: Model performance based on TSS.

tss_runs_env <- tss_env %>%
  rename(Model = algo, Run = run, TSS = validation) %>%
  filter(!is.na(TSS)) %>%
  mutate(Model = factor(Model, levels = c("RF", "MAXNET", "GLM", "GBM", "GAM")))

tss_summary_for_plot <- tss_runs_env %>%
  group_by(Model) %>%
  summarise(Mean_TSS = mean(TSS, na.rm = TRUE), .groups = "drop")

p_fig05 <- ggplot(tss_runs_env, aes(x = TSS, y = Model)) +
  geom_boxplot(
    aes(fill = Model),
    color = "grey25",
    linewidth = 0.45,
    outlier.shape = 16,
    outlier.size = 1.6,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_point(data = tss_summary_for_plot, aes(x = Mean_TSS, y = Model), color = "#b30000", size = 2.5) +
  scale_fill_grey(start = 0.85, end = 0.45) +
  labs(
    title = "Model performance (TSS per run) - ENV",
    x = "TSS",
    y = "Model"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

save_publication_figure(p_fig05, "Figure_05_Model_performance_TSS.png", "MAIN", width = 7, height = 4.8)


# 10.6. Figure 6: Relative importance of environmental predictors.

var_imp_df_env <- as.data.frame(get_variables_importance(myBiomodModelOut_env))

var_imp_summary_env <- var_imp_df_env %>%
  group_by(expl.var) %>%
  summarise(
    MeanImportance = mean(var.imp, na.rm = TRUE),
    SDImportance = sd(var.imp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(MeanImportance), MeanImportance > 0) %>%
  arrange(desc(MeanImportance))

p_fig06 <- ggplot(
  var_imp_summary_env,
  aes(
    x = reorder(expl.var, MeanImportance),
    y = MeanImportance,
    fill = MeanImportance
  )
) +
  geom_col(color = "grey25", linewidth = 0.25, width = 0.72) +
  coord_flip() +
  scale_fill_gradient(low = "grey80", high = "grey35", guide = "none") +
  labs(
    title = "Relative importance of environmental predictors",
    subtitle = "Mean permutation importance across models and runs",
    x = "Predictor",
    y = "Mean importance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10)
  )

save_publication_figure(p_fig06, "Figure_06_Variable_importance.png", "MAIN", width = 7, height = 4.8)

# 10.7. Figure 7: Climatic suitability map.
p_fig07 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_raster(
    data = suit_df,
    aes(x = x, y = y, fill = Suitability),
    interpolate = TRUE
  ) +
  geom_sf(data = europe, fill = NA, color = "grey70", linewidth = 0.06) +
  scale_fill_viridis_c(
    option = "D",
    begin = 0.08,
    end = 0.92,
    direction = 1,
    limits = c(0, 1),
    name = "Climatic suitability",
    guide = guide_colorbar(barheight = unit(5.5, "cm"))
  ) +
  labs(
    title = expression(bold("Climatic suitability for ") * bolditalic("Actinidia arguta") * bold(" in Europe")),
    subtitle = "Ensemble model (TSS-weighted mean, EMwmean)",
    caption = "Values represent continuous ensemble climatic suitability scaled from 0 to 1.",
    x = "Longitude",
    y = "Latitude"
  )

p_fig07 <- apply_map_style(p_fig07, EU_XLIM, EU_YLIM)
save_publication_figure(p_fig07, "Figure_07_Climatic_suitability_EMwmean.png", "MAIN")

# 10.8. Figure 8: Spatial uncertainty map.
p_fig08 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_raster(
    data = sd_df,
    aes(x = x, y = y, fill = SD),
    interpolate = TRUE
  ) +
  geom_sf(data = europe, fill = NA, color = "grey70", linewidth = 0.06) +
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    name = "Prediction uncertainty (SD)",
    guide = guide_colorbar(barheight = unit(5.5, "cm"))
  ) +
  labs(
    title = expression(bold("Spatial uncertainty of ") * bolditalic("Actinidia arguta") * bold(" projections in Europe")),
    subtitle = "Standard deviation across projected individual model predictions",
    caption = "Higher values indicate stronger disagreement among individual model projections.",
    x = "Longitude",
    y = "Latitude"
  )

p_fig08 <- apply_map_style(p_fig08, EU_XLIM, EU_YLIM)
save_publication_figure(p_fig08, "Figure_08_Spatial_uncertainty.png", "MAIN")

# 10.9. Figure 9: Bivariate suitability and uncertainty map.

# Shared draw for a classified (bi_class already applied) suitability x y_col bivariate map.
make_bivariate_map <- function(df, y_lab, title, caption,
                                border_color, border_linewidth, arrow_size) {
  map <- ggplot() +
    geom_sf(data = europe, fill = "grey97", color = "grey75", linewidth = 0.20) +
    geom_raster(data = df, aes(x = x, y = y, fill = bi_class), show.legend = FALSE) +
    biscale::bi_scale_fill(pal = "DkBlue2", dim = 3, flip_axes = FALSE, rotate_pal = FALSE) +
    geom_sf(data = europe, fill = NA, color = border_color, linewidth = border_linewidth) +
    coord_sf(xlim = EU_XLIM, ylim = EU_YLIM, expand = FALSE) +
    north_arrow_layer(arrow_size) +
    labs(
      title = title,
      subtitle = expression(italic("Actinidia arguta")),
      caption = caption,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_map_publication(10) +
    theme(legend.position = "none")

  legend <- biscale::bi_legend(
    pal = "DkBlue2", dim = 3,
    xlab = "Suitability →", ylab = paste0(y_lab, " ↑"),
    size = 8, flip_axes = FALSE, rotate_pal = FALSE
  )

  cowplot::ggdraw() +
    cowplot::draw_plot(map, 0, 0, 1, 1) +
    cowplot::draw_plot(legend, 0.075, 0.095, 0.19, 0.19)
}

biv_unc_df <- merge(
  suit_df[, c("x", "y", "Suitability")],
  sd_df[, c("x", "y", "SD")],
  by = c("x", "y"),
  all = FALSE
) %>%
  filter(!is.na(Suitability), !is.na(SD)) %>%
  distinct(x, y, .keep_all = TRUE)

biv_unc_df <- biscale::bi_class(biv_unc_df, x = Suitability, y = SD, style = "quantile", dim = 3)

p_fig09 <- make_bivariate_map(
  biv_unc_df, "Uncertainty",
  title = "Bivariate climatic suitability and prediction uncertainty in Europe",
  caption = "Bivariate classification based on EMwmean suitability and standard deviation across projected individual models.",
  border_color = "grey55", border_linewidth = 0.18, arrow_size = 0.60
)

save_publication_figure(p_fig09, "Figure_09_Bivariate_suitability_uncertainty.png", "MAIN")

# 10.10. Figure 10: Early spring frost risk map under the -1.0 °C threshold.
frost_cols_10 <- c(
  "High risk (≤ -1°C)" = "#b30000",
  "Moderate risk (-1 to 0°C)" = "#f46d43",
  "Low warning (0 to 1°C)" = "#fee8c8"
)

p_fig10 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_tile(data = frost_10_df, aes(x = x, y = y, fill = Risk), alpha = 0.88) +
  geom_sf(data = europe, fill = NA, color = "grey60", linewidth = 0.08) +
  scale_fill_manual(values = frost_cols_10, name = "Early spring frost risk") +
  labs(
    title = "Early spring frost risk in Europe",
    subtitle = "Less conservative threshold scheme: high risk ≤ -1°C",
    caption = "Frost-risk classes derived from mean March minimum temperature.",
    x = "Longitude",
    y = "Latitude"
  )

p_fig10 <- apply_map_style(p_fig10, EU_XLIM, EU_YLIM)
save_publication_figure(p_fig10, "Figure_10_Frost_risk_threshold_minus1.png", "MAIN")

# 10.11. Figure 11: Suitability and frost-risk overlay under the -1.5 °C threshold.
make_overlay_map <- function(threshold = -1.5, moderate_label, title_suffix) {
  frost_moderate <- frost_df %>%
    filter(!is.na(MarchTmin), MarchTmin > threshold, MarchTmin <= 0) %>%
    mutate(fill_group = "Moderate frost risk")
  
  highrisk_df <- frost_df %>%
    filter(!is.na(MarchTmin))

  p <- ggplot() +
    geom_sf(data = europe, fill = "grey96", color = NA) +
    geom_tile(data = suit_df, aes(x = x, y = y, fill = Suitability)) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      name = "Climatic suitability",
      guide = guide_colorbar(order = 1, barheight = unit(5.0, "cm"))
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = frost_moderate,
      aes(x = x, y = y, fill = fill_group),
      alpha = 0.18
    ) +
    scale_fill_manual(
      values = c("Moderate frost risk" = "#fdae61"),
      name = "Frost-risk zones",
      guide = guide_legend(order = 2, override.aes = list(alpha = 0.5))
    ) +
    geom_contour(
    data = highrisk_df,
    aes(x = x, y = y, z = MarchTmin),
    breaks = threshold,
    binwidth = NULL,
    color = "#cc0000",
    linewidth = 0.22,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
    
    geom_sf(data = europe, fill = NA, color = "grey70", linewidth = 0.06) +
    labs(
      title = expression(bold("Climatic suitability and early spring frost risk for ") * bolditalic("Actinidia arguta") * bold(" in Europe")),
      subtitle = title_suffix,
      caption = paste0(
        "Ensemble model (EMwmean). Moderate frost risk: ",
        moderate_label,
        "; high-risk contour: ≤ ",
        threshold,
        "°C."
      ),
      x = "Longitude",
      y = "Latitude"
    )

  apply_map_style(p, EU_XLIM, EU_YLIM)
}

p_fig11 <- make_overlay_map(
  threshold = -1.5,
  moderate_label = "-1.5 to 0°C",
  title_suffix = "Intermediate threshold scheme: high risk ≤ -1.5°C"
)

save_publication_figure(p_fig11, "Figure_11_Overlay_suitability_frost_minus1p5.png", "MAIN")

# 10.12. Figure 12: Suitability and frost-risk overlay under the -1.0 °C threshold.
p_fig12 <- make_overlay_map(
  threshold = -1,
  moderate_label = "-1 to 0°C",
  title_suffix = "Less conservative threshold scheme: high risk ≤ -1°C"
)

save_publication_figure(p_fig12, "Figure_12_Overlay_suitability_frost_minus1.png", "MAIN")

# 10.13. Figure 13: GLM response curves.
response_plot_created <- FALSE
try({
  glm_models <- grep("_GLM$", get_built_models(myBiomodModelOut_env), value = TRUE)
  glm_models <- glm_models[!grepl("allRun", glm_models)]

  rc_glm <- bm_PlotResponseCurves(
    bm.out = myBiomodModelOut_env,
    models.chosen = glm_models,
    fixed.var = "median",
    show.variables = names(bio_selected)
  )

  df_glm <- rc_glm$tab

  response_summary_glm <- df_glm %>%
    group_by(expl.name, expl.val) %>%
    summarise(
      mean_pred = mean(pred.val, na.rm = TRUE),
      sd_pred = sd(pred.val, na.rm = TRUE),
      n = dplyr::n(),
      se_pred = sd_pred / sqrt(n),
      lower = pmax(mean_pred - 1.96 * se_pred, 0),
      upper = pmin(mean_pred + 1.96 * se_pred, 1),
      .groups = "drop"
    ) %>%
    mutate(
      expl.name = factor(
        expl.name,
        levels = c("spring_min", "BIO4", "BIO5", "BIO15", "BIO18"),
        labels = c(
          "Spring minimum temperature",
          "Temperature seasonality (BIO4)",
          "Max temperature of warmest month (BIO5)",
          "Precipitation seasonality (BIO15)",
          "Precipitation of warmest quarter (BIO18)"
        )
      )
    )

  p_fig13 <- ggplot(response_summary_glm, aes(x = expl.val, y = mean_pred)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#8da0ff",
      alpha = 0.18
    ) +
    
    geom_line(
      color = "#d7191c",
      linewidth = 1.15,
      lineend = "round"
    ) +
    facet_wrap(~ expl.name, scales = "free_x", ncol = 2) +
    labs(
      title = "Response curves for selected environmental predictors",
      subtitle = "GLM environmental cross-validation runs: mean response with 95% confidence interval",
      x = NULL,
      y = "Predicted suitability"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )

  save_publication_figure(p_fig13, "Figure_13_Response_curves_GLM.png", "MAIN", height = 8)
  response_plot_created <- TRUE
}, silent = TRUE)

if (!response_plot_created) {
  message("Figure 13 response curves were skipped because bm_PlotResponseCurves returned an unexpected structure.")
}

# ---- 11. SUPPLEMENTARY FIGURES ----

# 11.1. Supplementary Figure S1: Bias surface from cleaned occurrence records before thinning.
asia_sf <- st_crop(
  world_sf,
  xmin = ASIA_XLIM[1], xmax = ASIA_XLIM[2],
  ymin = ASIA_YLIM[1], ymax = ASIA_YLIM[2]
)

bias_template <- rast(
  xmin = ASIA_XLIM[1], xmax = ASIA_XLIM[2],
  ymin = ASIA_YLIM[1], ymax = ASIA_YLIM[2],
  resolution = 0.25,
  crs = "EPSG:4326"
)

points_vect_bias <- terra::vect(points_region)
bias_count <- terra::rasterize(points_vect_bias, bias_template, fun = "sum", background = 0)

bias_smooth <- terra::focal(
  bias_count,
  w = matrix(1, 5, 5),
  fun = sum,
  na.policy = "omit"
)

bias_df <- as.data.frame(bias_smooth, xy = TRUE, na.rm = TRUE)
colnames(bias_df) <- c("x", "y", "Bias")
bias_df$Bias_plot <- log1p(bias_df$Bias)

p_s1 <- ggplot() +
  geom_sf(data = asia_sf, fill = "grey96", color = "grey70", linewidth = 0.20) +
  geom_raster(
    data = bias_df,
    aes(x = x, y = y, fill = Bias_plot),
    alpha = 0.75,
    interpolate = TRUE
  ) +
  geom_sf(data = asia_sf, fill = NA, color = "grey72", linewidth = 0.16) +
  scale_fill_gradientn(
    colors = c("white", "#fee8c8", "#fdbb84", "#e34a33"),
    name = "Log-transformed\nsampling density"
  ) +
  labs(
    title = "Spatial sampling density of occurrence records in East Asia",
    subtitle = "Bias surface derived from cleaned occurrence records before thinning",
    caption = "Rasterized occurrence counts were smoothed using a 5 × 5 moving window.",
    x = "Longitude",
    y = "Latitude"
  )

p_s1 <- apply_map_style(p_s1, ASIA_XLIM, ASIA_YLIM)
save_publication_figure(p_s1, "Figure_S1_Bias_surface_EastAsia.png", "SUPPLEMENTARY", width = 8.5, height = 6.5)

# 11.2. Supplementary Figure S2: EMmean suitability sensitivity map.
p_s2 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_tile(data = suit_df_emmean, aes(x = x, y = y, fill = Suitability)) +
  geom_sf(data = europe, fill = NA, color = "grey74", linewidth = 0.05) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0, 1),
    name = "Climatic suitability",
    guide = guide_colorbar(barheight = unit(5.2, "cm"))
  ) +
  labs(
    title = expression(bold("Climatic suitability for ") * bolditalic("Actinidia arguta") * bold(" in Europe")),
    subtitle = "Sensitivity analysis of ensemble predictions (EMmean)",
    caption = "This map is used as a sensitivity comparison against the main TSS-weighted ensemble map.",
    x = "Longitude",
    y = "Latitude"
  )

p_s2 <- apply_map_style(p_s2, EU_XLIM, EU_YLIM)
save_publication_figure(p_s2, "Figure_S2_Climatic_suitability_EMmean.png", "SUPPLEMENTARY")

# 11.3. Supplementary Figure S3: Bivariate suitability and frost risk.
biv_frost_df <- merge(
  suit_df[, c("x", "y", "Suitability")],
  frost_df[, c("x", "y", "MarchTmin")],
  by = c("x", "y"),
  all = FALSE
) %>%
  filter(!is.na(Suitability), !is.na(MarchTmin)) %>%
  mutate(FrostRisk = -MarchTmin)

biv_frost_df <- biscale::bi_class(biv_frost_df, x = Suitability, y = FrostRisk, style = "quantile", dim = 3)

p_s3 <- make_bivariate_map(
  biv_frost_df, "Frost risk",
  title = "Bivariate climatic suitability and spring frost risk in Europe",
  caption = "Bivariate classification based on EMwmean suitability and mean March minimum temperature.",
  border_color = "grey40", border_linewidth = 0.28, arrow_size = 0.70
)

save_publication_figure(p_s3, "Figure_S3_Bivariate_suitability_frost_risk.png", "SUPPLEMENTARY")

# 11.4. Supplementary Figure S4: Suitability > 0.5.
suit_df_05 <- suit_df %>% filter(Suitability > 0.5)

p_s4 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_tile(data = suit_df_05, aes(x = x, y = y, fill = Suitability)) +
  geom_sf(data = europe, fill = NA, color = "grey70", linewidth = 0.06) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0.5, 1),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = scales::number_format(accuracy = 0.1),
    oob = scales::squish,
    name = "Climatic suitability"
  ) +
  labs(
    title = "Areas with climatic suitability > 0.5 in Europe",
    subtitle = expression(italic("Actinidia arguta")),
    caption = "Only pixels with suitability > 0.5 are displayed.",
    x = "Longitude",
    y = "Latitude"
  )

p_s4 <- apply_map_style(p_s4, EU_XLIM, EU_YLIM)
save_publication_figure(p_s4, "Figure_S4_Suitability_greater_than_0p5.png", "SUPPLEMENTARY")

# 11.5. Supplementary Figure S5: Suitability > 0.8.
suit_df_08 <- suit_df %>% filter(Suitability > 0.8)

p_s5 <- ggplot() +
  geom_sf(data = europe, fill = "grey96", color = NA) +
  geom_tile(data = suit_df_08, aes(x = x, y = y, fill = Suitability)) +
  geom_sf(data = europe, fill = NA, color = "grey70", linewidth = 0.06) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0.8, 1),
    breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
    labels = scales::number_format(accuracy = 0.01),
    oob = scales::squish,
    name = "Climatic suitability"
  ) +
  labs(
    title = "Areas with climatic suitability > 0.8 in Europe",
    subtitle = expression(italic("Actinidia arguta")),
    caption = "Only pixels with suitability > 0.8 are displayed.",
    x = "Longitude",
    y = "Latitude"
  )

p_s5 <- apply_map_style(p_s5, EU_XLIM, EU_YLIM)
save_publication_figure(p_s5, "Figure_S5_Suitability_greater_than_0p8.png", "SUPPLEMENTARY")

# ---- 12. COUNTRY-LEVEL SUMMARIES AND SUPPLEMENTARY COUNTRY MAPS ----

# 12.1. Calculate country-level suitability summaries.
suit_rast <- suit_map_main
if (terra::global(suit_rast, "max", na.rm = TRUE)[1, 1] > 1) {
  suit_rast <- suit_rast / 1000
}

europe_sf <- europe[, c("admin")]
europe_vect <- terra::vect(europe_sf)

suit_eur <- mask(crop(suit_rast, europe_vect), europe_vect)
country_vals <- terra::extract(suit_eur, europe_vect)
country_vals <- as.data.frame(country_vals)

names(country_vals)[1] <- "ID"
names(country_vals)[2] <- "Suitability"

country_vals$country <- europe_sf$admin[country_vals$ID]

country_vals <- country_vals %>%
  filter(!is.na(Suitability))

country_summary <- country_vals %>%
  group_by(country) %>%
  summarise(
    max_suitability = max(Suitability, na.rm = TRUE),
    mean_suitability = mean(Suitability, na.rm = TRUE),
    cells_total = n(),
    cells_gt_05 = sum(Suitability > 0.5, na.rm = TRUE),
    cells_gt_08 = sum(Suitability > 0.8, na.rm = TRUE),
    pct_gt_05 = 100 * cells_gt_05 / cells_total,
    pct_gt_08 = 100 * cells_gt_08 / cells_total,
    .groups = "drop"
  ) %>%
  arrange(desc(max_suitability))

country_classification_05 <- country_summary %>%
  mutate(
    suitability_class_05 = case_when(
      pct_gt_05 > 30 ~ "High suitability extent",
      pct_gt_05 > 5 ~ "Moderate suitability extent",
      TRUE ~ "Marginal suitability extent"
    )
  ) %>%
  arrange(desc(pct_gt_05))

country_classification_08 <- country_summary %>%
  mutate(
    suitability_class_08 = case_when(
      pct_gt_08 > 10 ~ "High core suitability extent",
      pct_gt_08 > 1 ~ "Moderate core suitability extent",
      TRUE ~ "Marginal core suitability extent"
    )
  ) %>%
  arrange(desc(pct_gt_08))

save_csv_table(country_summary, "Country_suitability_summary.csv")
save_csv_table(country_classification_05, "Country_suitability_classification_gt_0p5.csv")
save_csv_table(country_classification_08, "Country_suitability_classification_gt_0p8.csv")

save_xlsx_table(
  list(
    Summary = country_summary,
    Classification_gt_0p5 = country_classification_05,
    Classification_gt_0p8 = country_classification_08
  ),
  "Country_suitability_results.xlsx"
)

# 12.2. Create country-level threshold maps.
europe_country <- europe %>%
  left_join(country_summary, by = c("admin" = "country"))

make_country_max_map <- function(threshold) {
  map_df <- europe_country %>%
    mutate(show_country = !is.na(max_suitability) & max_suitability >= threshold)

  p <- ggplot() +
    geom_sf(data = map_df, fill = "white", color = "grey70", linewidth = 0.15) +
    geom_sf(
      data = subset(map_df, show_country),
      aes(fill = max_suitability),
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = 1,
      limits = c(0.5, 1),
      breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      name = "Country suitability\n(max)",
      guide = guide_colorbar(barheight = unit(5.5, "cm"))
    ) +
    labs(
      title = paste0("Countries with maximum suitability ≥ ", threshold, " in Europe"),
      subtitle = expression(italic("Actinidia arguta")),
      caption = "Country suitability is mapped as the maximum predicted suitability within each country.",
      x = "Longitude",
      y = "Latitude"
    )

  apply_map_style(p, EU_XLIM, EU_YLIM)
}

p_s6 <- make_country_max_map(0.5)
p_s7 <- make_country_max_map(0.7)

save_publication_figure(p_s6, "Figure_S6_Countries_max_suitability_ge_0p5.png", "SUPPLEMENTARY")
save_publication_figure(p_s7, "Figure_S7_Countries_max_suitability_ge_0p7.png", "SUPPLEMENTARY")

p_s6_panel <- p_s6 +
  theme(legend.position = "none") +
  labs(title = "A. Maximum suitability ≥ 0.5")

p_s7_panel <- p_s7 +
  theme(legend.position = "right") +
  labs(title = "B. Maximum suitability ≥ 0.7")

p_s8 <- p_s6_panel + p_s7_panel + patchwork::plot_layout(ncol = 2, widths = c(1, 1.08))

save_publication_figure(
  p_s8,
  "Figure_S8_Country_suitability_thresholds_two_panel.png",
  "SUPPLEMENTARY",
  width = 12.5,
  height = 6.2
)

# ---- 12.3. Supplementary Figure S9: Maximum suitability by country ----

country_summary_plot <- country_summary %>%
  filter(!is.na(max_suitability)) %>%
  arrange(desc(max_suitability))

p_s9 <- ggplot(
  country_summary_plot,
  aes(
    x = reorder(country, max_suitability),
    y = max_suitability
  )
) +
  geom_col(
    fill = "#3b6fb6",
    width = 0.72
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = expression(
      bold("Maximum climatic suitability by country for ") *
        bolditalic("Actinidia arguta")
    ),
    subtitle = "Maximum predicted suitability extracted from EMwmean raster",
    x = NULL,
    y = "Maximum suitability",
    caption = "Countries are ranked by maximum predicted climatic suitability."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_publication_figure(
  p_s9,
  "Figure_S9_Max_suitability_by_country.png",
  "SUPPLEMENTARY",
  width = 8.5,
  height = 7
)

# ---- 12.4. Supplementary Figure S10: Mean suitability by country ----

country_summary_plot2 <- country_summary %>%
  filter(!is.na(mean_suitability)) %>%
  arrange(desc(mean_suitability))

p_s10 <- ggplot(
  country_summary_plot2,
  aes(
    x = reorder(country, mean_suitability),
    y = mean_suitability
  )
) +
  geom_col(
    fill = "#c00000",
    width = 0.72
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = expression(
      bold("Mean climatic suitability by country for ") *
        bolditalic("Actinidia arguta")
    ),
    subtitle = "Mean predicted suitability extracted from EMwmean raster",
    x = NULL,
    y = "Mean suitability",
    caption = "Countries are ranked by mean predicted climatic suitability."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_publication_figure(
  p_s10,
  "Figure_S10_Mean_suitability_by_country.png",
  "SUPPLEMENTARY",
  width = 8.5,
  height = 7
)

# ---- 13. TABLE EXPORTS ----

csv_exports <- list(
  TSS_per_run_ENV.csv = tss_env,
  TSS_summary_ENV.csv = tss_summary_env,
  ROC_per_run_ENV.csv = roc_env,
  ROC_summary_ENV.csv = roc_summary_env,
  Variable_importance_ENV.csv = var_imp_summary_env,
  VIF_candidate_predictors.csv = tibble::as_tibble(vif_res@results),
  VIF_final_predictors.csv = tibble::as_tibble(vif_final@results),
  Correlation_matrix_candidate_predictors.csv =
    tibble::rownames_to_column(as.data.frame(cor_mat_candidates), "Variable"),
  Correlation_matrix_final_predictors.csv =
    tibble::rownames_to_column(as.data.frame(cor_mat_final), "Variable")
)
for (nm in names(csv_exports)) save_csv_table(csv_exports[[nm]], nm)

# ---- 14. FINAL STATUS MESSAGE ----

message("=================================================================")
message("FINAL EXPORT SCRIPT COMPLETED SUCCESSFULLY")
message("Main figures:          ", file.path(output_root, "MAIN"))
message("Supplementary figures: ", file.path(output_root, "SUPPLEMENTARY"))
message("Tables:                ", file.path(output_root, "TABLES"))
message("Rasters:               ", file.path(output_root, "RASTERS"))
message("=================================================================")
