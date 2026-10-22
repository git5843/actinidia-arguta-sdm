# ==============================================================================
# RUN_THIS_TOMORROW_Actinidia_revision.R
# Actinidia arguta - reviewer revision analysis
#
# PURPOSE
#   1) Reproduce occurrence cleaning with transparent counts.
#   2) Correct WorldClim temperature scaling (NO /10 for WorldClim 2.1 tmin/tavg).
#   3) Use the VIF-consistent 4-predictor main SDM: BIO4, BIO5, BIO15, BIO18.
#   4) Keep spring temperature OUT of the main SDM and use March-May monthly Tmin
#      as a separate spring temperature-exposure assessment.
#   5) Run a sensitivity analysis with occurrence-density-weighted background
#      sampling versus the uniform-background main model.
#   6) Export clamping/extrapolation, area statistics, diagnostics, rasters and
#      reviewer-ready tables.
#
# IMPORTANT
#   - Do NOT overwrite/delete the old analysis. This script writes to a NEW folder:
#       REVISION_2026_EXPORTS
#   - Put this script in the project folder that contains occurrence.txt and data_env/.
#   - In RStudio: Session > Set Working Directory > To Source File Location
#     and then Source this entire script.
# ==============================================================================

rm(list = ls())
graphics.off()
cat("\014")

set.seed(123)
options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. PACKAGES AND USER SETTINGS
# ==============================================================================

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "sf", "terra", "geodata", "usdm",
  "biomod2", "ggplot2", "ggspatial", "scales", "viridis", "rnaturalearth",
  "rnaturalearthdata", "pROC", "spThin", "maxnet", "writexl", "reshape2"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, then rerun the script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Main settings
SEED <- 123
N_BACKGROUND <- 10000
N_DOMAIN_CORR <- 100000
SUIT_THRESHOLDS <- c(0.50, 0.60)
TEMP_THRESHOLDS <- c(-1.0, -1.5)

# Geographic extents (same as submitted analysis)
EU_XLIM <- c(-10, 40)
EU_YLIM <- c(34, 70)
ASIA_XLIM <- c(100, 160)
ASIA_YLIM <- c(20, 60)

extent_asia <- terra::ext(ASIA_XLIM[1], ASIA_XLIM[2], ASIA_YLIM[1], ASIA_YLIM[2])
extent_eu   <- terra::ext(EU_XLIM[1], EU_XLIM[2], 30, EU_YLIM[2])

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_root <- file.path(project_root, "REVISION_2026_EXPORTS")

for (d in c("MAIN", "SUPPLEMENTARY", "TABLES", "RASTERS", "BIOMOD_WORK", "TEMP")) {
  dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
}

# Force terra temporary files into the revision folder, not the system drive.
terra::terraOptions(
  tempdir = file.path(output_root, "TEMP"),
  progress = 1,
  memfrac = 0.60
)

# Locate occurrence file.
occurrence_candidates <- c(
  file.path(project_root, "occurrence.txt"),
  file.path(project_root, "occurrence(1).txt")
)
occurrence_file <- occurrence_candidates[file.exists(occurrence_candidates)][1]
if (is.na(occurrence_file) || length(occurrence_file) == 0) {
  stop("occurrence.txt was not found. Put it in the same project folder as this script.")
}

# WorldClim cache folder.
data_env_dir <- file.path(project_root, "data_env")
dir.create(data_env_dir, recursive = TRUE, showWarnings = FALSE)

biomod_root <- file.path(output_root, "BIOMOD_WORK")

# ==============================================================================
# 2. HELPERS
# ==============================================================================

save_csv <- function(x, filename) {
  readr::write_csv(x, file.path(output_root, "TABLES", filename))
}

save_xlsx <- function(x, filename) {
  writexl::write_xlsx(x, file.path(output_root, "TABLES", filename))
}

save_fig <- function(p, filename, folder = "MAIN", width = 9, height = 6.5, dpi = 450) {
  ggsave(
    filename = file.path(output_root, folder, filename),
    plot = p, width = width, height = height, units = "in", dpi = dpi, bg = "white"
  )
}

normalize_prediction <- function(r) {
  mx <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (is.finite(mx) && mx > 1) r <- r / 1000
  r
}

vif_from_df <- function(df) {
  cm <- cor(df, use = "complete.obs", method = "pearson")
  inv <- tryCatch(solve(cm), error = function(e) NULL)
  if (is.null(inv)) {
    return(tibble::tibble(Variable = colnames(cm), VIF = NA_real_))
  }
  tibble::tibble(Variable = colnames(cm), VIF = diag(inv))
}

extract_em_layer <- function(em_projection, ensemble_type = c("EMwmean", "EMmean")) {
  ensemble_type <- match.arg(ensemble_type)
  preds <- get_predictions(em_projection)
  nm <- names(preds)
  idx <- grep(ensemble_type, nm, fixed = TRUE)
  if (length(idx) == 0) stop("Could not find ", ensemble_type, " layer in ensemble predictions.")

  # Prefer the continuous merged ensemble layer.
  preferred <- idx[grepl("mergedData_mergedRun_mergedAlgo", nm[idx], fixed = TRUE)]
  if (length(preferred) > 0) idx <- preferred

  out <- preds[[idx[1]]]
  normalize_prediction(out)
}

sample_background <- function(sampling_raster, n, pres_xy, template,
                              method = c("random", "weights"), seed = 123) {
  method <- match.arg(method)
  set.seed(seed)
  pres_cells <- unique(terra::cellFromXY(template, as.matrix(pres_xy[, c("x", "y")])))
  pres_cells <- pres_cells[!is.na(pres_cells)]

  for (mult in c(1.25, 1.5, 2, 3)) {
    raw <- terra::spatSample(
      sampling_raster,
      size = ceiling(n * mult),
      method = method,
      xy = TRUE,
      na.rm = TRUE,
      as.df = TRUE
    )

    if (!all(c("x", "y") %in% names(raw))) {
      stop("terra::spatSample did not return x/y columns as expected.")
    }

    bg <- raw[, c("x", "y"), drop = FALSE]
    bg$cell <- terra::cellFromXY(template, as.matrix(bg[, c("x", "y")]))
    bg <- bg[!is.na(bg$cell), , drop = FALSE]
    bg <- bg[!bg$cell %in% pres_cells, , drop = FALSE]
    bg <- bg[!duplicated(bg$cell), , drop = FALSE]

    if (nrow(bg) >= n) {
      bg <- bg[seq_len(n), c("x", "y"), drop = FALSE]
      rownames(bg) <- NULL
      return(bg)
    }
  }

  stop("Could not obtain the requested number of unique background cells.")
}

fit_biomod_scenario <- function(resp_name, bg_points, bio_selected, env_eu,
                                var_import_n = 10, seed = 123,
                                build_ensemble = TRUE) {
  resp_xy <- rbind(pres_final, bg_points)
  resp_var <- c(rep(1, nrow(pres_final)), rep(0, nrow(bg_points)))

  env_vals <- terra::extract(
    bio_selected,
    terra::vect(resp_xy, geom = c("x", "y"), crs = "EPSG:4326")
  )[, -1, drop = FALSE]
  env_vals <- as.data.frame(env_vals)
  keep_complete <- complete.cases(env_vals)

  resp_xy2 <- resp_xy[keep_complete, , drop = FALSE]
  resp_var2 <- resp_var[keep_complete]

  cat("\n", resp_name, "\n", sep = "")
  cat("Presences entering BIOMOD2:", sum(resp_var2 == 1), "\n")
  cat("Background entering BIOMOD2:", sum(resp_var2 == 0), "\n")

  bm_data <- BIOMOD_FormatingData(
    resp.var = resp_var2,
    expl.var = bio_selected,
    resp.xy = resp_xy2,
    resp.name = resp_name,
    dir.name = biomod_root,
    PA.nb.rep = 0,
    filter.raster = TRUE
  )

  bm_model <- BIOMOD_Modeling(
    bm.format = bm_data,
    models = c("GLM", "GAM", "GBM", "RF", "MAXNET"),
    CV.strategy = "env",
    CV.k = 5,
    CV.balance = "presences",
    metric.eval = c("TSS", "ROC"),
    var.import = var_import_n,
    seed.val = seed
  )

  # ---------------------------------------------------------------------------
  # Projection of individual models
  # ---------------------------------------------------------------------------
  
  proj_name <- paste0("Europe_", resp_name)
  
  bm_proj <- BIOMOD_Projection(
    bm.mod = bm_model,
    new.env = env_eu,
    proj.name = proj_name,
    models.chosen = "all",
    metric.binary = "TSS",
    compress = "xz",
    build.clamping.mask = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Ensemble modeling and forecasting
  # Main analysis only.
  # For sensitivity runs where no model reaches the prespecified TSS threshold,
  # build_ensemble = FALSE prevents construction of an unsupported ensemble.
  # ---------------------------------------------------------------------------
  
  bm_em <- NULL
  bm_emproj <- NULL
  
  if (isTRUE(build_ensemble)) {
    
    bm_em <- BIOMOD_EnsembleModeling(
      bm.mod = bm_model,
      models.chosen = "all",
      em.by = "all",
      em.algo = c("EMmean", "EMwmean"),
      metric.select = "TSS",
      metric.select.thresh = 0.6,
      metric.eval = c("TSS", "ROC"),
      var.import = 3
    )
    
    bm_emproj <- tryCatch(
      {
        BIOMOD_EnsembleForecasting(
          bm.em = bm_em,
          bm.proj = bm_proj,
          models.chosen = "all"
        )
      },
      error = function(e) {
        
        warning(
          "BIOMOD_EnsembleForecasting failed: ",
          conditionMessage(e),
          "\nAttempting to recover the ensemble projection already written to disk."
        )
        
        emproj_files <- list.files(
          file.path(biomod_root, resp_name),
          pattern = "\\.ensemble\\.projection\\.out$",
          recursive = TRUE,
          full.names = TRUE
        )
        
        if (length(emproj_files) == 0) {
          stop(
            "Ensemble forecasting failed and no saved ensemble projection was found."
          )
        }
        
        emproj_info <- file.info(emproj_files)
        emproj_file <- rownames(emproj_info)[which.max(emproj_info$mtime)]
        
        loaded_object_name <- load(emproj_file)
        recovered <- get(loaded_object_name)
        
        if (!inherits(recovered, "BIOMOD.projection.out")) {
          stop(
            "Recovered ensemble projection is not a BIOMOD.projection.out object."
          )
        }
        
        message(
          "Recovered saved ensemble projection: ",
          emproj_file
        )
        
        recovered
      }
    )
  }

  list(
    data = bm_data,
    model = bm_model,
    ensemble = bm_em,
    projection = bm_proj,
    ensemble_projection = bm_emproj,
    resp_xy = resp_xy2,
    resp_var = resp_var2,
    proj_name = proj_name
  )
}

comparison_stats <- function(r1, r2, label1, label2) {
  if (!terra::compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- terra::resample(r2, r1, method = "bilinear")
  }
  d <- r2 - r1
  smp <- terra::spatSample(c(r1, r2), size = 100000, method = "random",
                           na.rm = TRUE, as.df = TRUE)
  rho <- if (ncol(smp) >= 2) cor(smp[[1]], smp[[2]], use = "complete.obs") else NA_real_

  tibble::tibble(
    Raster_1 = label1,
    Raster_2 = label2,
    Pearson_r_sample = rho,
    Mean_abs_difference = terra::global(abs(d), "mean", na.rm = TRUE)[1, 1],
    RMSE = sqrt(terra::global(d^2, "mean", na.rm = TRUE)[1, 1]),
    Mean_difference_R2_minus_R1 = terra::global(d, "mean", na.rm = TRUE)[1, 1]
  )
}

# ==============================================================================
# 3. OCCURRENCE CLEANING - TRANSPARENT COUNTS
# ==============================================================================

world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_vect <- terra::vect(world_sf)

gbif_data <- readr::read_tsv(occurrence_file, show_col_types = FALSE)
n_raw <- nrow(gbif_data)

points_coords <- gbif_data %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
n_coords <- nrow(points_coords)

points_unique_df <- points_coords %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
n_unique <- nrow(points_unique_df)

points_sf <- st_as_sf(
  points_unique_df,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

world_union <- st_union(world_sf)
land_index <- lengths(st_intersects(points_sf, world_union)) > 0
points_land <- points_sf[land_index, ]
points_land <- points_land[!duplicated(points_land$geometry), ]
n_land <- nrow(points_land)

study_region_sf <- st_crop(
  world_sf,
  xmin = ASIA_XLIM[1], xmax = ASIA_XLIM[2],
  ymin = ASIA_YLIM[1], ymax = ASIA_YLIM[2]
)

points_region <- points_land[
  st_within(points_land, st_union(study_region_sf), sparse = FALSE),
]
n_region <- nrow(points_region)

thin_data <- data.frame(
  species = "MiniKiwi",
  lat = st_coordinates(points_region)[, 2],
  long = st_coordinates(points_region)[, 1]
)

set.seed(SEED)
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

thin_sizes <- vapply(thinned, nrow, integer(1))
best_thin_index <- which.max(thin_sizes)
pres_thinned <- as.data.frame(thinned[[best_thin_index]])[, c("Longitude", "Latitude")]
colnames(pres_thinned) <- c("x", "y")
n_thinned <- nrow(pres_thinned)

thin_rep_summary <- tibble::tibble(
  replicate = seq_along(thin_sizes),
  n_retained = thin_sizes,
  selected = seq_along(thin_sizes) == best_thin_index
)
save_csv(thin_rep_summary, "Occurrence_thinning_replicates.csv")

# ==============================================================================
# 4. WORLDCLIM 2.1 AND CORRECT TEMPERATURE UNITS
# ==============================================================================

bio_global  <- geodata::worldclim_global(var = "bio",  res = 2.5, path = data_env_dir)
tmin_global <- geodata::worldclim_global(var = "tmin", res = 2.5, path = data_env_dir)
tavg_global <- geodata::worldclim_global(var = "tavg", res = 2.5, path = data_env_dir)

bio_asia <- terra::crop(bio_global, extent_asia)
names(bio_asia) <- paste0("BIO", 1:19)

# IMPORTANT REVISION CORRECTION:
# WorldClim 2.1 monthly temperature GeoTIFFs returned by geodata are already in °C.
# Therefore NO division by 10 is applied here.
spring_min <- terra::app(terra::crop(tmin_global[[3:4]], extent_asia), fun = min)
names(spring_min) <- "spring_min"

chilling_proxy <- terra::app(
  terra::crop(tavg_global[[c(1, 2, 3, 11, 12)]], extent_asia),
  fun = mean
)
names(chilling_proxy) <- "chilling_proxy"

candidate_stack <- c(bio_asia, spring_min, chilling_proxy)

# Raster-cell filtering is made explicit BEFORE BIOMOD2.
cell_template <- bio_asia[[1]]
thin_cells <- terra::cellFromXY(cell_template, as.matrix(pres_thinned[, c("x", "y")]))
keep_cell <- !is.na(thin_cells) & !duplicated(thin_cells)
pres_raster_filtered <- pres_thinned[keep_cell, , drop = FALSE]
n_raster_filtered <- nrow(pres_raster_filtered)

# Main revised model uses the VIF-consistent 4-variable set.
legacy5_keep <- c("spring_min", "BIO4", "BIO5", "BIO15", "BIO18")
final_keep <- c("BIO4", "BIO5", "BIO15", "BIO18")

legacy5_stack <- candidate_stack[[legacy5_keep]]
bio_selected <- candidate_stack[[final_keep]]

# Remove any presence lacking values for revised final predictors.
pres_env <- terra::extract(
  bio_selected,
  terra::vect(pres_raster_filtered, geom = c("x", "y"), crs = "EPSG:4326")
)[, -1, drop = FALSE]
keep_env <- complete.cases(as.data.frame(pres_env))
pres_final <- pres_raster_filtered[keep_env, , drop = FALSE]
n_final <- nrow(pres_final)
pres_vect <- terra::vect(pres_final, geom = c("x", "y"), crs = "EPSG:4326")

filtering_summary <- tibble::tibble(
  Step = c(
    "Raw GBIF records",
    "Records with non-missing coordinates",
    "Unique coordinate pairs",
    "Terrestrial records after land mask",
    "Records inside East Asian calibration region",
    "After 5-km thinning (largest deterministic replicate)",
    "After one-record-per-WorldClim-cell filtering",
    "Final presences with complete revised predictor values"
  ),
  N = c(n_raw, n_coords, n_unique, n_land, n_region, n_thinned, n_raster_filtered, n_final)
)
save_csv(filtering_summary, "Occurrence_filtering_summary.csv")

# Country distribution for final modelling presences.
pres_final_sf <- st_as_sf(pres_final, coords = c("x", "y"), crs = 4326, remove = FALSE)
pres_country <- st_join(
  pres_final_sf,
  world_sf[, c("admin", "iso_a2")],
  join = st_within,
  left = TRUE
) %>%
  st_drop_geometry() %>%
  mutate(
    Reviewer_region = case_when(
      admin == "China" ~ "China",
      admin == "Japan" ~ "Japan",
      grepl("Korea", admin, ignore.case = TRUE) ~ "Korea",
      TRUE ~ "Other"
    )
  ) %>%
  count(Reviewer_region, name = "N") %>%
  arrange(desc(N))
save_csv(pres_country, "Final_presence_distribution_China_Japan_Korea.csv")

# Also export raw GBIF countryCode distribution for transparency.
raw_country_codes <- gbif_data %>% count(countryCode, sort = TRUE, name = "N")
save_csv(raw_country_codes, "Raw_GBIF_countryCode_counts.csv")

# ==============================================================================
# 5. MULTICOLLINEARITY DIAGNOSTICS - OCCURRENCE AND DOMAIN
# ==============================================================================

candidate_vals <- terra::extract(candidate_stack, pres_vect)[, -1, drop = FALSE]
candidate_vals <- as.data.frame(candidate_vals)
candidate_vals <- candidate_vals[complete.cases(candidate_vals), , drop = FALSE]

cor_mat_candidates_occ <- cor(candidate_vals, use = "complete.obs", method = "pearson")
vif_candidates_step <- usdm::vifstep(candidate_vals, th = 10)

legacy5_vals <- candidate_vals[, legacy5_keep, drop = FALSE]
final4_vals <- candidate_vals[, final_keep, drop = FALSE]

cor_legacy5_occ <- cor(legacy5_vals, use = "complete.obs", method = "pearson")
cor_final4_occ <- cor(final4_vals, use = "complete.obs", method = "pearson")

vif_legacy5_explicit <- vif_from_df(legacy5_vals)
vif_final4_explicit <- vif_from_df(final4_vals)
vif_final4_step <- usdm::vifstep(final4_vals, th = 10)

# Domain-based correlations: random sample of complete land cells in calibration domain.
set.seed(SEED)
domain_sample5 <- terra::spatSample(
  legacy5_stack,
  size = N_DOMAIN_CORR,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)
domain_sample5 <- domain_sample5[complete.cases(domain_sample5), , drop = FALSE]
cor_legacy5_domain <- cor(domain_sample5[, legacy5_keep, drop = FALSE], use = "complete.obs")
cor_final4_domain <- cor(domain_sample5[, final_keep, drop = FALSE], use = "complete.obs")

save_csv(as.data.frame(cor_mat_candidates_occ) %>% tibble::rownames_to_column("Variable"),
         "Correlation_candidates_occurrence_points.csv")
save_csv(as.data.frame(cor_legacy5_occ) %>% tibble::rownames_to_column("Variable"),
         "Correlation_legacy5_occurrence_points.csv")
save_csv(as.data.frame(cor_final4_occ) %>% tibble::rownames_to_column("Variable"),
         "Correlation_revised4_occurrence_points.csv")
save_csv(as.data.frame(cor_legacy5_domain) %>% tibble::rownames_to_column("Variable"),
         "Correlation_legacy5_calibration_domain_sample.csv")
save_csv(as.data.frame(cor_final4_domain) %>% tibble::rownames_to_column("Variable"),
         "Correlation_revised4_calibration_domain_sample.csv")
save_csv(vif_legacy5_explicit, "VIF_legacy5_explicit_before_removal.csv")
save_csv(vif_final4_explicit, "VIF_revised4_explicit.csv")
save_csv(tibble::as_tibble(vif_final4_step@results), "VIF_revised4_vifstep.csv")
save_csv(tibble::as_tibble(vif_candidates_step@results), "VIF_candidate_predictors_vifstep.csv")

cat("\nLegacy 5-variable pairwise correlation at presences (spring_min vs BIO4): ",
    round(cor_legacy5_occ["spring_min", "BIO4"], 3), "\n", sep = "")
cat("Legacy 5-variable VIF for spring_min: ",
    round(vif_legacy5_explicit$VIF[vif_legacy5_explicit$Variable == "spring_min"], 3), "\n", sep = "")

# ==============================================================================
# 6. BACKGROUND SAMPLING
# ==============================================================================

asia_land_vect <- terra::crop(world_vect, extent_asia)
asia_land_mask <- terra::rasterize(
  asia_land_vect,
  bio_selected[[1]],
  field = 1,
  background = NA
)

# Uniform background = main revised model (same baseline strategy as submitted paper).
bg_uniform <- sample_background(
  sampling_raster = asia_land_mask,
  n = N_BACKGROUND,
  pres_xy = pres_final,
  template = bio_selected[[1]],
  method = "random",
  seed = SEED
)
save_csv(bg_uniform, "Background_uniform_10000.csv")

# Occurrence-density surface for quantitative background sensitivity analysis.
# Uses CLEANED records before thinning, aligned to the modelling raster grid.
points_region_vect <- terra::vect(points_region)
bias_count <- terra::rasterize(
  points_region_vect,
  bio_selected[[1]],
  field = 1,
  fun = "sum",
  background = 0
)

# ~1 degree moving window at 2.5 arc-min resolution: 25 x 25 cells.
bias_smooth <- terra::focal(
  bias_count,
  w = matrix(1, 25, 25),
  fun = "sum",
  na.policy = "omit",
  fillvalue = 0
)

# Small positive floor preserves the full accessible land domain while strongly
# weighting cells near higher occurrence-sampling density.
bias_weights <- terra::ifel(!is.na(asia_land_mask), bias_smooth + 0.01, NA)
names(bias_weights) <- "OccurrenceDensityWeight"

terra::writeRaster(
  bias_weights,
  file.path(output_root, "RASTERS", "background_occurrence_density_weights.tif"),
  overwrite = TRUE
)

bg_bias <- sample_background(
  sampling_raster = bias_weights,
  n = N_BACKGROUND,
  pres_xy = pres_final,
  template = bio_selected[[1]],
  method = "weights",
  seed = SEED
)
save_csv(bg_bias, "Background_occurrence_density_weighted_10000.csv")

# ==============================================================================
# 7. EUROPEAN PREDICTORS + SPRING TMIN EXPOSURE LAYERS
# ==============================================================================

bio_eu_base <- terra::crop(bio_global, extent_eu)
names(bio_eu_base) <- paste0("BIO", 1:19)
env_eu <- bio_eu_base[[final_keep]]

if (!identical(names(bio_selected), names(env_eu))) {
  stop("Training and European projection predictor names do not match.")
}

# IMPORTANT: WorldClim 2.1 monthly tmin values returned by geodata
# are already in °C. NO division by 10 is applied.
tmin_spring_eu <- terra::crop(tmin_global[[3:5]], extent_eu)
names(tmin_spring_eu) <- c("MarchTmin_C", "AprilTmin_C", "MayTmin_C")

tmin_march_eu <- tmin_spring_eu[["MarchTmin_C"]]
tmin_april_eu <- tmin_spring_eu[["AprilTmin_C"]]
tmin_may_eu   <- tmin_spring_eu[["MayTmin_C"]]
names(tmin_march_eu) <- "MarchTmin_C"

# ==============================================================================
# 8. MAIN REVISED SDM - FOUR PREDICTORS + UNIFORM BACKGROUND
# ==============================================================================

main_run <- fit_biomod_scenario(
  resp_name = "MiniKiwi_REV4_UNIFORM",
  bg_points = bg_uniform,
  bio_selected = bio_selected,
  env_eu = env_eu,
  var_import_n = 10,
  seed = SEED
)

suit_main <- extract_em_layer(main_run$ensemble_projection, "EMwmean")
suit_mean <- extract_em_layer(main_run$ensemble_projection, "EMmean")

terra::writeRaster(suit_main,
                   file.path(output_root, "RASTERS", "suitability_REV4_uniform_EMwmean.tif"),
                   overwrite = TRUE)
terra::writeRaster(suit_mean,
                   file.path(output_root, "RASTERS", "suitability_REV4_uniform_EMmean.tif"),
                   overwrite = TRUE)
terra::writeRaster(tmin_march_eu,
                   file.path(output_root, "RASTERS", "march_mean_minimum_temperature_C.tif"),
                   overwrite = TRUE)

terra::writeRaster(tmin_april_eu,
                   file.path(output_root, "RASTERS", "april_mean_minimum_temperature_C.tif"),
                   overwrite = TRUE)

terra::writeRaster(tmin_may_eu,
                   file.path(output_root, "RASTERS", "may_mean_minimum_temperature_C.tif"),
                   overwrite = TRUE)

# Main model performance and variable importance.
evals_main <- get_evaluations(main_run$model)

tss_main <- evals_main %>%
  filter(metric.eval == "TSS") %>%
  select(algo, run, calibration, validation) %>%
  as_tibble()

tss_summary_main <- tss_main %>%
  group_by(algo) %>%
  summarise(Mean_TSS = mean(validation, na.rm = TRUE),
            SD_TSS = sd(validation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Mean_TSS))

roc_main <- evals_main %>%
  filter(metric.eval == "AUCroc") %>%
  select(algo, run, calibration, validation) %>%
  as_tibble()

roc_summary_main <- roc_main %>%
  group_by(algo) %>%
  summarise(Mean_ROC = mean(validation, na.rm = TRUE),
            SD_ROC = sd(validation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Mean_ROC))

var_imp_main <- as.data.frame(get_variables_importance(main_run$model)) %>%
  group_by(expl.var) %>%
  summarise(MeanImportance = mean(var.imp, na.rm = TRUE),
            SDImportance = sd(var.imp, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MeanImportance))

save_csv(tss_main, "TSS_per_run_REV4_uniform.csv")
save_csv(tss_summary_main, "TSS_summary_REV4_uniform.csv")
save_csv(roc_main, "ROC_per_run_REV4_uniform.csv")
save_csv(roc_summary_main, "ROC_summary_REV4_uniform.csv")
save_csv(var_imp_main, "Variable_importance_REV4_uniform.csv")

# Cross-validation design summary (answers reviewer reproducibility question).
calib_lines <- as.data.frame(get_calib_lines(main_run$model))
cv_summary <- lapply(names(calib_lines), function(nm) {
  x <- calib_lines[[nm]]
  tibble::tibble(
    Fold = nm,
    N_calibration = sum(x %in% TRUE, na.rm = TRUE),
    N_validation = sum(x %in% FALSE, na.rm = TRUE),
    N_missing = sum(is.na(x))
  )
}) %>% bind_rows()
save_csv(cv_summary, "Environmental_CV_fold_summary_REV4_uniform.csv")

# ==============================================================================
# 9. MODEL DISAGREEMENT + CLAMPING/EXTRAPOLATION
# ==============================================================================

model_stack_main <- get_predictions(main_run$projection)
model_stack_main <- normalize_prediction(model_stack_main)
sd_map_main <- terra::app(model_stack_main, fun = sd, na.rm = TRUE)
names(sd_map_main) <- "InterModel_SD"
terra::writeRaster(sd_map_main,
                   file.path(output_root, "RASTERS", "intermodel_disagreement_SD_REV4_uniform.tif"),
                   overwrite = TRUE)

# BIOMOD2 saves the clamping mask inside the projection folder when
# build.clamping.mask = TRUE. Find it robustly rather than assuming a fixed name.
clamp_files <- list.files(
  biomod_root,
  pattern = "clamp.*\\.(tif|grd)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Prefer a clamping file belonging to the main revised scenario.
clamp_main_candidates <- clamp_files[grepl("MiniKiwi_REV4_UNIFORM", clamp_files, fixed = TRUE)]
if (length(clamp_main_candidates) == 0) clamp_main_candidates <- clamp_files

clamp_main <- NULL
if (length(clamp_main_candidates) > 0) {
  clamp_main <- terra::rast(clamp_main_candidates[1])
  terra::writeRaster(
    clamp_main,
    file.path(output_root, "RASTERS", "clamping_mask_REV4_uniform.tif"),
    overwrite = TRUE
  )
  writeLines(clamp_main_candidates[1], file.path(output_root, "TABLES", "Clamping_source_file.txt"))
} else {
  warning("No clamping-mask file was found automatically. BIOMOD2 projection still ran with build.clamping.mask=TRUE.")
}

# ==============================================================================
# 10. EUROPE MASK + SPRING TEMPERATURE-EXPOSURE / AREA STATISTICS
# ==============================================================================

europe <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe <- europe[europe$continent == "Europe", ]
europe_vect <- terra::vect(europe)

suit_eur <- terra::mask(terra::crop(suit_main, europe_vect), europe_vect)

# Align March, April and May Tmin layers with the European suitability raster.
align_spring_layer <- function(r, template, mask_vect) {
  x <- terra::crop(r, template)
  if (!terra::compareGeom(x, template, stopOnError = FALSE)) {
    x <- terra::resample(x, template, method = "bilinear")
  }
  terra::mask(x, mask_vect)
}

march_eur <- align_spring_layer(tmin_march_eu, suit_eur, europe_vect)
april_eur <- align_spring_layer(tmin_april_eu, suit_eur, europe_vect)
may_eur   <- align_spring_layer(tmin_may_eu, suit_eur, europe_vect)

spring_tmin_eur <- c(march_eur, april_eur, may_eur)
names(spring_tmin_eur) <- c("March", "April", "May")

cell_area_km2 <- terra::cellSize(suit_eur, unit = "km")

area_where <- function(condition_raster) {
  x <- cell_area_km2 * condition_raster
  as.numeric(terra::global(x, "sum", na.rm = TRUE)[1, 1])
}

# --------------------------------------------------------------------------
# Area statistics by month, suitability threshold and temperature threshold.
# Temperature thresholds are descriptive exposure scenarios, NOT estimates
# of frost-damage probability.
# --------------------------------------------------------------------------

area_rows <- list()
k <- 1

for (month_name in names(spring_tmin_eur)) {
  
  temp_layer <- spring_tmin_eur[[month_name]]
  
  for (s_thr in SUIT_THRESHOLDS) {
    
    suitable <- suit_eur >= s_thr
    suitable_area <- area_where(suitable)
    
    for (t_thr in TEMP_THRESHOLDS) {
      
      exposed <- temp_layer <= t_thr
      exposed_suitable_area <- area_where(suitable & exposed)
      outside_exposure_area <- area_where(suitable & !exposed)
      
      area_rows[[k]] <- tibble::tibble(
        Month = month_name,
        Suitability_threshold = s_thr,
        Tmin_exposure_threshold_C = t_thr,
        Suitable_area_km2 = suitable_area,
        Suitable_area_exposed_km2 = exposed_suitable_area,
        Suitable_area_outside_exposure_km2 = outside_exposure_area,
        Percent_of_suitable_area_exposed =
          100 * exposed_suitable_area / suitable_area,
        Percent_of_suitable_area_outside_exposure =
          100 * outside_exposure_area / suitable_area
      )
      
      k <- k + 1
    }
  }
}

area_stats <- dplyr::bind_rows(area_rows)

save_csv(
  area_stats,
  "European_suitability_spring_temperature_exposure_area_statistics.csv"
)

# --------------------------------------------------------------------------
# Descriptive Tmin exposure classes for March, April and May.
# These are climatological exposure classes, not frost-risk classes.
# --------------------------------------------------------------------------

class_defs <- tibble::tibble(
  ExposureClass = 1:5,
  Label = c(
    "≤ -1.5°C",
    "> -1.5 to -1.0°C",
    "> -1.0 to 0°C",
    "> 0 to 1°C",
    "> 1°C"
  )
)

class_area_list <- list()

for (month_name in names(spring_tmin_eur)) {
  
  temp_layer <- spring_tmin_eur[[month_name]]
  
  exposure_class <- terra::ifel(
    temp_layer <= -1.5, 1,
    terra::ifel(
      temp_layer <= -1.0, 2,
      terra::ifel(
        temp_layer <= 0, 3,
        terra::ifel(temp_layer <= 1, 4, 5)
      )
    )
  )
  
  names(exposure_class) <- paste0(month_name, "_ExposureClass")
  
  terra::writeRaster(
    exposure_class,
    file.path(
      output_root,
      "RASTERS",
      paste0(month_name, "_Tmin_exposure_classes.tif")
    ),
    overwrite = TRUE
  )
  
  month_class_area <- lapply(1:5, function(cl) {
    tibble::tibble(
      Month = month_name,
      ExposureClass = cl,
      Area_km2 = area_where(exposure_class == cl)
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(class_defs, by = "ExposureClass") %>%
    dplyr::mutate(
      Percent_Europe_modelled_area =
        100 * Area_km2 / sum(Area_km2, na.rm = TRUE)
    )
  
  class_area_list[[month_name]] <- month_class_area
}

class_area <- dplyr::bind_rows(class_area_list)

save_csv(
  class_area,
  "Spring_Tmin_exposure_class_area_Europe.csv"
)
# ==============================================================================
# 11. BIAS-SENSITIVITY RUN - SAME 4 PREDICTORS, DENSITY-WEIGHTED BACKGROUND
# Bias-weighted background sensitivity analysis.
# Models are fitted from scratch using the same four predictors and the
# occurrence-density-weighted background generated in Section 6.
# The prespecified TSS >= 0.60 ensemble-selection criterion is retained.
# No ensemble is constructed unless supported by this criterion.
# ==============================================================================

set.seed(SEED)

bias_run <- fit_biomod_scenario(
  resp_name = "MiniKiwi_REV4_BIASWEIGHT",
  bg_points = bg_bias,
  bio_selected = bio_selected,
  env_eu = env_eu,
  var_import_n = 10,
  seed = SEED,
  build_ensemble = FALSE
)

# ------------------------------------------------------------------------------
# Model evaluation
# ------------------------------------------------------------------------------

evals_bias <- get_evaluations(
  bias_run$model
)

tss_bias_summary <- evals_bias %>%
  filter(metric.eval == "TSS") %>%
  group_by(algo) %>%
  summarise(
    Mean_TSS = mean(validation, na.rm = TRUE),
    SD_TSS = sd(validation, na.rm = TRUE),
    Max_TSS = max(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_TSS))

roc_bias_summary <- evals_bias %>%
  filter(metric.eval == "AUCroc") %>%
  group_by(algo) %>%
  summarise(
    Mean_AUC = mean(validation, na.rm = TRUE),
    SD_AUC = sd(validation, na.rm = TRUE),
    Max_AUC = max(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_AUC))

save_csv(
  tss_bias_summary,
  "TSS_summary_REV4_biasweighted.csv"
)

save_csv(
  roc_bias_summary,
  "ROC_summary_REV4_biasweighted.csv"
)

# ------------------------------------------------------------------------------
# Prespecified ensemble-selection criterion
# ------------------------------------------------------------------------------

max_bias_tss <- max(
  evals_bias$validation[
    evals_bias$metric.eval == "TSS"
  ],
  na.rm = TRUE
)

bias_ensemble_supported <- is.finite(max_bias_tss) &&
  max_bias_tss >= 0.60

bias_threshold_result <- tibble::tibble(
  Ensemble_selection_metric = "TSS",
  Prespecified_threshold = 0.60,
  Maximum_validation_TSS = max_bias_tss,
  Ensemble_constructed = FALSE,
  Interpretation = if (
    !bias_ensemble_supported
  ) {
    paste(
      "No bias-weighted model met the prespecified TSS >= 0.60",
      "ensemble-selection criterion; therefore no bias-weighted ensemble",
      "projection or threshold-area comparison was produced."
    )
  } else {
    paste(
      "At least one bias-weighted model met the prespecified TSS >= 0.60",
      "criterion. Ensemble construction was not performed automatically",
      "because this sensitivity analysis was prespecified as a diagnostic run."
    )
  }
)

save_csv(
  bias_threshold_result,
  "Bias_sensitivity_ensemble_threshold_result.csv"
)

print(tss_bias_summary)
print(roc_bias_summary)
print(bias_threshold_result)

cat(
  "\nBias-weighted sensitivity analysis completed from scratch.\n",
  "Maximum validation TSS: ",
  round(max_bias_tss, 3),
  "\n",
  "TSS >= 0.60 criterion met: ",
  bias_ensemble_supported,
  "\n",
  sep = ""
)
# ==============================================================================
# 12. OPTIONAL COMPARISON WITH THE SUBMITTED 5-PREDICTOR RASTER
# ==============================================================================

legacy_candidates <- c(
  file.path(project_root, "FINAL_EXPORTS", "RASTERS", "suitability_EMwmean.tif"),
  file.path(project_root, "suitability_EMwmean.tif")
)
legacy_file <- legacy_candidates[file.exists(legacy_candidates)][1]

if (!is.na(legacy_file) && length(legacy_file) > 0) {
  legacy_suit <- normalize_prediction(terra::rast(legacy_file))
  legacy_compare <- comparison_stats(
    legacy_suit, suit_main,
    "Submitted 5-predictor model", "Revised 4-predictor model"
  )
  save_csv(legacy_compare, "Sensitivity_submitted5_vs_revised4_summary.csv")

  if (!terra::compareGeom(legacy_suit, suit_main, stopOnError = FALSE)) {
    legacy_suit <- terra::resample(legacy_suit, suit_main, method = "bilinear")
  }
  legacy_diff <- suit_main - legacy_suit
  names(legacy_diff) <- "Revised4MinusSubmitted5"
  terra::writeRaster(
    legacy_diff,
    file.path(output_root, "RASTERS", "sensitivity_difference_revised4_minus_submitted5.tif"),
    overwrite = TRUE
  )
} else {
  writeLines(
    "Legacy submitted raster was not found automatically. Copy FINAL_EXPORTS/RASTERS/suitability_EMwmean.tif into the project if you want automatic old-vs-new comparison.",
    file.path(output_root, "TABLES", "Legacy_raster_comparison_NOT_RUN.txt")
  )
}

# ==============================================================================
# 13. PUBLICATION FIGURES
# ==============================================================================

# ------------------------------------------------------------------------------
# COMMON PUBLICATION THEME
# ------------------------------------------------------------------------------

theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = 13,
        face = "plain",
        hjust = 0,
        margin = margin(b = 10)
      ),
      plot.caption = element_text(
        size = 10.5,
        face = "plain",
        hjust = 0,
        margin = margin(t = 8)
      ),
      axis.title = element_text(
        size = 13,
        face = "bold"
      ),
      axis.text = element_text(
        size = 11
      ),
      legend.title = element_text(
        size = 12,
        face = "bold"
      ),
      legend.text = element_text(
        size = 11
      ),
      strip.text = element_text(
        size = 13,
        face = "bold"
      ),
      strip.background = element_rect(
        fill = "white",
        color = "grey70",
        linewidth = 0.3
      ),
      panel.grid = element_blank(),
      plot.margin = margin(12, 15, 12, 12)
    )
}


# ==============================================================================
# FIGURE 01 — OCCURRENCE RECORDS
# ==============================================================================

pres_plot <- st_as_sf(
  pres_final,
  coords = c("x", "y"),
  crs = 4326
)

p_occ <- ggplot() +
  geom_sf(
    data = study_region_sf,
    fill = "grey96",
    color = "grey70",
    linewidth = 0.25
  ) +
  geom_sf(
    data = points_region,
    color = "grey55",
    size = 0.55,
    alpha = 0.40
  ) +
  geom_sf(
    data = pres_plot,
    color = "#B2180B",
    size = 0.75,
    alpha = 0.85
  ) +
  coord_sf(
    xlim = ASIA_XLIM,
    ylim = ASIA_YLIM,
    expand = FALSE
  ) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.20, "in"),
    pad_y = grid::unit(0.20, "in"),
    height = grid::unit(0.55, "in"),
    width = grid::unit(0.45, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  labs(
    title = "Occurrence records used for model calibration",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_publication()

save_fig(
  p_occ,
  "Figure_01_Occurrence_records.png",
  width = 8.5,
  height = 6.5
)


# ==============================================================================
# FIGURE 02 — CANDIDATE-PREDICTOR CORRELATION HEATMAP
# BIO1–BIO19 ONLY
# ==============================================================================

# Candidate SDM predictors presented in the manuscript are the
# 19 WorldClim bioclimatic variables.
bio_candidate_names <- paste0("BIO", 1:19)

candidate_bio_vals <- candidate_vals[
  ,
  bio_candidate_names,
  drop = FALSE
]

cor_bio_candidates <- cor(
  candidate_bio_vals,
  use = "complete.obs",
  method = "pearson"
)
# ------------------------------------------------------------------------------
# Statistical significance of Pearson correlations — Figure 2
# ------------------------------------------------------------------------------

n_cor <- nrow(candidate_bio_vals)

# Two-sided p-values corresponding to each Pearson correlation coefficient
t_cor <- cor_bio_candidates * sqrt(
  (n_cor - 2) / (1 - cor_bio_candidates^2)
)

p_bio_candidates <- 2 * pt(
  -abs(t_cor),
  df = n_cor - 2
)

# Diagonal correlations (r = 1) are not tested
diag(p_bio_candidates) <- NA

# Export p-value matrix
save_csv(
  as.data.frame(p_bio_candidates) %>%
    tibble::rownames_to_column("Variable"),
  "Pvalues_candidate_predictor_correlations.csv"
)

# Check significance of the strong correlations displayed in Figure 2
strong_cor_mask <- abs(cor_bio_candidates) >= 0.75 &
  abs(cor_bio_candidates) < 0.999

strong_p_values <- p_bio_candidates[strong_cor_mask]

cat(
  "\nFigure 2 strong correlations (|r| >= 0.75):\n",
  "Number of tested strong correlations: ", length(strong_p_values), "\n",
  "Maximum p-value: ", max(strong_p_values, na.rm = TRUE), "\n",
  "All p < 0.001: ", all(strong_p_values < 0.001, na.rm = TRUE), "\n",
  sep = ""
)
# Keep lower triangle + diagonal.
cor_candidates_plot <- cor_bio_candidates
cor_candidates_plot[upper.tri(cor_candidates_plot)] <- NA

cor_candidates_long <- reshape2::melt(
  cor_candidates_plot,
  na.rm = TRUE
)

colnames(cor_candidates_long) <- c(
  "Var1",
  "Var2",
  "Correlation"
)

candidate_order <- colnames(cor_bio_candidates)

cor_candidates_long$Var1 <- factor(
  cor_candidates_long$Var1,
  levels = candidate_order
)

cor_candidates_long$Var2 <- factor(
  cor_candidates_long$Var2,
  levels = rev(candidate_order)
)

# Annotate strong pairwise correlations only.
cor_candidates_long$Label <- ifelse(
  abs(cor_candidates_long$Correlation) >= 0.75 &
    abs(cor_candidates_long$Correlation) < 0.999,
  sprintf("%.2f", cor_candidates_long$Correlation),
  ""
)

p_cor_candidates <- ggplot(
  cor_candidates_long,
  aes(
    x = Var1,
    y = Var2,
    fill = Correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  geom_text(
    aes(label = Label),
    size = 3.4,
    family = "Arial",
    fontface = "plain",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#E3362A",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Pearson r"
  ) +
  coord_fixed() +
  labs(
    title = "Correlation heatmap of candidate predictors",
    x = NULL,
    y = NULL
  ) +
  theme_publication(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10.5
    ),
    axis.text.y = element_text(
      size = 10.5
    ),
    axis.title = element_blank(),
    legend.position = "right"
  )

save_fig(
  p_cor_candidates,
  "Figure_02_Candidate_predictor_correlations.png",
  width = 9.5,
  height = 8.5
)


# ==============================================================================
# FIGURE 03 — FINAL FOUR-PREDICTOR CORRELATION MATRIX
# ==============================================================================

cor4 <- cor_final4_occ
cor4[upper.tri(cor4)] <- NA

cor4_long <- reshape2::melt(
  cor4,
  na.rm = TRUE
)

colnames(cor4_long) <- c(
  "Var1",
  "Var2",
  "Correlation"
)

final_order <- colnames(cor_final4_occ)

cor4_long$Var1 <- factor(
  cor4_long$Var1,
  levels = final_order
)

cor4_long$Var2 <- factor(
  cor4_long$Var2,
  levels = rev(final_order)
)

p_cor <- ggplot(
  cor4_long,
  aes(
    x = Var1,
    y = Var2,
    fill = Correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Correlation)),
    size = 4.3,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#C95B5B",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Pearson r"
  ) +
  coord_fixed() +
  labs(
    title = "Correlation matrix of the final predictor set",
    subtitle = "Reduced multicollinearity after predictor selection",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 13,
      face = "plain",
      hjust = 0,
      margin = margin(b = 8)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11
    ),
    axis.text.y = element_text(
      size = 11
    ),
    axis.title = element_blank(),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      size = 11
    ),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

save_fig(
  p_cor,
  "Figure_03_Final_predictor_correlations.png",
  width = 8.5,
  height = 8.0
)

# ==============================================================================
# PREPARE EUROPEAN MAP DATA
# ==============================================================================

suit_eur_df <- as.data.frame(
  suit_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(suit_eur_df) <- c(
  "x",
  "y",
  "Suitability"
)

spring_df <- as.data.frame(
  spring_tmin_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(spring_df) <- c(
  "x",
  "y",
  "March",
  "April",
  "May"
)

spring_long <- tidyr::pivot_longer(
  spring_df,
  cols = c("March", "April", "May"),
  names_to = "Month",
  values_to = "Tmin_C"
)

spring_long$Month <- factor(
  spring_long$Month,
  levels = c("March", "April", "May")
)

# ==============================================================================
# FIGURE 04 — GLM ROC CURVE FROM ENVIRONMENTAL CROSS-VALIDATION RUNS
# ==============================================================================

# Observed response used by BIOMOD2
glm_obs <- biomod2::get_formal_data(
  main_run$model,
  subinfo = "resp.var"
)

glm_obs <- as.numeric(glm_obs)

# Stored GLM predictions for the 20 environmental cross-validation runs.
# The allRun model is excluded because it has no independent validation subset.
glm_pred <- biomod2::get_predictions(
  main_run$model,
  algo = "GLM",
  run = paste0("RUN", 1:20),
  model.as.col = TRUE
)

glm_pred <- as.data.frame(glm_pred)

# Environmental cross-validation membership
glm_cv <- as.data.frame(
  biomod2::get_calib_lines(main_run$model)
)

# Keep only the 20 environmental CV runs
glm_cv <- glm_cv[
  ,
  paste0("_allData_RUN", 1:20),
  drop = FALSE
]

# Safety checks
if (nrow(glm_pred) != length(glm_obs)) {
  stop(
    "Prediction rows (", nrow(glm_pred),
    ") do not match response length (",
    length(glm_obs), ")."
  )
}

if (nrow(glm_cv) != length(glm_obs)) {
  stop(
    "CV rows (", nrow(glm_cv),
    ") do not match response length (",
    length(glm_obs), ")."
  )
}

if (ncol(glm_pred) != 20) {
  stop(
    "Expected 20 GLM ENV prediction columns but found ",
    ncol(glm_pred), "."
  )
}

# Common false-positive-rate grid
fpr_grid <- seq(
  0,
  1,
  length.out = 501
)

roc_list <- vector(
  "list",
  20
)

auc_values <- numeric(20)

# Calculate ROC separately for each validation run
for (i in 1:20) {
  
  # FALSE identifies observations withheld for validation
  validation_index <- which(
    glm_cv[[i]] %in% FALSE
  )
  
  obs_i <- glm_obs[
    validation_index
  ]
  
  pred_i <- glm_pred[
    validation_index,
    i
  ]
  
  keep_i <- is.finite(obs_i) &
    is.finite(pred_i)
  
  obs_i <- obs_i[
    keep_i
  ]
  
  pred_i <- pred_i[
    keep_i
  ]
  
  if (length(unique(obs_i)) < 2) {
    stop(
      "RUN", i,
      " does not contain both response classes."
    )
  }
  
  roc_i <- pROC::roc(
    response = obs_i,
    predictor = pred_i,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_values[i] <- as.numeric(
    pROC::auc(roc_i)
  )
  
  roc_tmp <- data.frame(
    FPR = 1 - roc_i$specificities,
    TPR = roc_i$sensitivities
  )
  
  # Collapse duplicate FPR values
  roc_tmp <- roc_tmp %>%
    dplyr::group_by(FPR) %>%
    dplyr::summarise(
      TPR = max(
        TPR,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(FPR)
  
  # Guarantee ROC endpoints
  roc_tmp <- dplyr::bind_rows(
    data.frame(
      FPR = 0,
      TPR = 0
    ),
    roc_tmp,
    data.frame(
      FPR = 1,
      TPR = 1
    )
  ) %>%
    dplyr::group_by(FPR) %>%
    dplyr::summarise(
      TPR = max(
        TPR,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(FPR)
  
  # Interpolate each ROC curve to the same FPR grid
  interp_tpr <- approx(
    x = roc_tmp$FPR,
    y = roc_tmp$TPR,
    xout = fpr_grid,
    method = "linear",
    ties = max,
    rule = 2
  )$y
  
  roc_list[[i]] <- data.frame(
    Run = paste0(
      "RUN",
      i
    ),
    FPR = fpr_grid,
    TPR = interp_tpr
  )
}

# Combine the 20 validation ROC curves
roc_all <- dplyr::bind_rows(
  roc_list
)

# Mean ROC curve and pointwise 2.5th–97.5th percentile envelope
# This follows the same interpretation used for the submitted Figure 4.
roc_summary_glm <- roc_all %>%
  dplyr::group_by(FPR) %>%
  dplyr::summarise(
    Mean_TPR = mean(
      TPR,
      na.rm = TRUE
    ),
    Lower_TPR = stats::quantile(
      TPR,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Upper_TPR = stats::quantile(
      TPR,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

# BIOMOD2 validation AUC values for GLM
glm_auc_biomod <- roc_main %>%
  dplyr::filter(
    algo == "GLM",
    is.finite(validation)
  )

glm_auc_mean <- mean(
  glm_auc_biomod$validation,
  na.rm = TRUE
)

glm_auc_sd <- sd(
  glm_auc_biomod$validation,
  na.rm = TRUE
)

# Independent check against AUC calculated directly from predictions
cat(
  "\nGLM AUC from BIOMOD2:",
  sprintf(
    "%.3f ± %.3f",
    glm_auc_mean,
    glm_auc_sd
  ),
  "\n"
)

cat(
  "GLM AUC recalculated from validation predictions:",
  sprintf(
    "%.3f ± %.3f",
    mean(auc_values),
    sd(auc_values)
  ),
  "\n"
)

# ------------------------------------------------------------------------------
# FINAL PUBLICATION FIGURE
# ------------------------------------------------------------------------------

p_roc <- ggplot(
  roc_summary_glm,
  aes(
    x = FPR,
    y = Mean_TPR
  )
) +
  
  # Pointwise 2.5th–97.5th percentile envelope across validation runs
  geom_ribbon(
    aes(
      ymin = Lower_TPR,
      ymax = Upper_TPR
    ),
    fill = "#3182BD",
    alpha = 0.25
  ) +
  
  # Random-classification diagonal
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey35"
  ) +
  
  # Mean ROC curve
  geom_line(
    color = "#B30000",
    linewidth = 1.20
  ) +
  
  # Mean validation AUC ± SD
  annotate(
    "text",
    x = 0.64,
    y = 0.14,
    label = sprintf(
      "Mean validation AUC = %.3f ± %.3f",
      glm_auc_mean,
      glm_auc_sd
    ),
    size = 4.5
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    labels = sprintf(
      "%.2f",
      c(
        0,
        0.25,
        0.50,
        0.75,
        1
      )
    ),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    labels = sprintf(
      "%.2f",
      c(
        0,
        0.25,
        0.50,
        0.75,
        1
      )
    ),
    expand = c(0, 0)
  ) +
  
  coord_equal() +
  
  labs(
    title = "ROC curve for the GLM model",
    subtitle = "Environmental cross-validation runs",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 13,
      face = "plain",
      hjust = 0,
      margin = margin(b = 10)
    ),
    axis.title = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text = element_text(
      size = 11
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.30
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      12,
      15,
      12,
      12
    )
  )

# Save Figure 4.
# save_fig() already saves to REVISION_2026_EXPORTS/MAIN by default.
save_fig(
  p_roc,
  "Figure_04_ROC_curve_GLM.png",
  folder = "MAIN",
  width = 7.5,
  height = 7.5,
  dpi = 450
)

print(p_roc)
# Confirm exact output location
figure4_file <- file.path(
  output_root,
  "MAIN",
  "Figure_04_ROC_curve_GLM.png"
)

cat(
  "\nFINAL Figure 4 saved here:\n",
  normalizePath(
    figure4_file,
    mustWork = TRUE
  ),
  "\n"
)
# ==============================================================================
# FIGURE 05 — MODEL PERFORMANCE (TSS)
# ==============================================================================

# Use validation TSS values from the final uniform-background SDM.
# Exclude allRun because Figure 5 summarizes environmental CV validation runs.

tss_fig5 <- tss_main %>%
  dplyr::filter(
    is.finite(validation),
    run != "allRun"
  ) %>%
  dplyr::transmute(
    Model = algo,
    Run   = run,
    TSS   = validation
  )

# Fixed algorithm order
tss_fig5$Model <- factor(
  tss_fig5$Model,
  levels = c(
    "GLM",
    "MAXNET",
    "RF",
    "GBM",
    "GAM"
  )
)

# Safety check:
# 5 algorithms × 20 ENV validation runs = 100 observations
if (nrow(tss_fig5) != 100) {
  stop(
    "Figure 5 expected 100 validation TSS values ",
    "(5 algorithms × 20 ENV runs), but found ",
    nrow(tss_fig5), "."
  )
}

# Mean TSS for each algorithm
tss_fig5_mean <- tss_fig5 %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    Mean_TSS = mean(TSS, na.rm = TRUE),
    .groups = "drop"
  )

# Print means as a numerical check
print(tss_fig5_mean)

# ==============================================================================
# FIGURE 05 — MODEL PERFORMANCE BASED ON TSS
# ==============================================================================

# Match original figure order
tss_fig5$Model <- factor(
  tss_fig5$Model,
  levels = c("RF", "MAXNET", "GLM", "GBM", "GAM")
)

p_tss <- ggplot(
  tss_fig5,
  aes(x = TSS, y = Model)
) +
  
  # Grey boxplots as in the original publication figure
  geom_boxplot(
    fill = "grey65",
    color = "grey25",
    linewidth = 0.65,
    outlier.shape = 16,
    outlier.size = 2.2,
    outlier.colour = "grey45",
    width = 0.75
  ) +
  
  # TSS = 0.60 ensemble-selection threshold
  geom_vline(
    xintercept = 0.60,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.75
  ) +
  
  # Mean TSS shown as red point
  geom_point(
    data = tss_fig5_mean,
    aes(x = Mean_TSS, y = Model),
    inherit.aes = FALSE,
    color = "#C00000",
    size = 3.3
  ) +
  
  scale_x_continuous(
    limits = c(-0.045, 0.90),
    breaks = c(0.00, 0.25, 0.50, 0.75),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Model performance (TSS per run) - ENV",
    x = "TSS",
    y = "Model"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text.x = element_text(
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 14
    ),
    
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.45
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey30",
      fill = NA,
      linewidth = 0.6
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

print(p_tss)

save_fig(
  p_tss,
  "Figure_05_Model_performance_TSS.png",
  folder = "MAIN",
  width = 10,
  height = 6.5,
  dpi = 450
)

# Confirm exact output
figure5_file <- file.path(
  output_root,
  "MAIN",
  "Figure_05_Model_performance_TSS.png"
)

cat(
  "\nFigure 5 saved in MAIN:\n",
  normalizePath(
    figure5_file,
    mustWork = TRUE
  ),
  "\n"
)

# ==============================================================================
# FIGURE 06 — RELATIVE IMPORTANCE OF ENVIRONMENTAL PREDICTORS
# ==============================================================================

# Use variable importance already calculated from the final revised model
var_imp_fig6 <- var_imp_main %>%
  dplyr::filter(
    expl.var %in% c("BIO4", "BIO5", "BIO15", "BIO18")
  )

# Check values
print(var_imp_fig6)

# Order identical in style to the original Figure 6:
# highest importance at the top
var_imp_fig6$expl.var <- factor(
  var_imp_fig6$expl.var,
  levels = c("BIO15", "BIO5", "BIO18", "BIO4")
)

p_varimp <- ggplot(
  var_imp_fig6,
  aes(
    x = expl.var,
    y = MeanImportance,
    fill = MeanImportance
  )
) +
  
  geom_col(
    color = "grey25",
    linewidth = 0.35,
    width = 0.72
  ) +
  
  coord_flip() +
  
  scale_fill_gradient(
    low = "grey85",
    high = "grey30",
    guide = "none"
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.55),
    breaks = seq(0, 0.5, by = 0.1),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Relative importance of environmental predictors",
    subtitle = "Mean permutation importance across models and runs",
    x = "Predictor",
    y = "Mean importance"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text.x = element_text(
      size = 12
    ),
    
    axis.text.y = element_text(
      size = 13
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.4
    )
  )

print(p_varimp)

save_fig(
  p_varimp,
  "Figure_06_Variable_importance.png",
  folder = "MAIN",
  width = 9,
  height = 6,
  dpi = 450
)

cat(
  "\nFigure 6 saved in MAIN:\n",
  normalizePath(
    file.path(
      output_root,
      "MAIN",
      "Figure_06_Variable_importance.png"
    ),
    mustWork = TRUE
  ),
  "\n"
)
# ==============================================================================
# FIGURE 07 — CLIMATIC SUITABILITY ACROSS EUROPE
# ==============================================================================

p_suit <- ggplot() +
  
  geom_raster(
    data = suit_eur_df,
    aes(
      x = x,
      y = y,
      fill = Suitability
    )
  ) +
  
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    name = "Climatic\nsuitability"
  ) +
  
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.20
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = c(34, 70),
    expand = FALSE
  ) +
  
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.18, "in"),
    pad_y = grid::unit(0.18, "in"),
    height = grid::unit(0.50, "in"),
    width = grid::unit(0.40, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = expression(
      bold("Climatic suitability for ") *
        bolditalic("Actinidia arguta")
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    plot.title = element_text(
      size = 18,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(p_suit)

save_fig(
  p_suit,
  "Figure_07_Climatic_suitability_Europe.png",
  folder = "MAIN",
  width = 9,
  height = 8,
  dpi = 450
)

cat(
  "\nFigure 7 saved in MAIN:\n",
  normalizePath(
    file.path(
      output_root,
      "MAIN",
      "Figure_07_Climatic_suitability_Europe.png"
    ),
    mustWork = TRUE
  ),
  "\n"
)
# ==============================================================================
# FIGURE 08 — INTER-MODEL DISAGREEMENT ACROSS EUROPE
# ==============================================================================

sd_eur <- terra::mask(
  terra::crop(sd_map_main, europe_vect),
  europe_vect
)

sd_eur_df <- as.data.frame(
  sd_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(sd_eur_df) <- c("x", "y", "InterModel_SD")

sd_min <- min(sd_eur_df$InterModel_SD, na.rm = TRUE)
sd_max <- max(sd_eur_df$InterModel_SD, na.rm = TRUE)

p_disagreement <- ggplot() +
  
  # Background
  geom_sf(
    data = europe,
    fill = "grey96",
    color = NA
  ) +
  
  # Inter-model SD
  geom_raster(
    data = sd_eur_df,
    aes(x = x, y = y, fill = InterModel_SD)
  ) +
  
  # Country borders
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.08
  ) +
  
  # Colour scale
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(sd_min, sd_max),
    breaks = c(sd_min, 0.100, 0.200, sd_max),
    labels = c(
      sprintf("%.3f", sd_min),
      "0.100",
      "0.200",
      sprintf("%.3f", sd_max)
    ),
    name = "Inter-model\nSD"
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = EU_YLIM,
    expand = FALSE
  ) +
  
  # North arrow — same logic/style as Figure 7
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.18, "in"),
    pad_y = grid::unit(0.18, "in"),
    height = grid::unit(0.50, "in"),
    width = grid::unit(0.40, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = expression(
      bold("Inter-model disagreement for ") *
        bolditalic("Actinidia arguta") *
        bold(" projections in Europe")
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_publication(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    panel.grid = element_blank()
  )


print(p_disagreement)


save_fig(
  p_disagreement,
  "Figure_08_Intermodel_disagreement_Europe.png",
  folder = "MAIN",
  width = 10,
  height = 8,
  dpi = 450
)

# ==============================================================================
# FIGURE 09 — CLIMATIC SUITABILITY + ENVIRONMENTAL EXTRAPOLATION
# ==============================================================================

# Prepare the main-model clamping raster on the same grid as the
# European climatic-suitability projection.

clamp9 <- terra::crop(
  clamp_main,
  suit_eur
)

if (!terra::compareGeom(
  clamp9,
  suit_eur,
  stopOnError = FALSE
)) {
  
  clamp9 <- terra::resample(
    clamp9,
    suit_eur,
    method = "near"
  )
}

clamp9 <- terra::mask(
  clamp9,
  europe_vect
)

# ------------------------------------------------------------------------------
# Continuous climatic suitability
# ------------------------------------------------------------------------------

suit9_df <- as.data.frame(
  suit_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(suit9_df) <- c(
  "x",
  "y",
  "Suitability"
)

# ------------------------------------------------------------------------------
# Environmental extrapolation overlay
#
# Clamping values:
# 1 = one predictor outside the East Asian calibration range
# 2 = two predictors outside the East Asian calibration range
#
# Cells with no clamping are omitted from the overlay so that the
# continuous climatic-suitability surface remains visible.
# ------------------------------------------------------------------------------

clamp9_overlay <- terra::ifel(
  clamp9 > 0,
  clamp9,
  NA
)

clamp9_df <- as.data.frame(
  clamp9_overlay,
  xy = TRUE,
  na.rm = TRUE
)

colnames(clamp9_df) <- c(
  "x",
  "y",
  "Clamping"
)

clamp9_df$Clamping <- factor(
  round(clamp9_df$Clamping),
  levels = c(1, 2),
  labels = c(
    "1 predictor outside range",
    "2 predictors outside range"
  )
)

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p9 <- ggplot() +
  
  # Continuous climatic suitability
  geom_raster(
    data = suit9_df,
    aes(
      x = x,
      y = y,
      fill = Suitability
    )
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    name = "Climatic\nsuitability"
  ) +
  
  # Country borders
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey55",
    linewidth = 0.12
  ) +
  
  # Start a second fill scale for extrapolation
  ggnewscale::new_scale_fill() +
  
  # Environmental extrapolation
  geom_raster(
    data = clamp9_df,
    aes(
      x = x,
      y = y,
      fill = Clamping
    ),
    alpha = 0.70
  ) +
  
  scale_fill_manual(
    values = c(
      "1 predictor outside range" = "#FFB347",
      "2 predictors outside range" = "#C51B7D"
    ),
    name = "Environmental\nextrapolation",
    drop = FALSE
  ) +
  
  # Draw borders again above the extrapolation layer
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey45",
    linewidth = 0.12
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = EU_YLIM,
    expand = FALSE
  ) +
  
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.10, "in"),
    pad_y = grid::unit(0.10, "in"),
    height = grid::unit(0.32, "in"),
    width = grid::unit(0.25, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = expression(
      bold("Climatic suitability and environmental extrapolation for ") *
        bolditalic("Actinidia arguta")
    ),
    subtitle = paste0(
      "Clamping indicates predictors outside the East Asian calibration range",
      "\n",
      "Clamped area: 10.2% of Europe; 0.4% within suitability \u2265 0.50; ",
      "0% within suitability \u2265 0.60"
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_publication(base_size = 13) +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 12,
      face = "plain",
      hjust = 0,
      lineheight = 1.15
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )

# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

save_fig(
  p9,
  "Figure_09_Climatic_suitability_environmental_extrapolation.png",
  folder = "MAIN",
  width = 12,
  height = 8,
  dpi = 450
)

print(p9)

cat(
  "\nFigure 9 saved successfully:\n",
  file.path(
    output_root,
    "MAIN",
    "Figure_09_Climatic_suitability_environmental_extrapolation.png"
  ),
  "\n"
)
# ==============================================================================
# FIGURE 10 — MARCH–APRIL–MAY Tmin EXPOSURE
# ==============================================================================

tmin_min <- min(spring_long$Tmin_C, na.rm = TRUE)
tmin_max <- max(spring_long$Tmin_C, na.rm = TRUE)

tmin_breaks <- sort(unique(c(
  tmin_min,
  -10,
  0,
  10,
  tmin_max
)))

tmin_labels <- ifelse(
  tmin_breaks < 0,
  paste0("\u2212", sprintf("%.1f", abs(tmin_breaks))),
  sprintf("%.1f", tmin_breaks)
)
p_tmin <- ggplot() +
  
  geom_raster(
    data = spring_long,
    aes(
      x = x,
      y = y,
      fill = Tmin_C
    )
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(tmin_min, tmin_max),
    breaks = tmin_breaks,
    labels = tmin_labels,
    name = "Mean minimum\ntemperature (°C)"
  ) +
  
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.20
  ) +
  
  facet_wrap(
    ~ Month,
    ncol = 3
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = c(34, 70),
    expand = FALSE
  ) +
  
  labs(
    title = "Spring minimum-temperature exposure across Europe",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    # Clean scientific panel headings
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 5)
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    panel.grid = element_blank(),
    
    # Prevent longitude labels of adjacent panels from overlapping
    panel.spacing.x = grid::unit(1.2, "cm"),
    
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(p_tmin)

save_fig(
  p_tmin,
  "Figure_10_Spring_Tmin_exposure_Europe.png",
  folder = "MAIN",
  width = 12,
  height = 5.8,
  dpi = 450
)
# ==============================================================================
# FIGURE 11
# CLIMATIC SUITABILITY × SPRING LOW-TEMPERATURE EXPOSURE
#
# Columns: March | April | May
# Rows:    Tmin <= -1.0 C | Tmin <= -1.5 C
# ==============================================================================

library(terra)
library(ggplot2)
library(dplyr)
library(ggspatial)

# ------------------------------------------------------------------------------
# 1. PREPARE MONTHLY Tmin RASTERS
# ------------------------------------------------------------------------------

march_fig11 <- march_eur


april_fig11 <- terra::crop(
  tmin_april_eu,
  suit_eur
)

if (!terra::compareGeom(
  april_fig11,
  suit_eur,
  stopOnError = FALSE
)) {
  
  april_fig11 <- terra::resample(
    april_fig11,
    suit_eur,
    method = "bilinear"
  )
}

april_fig11 <- terra::mask(
  april_fig11,
  europe_vect
)


may_fig11 <- terra::crop(
  tmin_may_eu,
  suit_eur
)

if (!terra::compareGeom(
  may_fig11,
  suit_eur,
  stopOnError = FALSE
)) {
  
  may_fig11 <- terra::resample(
    may_fig11,
    suit_eur,
    method = "bilinear"
  )
}

may_fig11 <- terra::mask(
  may_fig11,
  europe_vect
)


# ------------------------------------------------------------------------------
# 2. SUITABILITY DATA
# ------------------------------------------------------------------------------

suit_fig11_df <- as.data.frame(
  suit_eur,
  xy = TRUE,
  na.rm = TRUE
)

names(suit_fig11_df)[1:3] <- c(
  "x",
  "y",
  "Suitability"
)


# ------------------------------------------------------------------------------
# 3. FUNCTION FOR Tmin CONTOURS
# ------------------------------------------------------------------------------

make_contour_df <- function(
    r,
    threshold,
    month_name,
    threshold_name
) {
  
  z <- r - threshold
  
  df <- as.data.frame(
    z,
    xy = TRUE,
    na.rm = TRUE
  )
  
  names(df)[1:3] <- c(
    "x",
    "y",
    "z"
  )
  
  df$Month <- month_name
  df$Threshold <- threshold_name
  
  return(df)
}


# ------------------------------------------------------------------------------
# 4. SIX MONTH × THRESHOLD COMBINATIONS
# ------------------------------------------------------------------------------

contour_march_10 <- make_contour_df(
  march_fig11,
  -1.0,
  "March",
  "Tmin ≤ −1.0°C"
)

contour_april_10 <- make_contour_df(
  april_fig11,
  -1.0,
  "April",
  "Tmin ≤ −1.0°C"
)

contour_may_10 <- make_contour_df(
  may_fig11,
  -1.0,
  "May",
  "Tmin ≤ −1.0°C"
)


contour_march_15 <- make_contour_df(
  march_fig11,
  -1.5,
  "March",
  "Tmin ≤ −1.5°C"
)

contour_april_15 <- make_contour_df(
  april_fig11,
  -1.5,
  "April",
  "Tmin ≤ −1.5°C"
)

contour_may_15 <- make_contour_df(
  may_fig11,
  -1.5,
  "May",
  "Tmin ≤ −1.5°C"
)


contour_fig11_df <- bind_rows(
  contour_march_10,
  contour_april_10,
  contour_may_10,
  contour_march_15,
  contour_april_15,
  contour_may_15
)


# ------------------------------------------------------------------------------
# 5. PANEL ORDER
# ------------------------------------------------------------------------------

contour_fig11_df$Month <- factor(
  contour_fig11_df$Month,
  levels = c(
    "March",
    "April",
    "May"
  )
)

contour_fig11_df$Threshold <- factor(
  contour_fig11_df$Threshold,
  levels = c(
    "Tmin ≤ −1.0°C",
    "Tmin ≤ −1.5°C"
  )
)


# ------------------------------------------------------------------------------
# 6. DUPLICATE SUITABILITY BACKGROUND FOR ALL SIX PANELS
# ------------------------------------------------------------------------------

panel_grid <- expand.grid(
  Month = c(
    "March",
    "April",
    "May"
  ),
  Threshold = c(
    "Tmin ≤ −1.0°C",
    "Tmin ≤ −1.5°C"
  ),
  stringsAsFactors = FALSE
)


suit_panels <- lapply(
  seq_len(nrow(panel_grid)),
  function(i) {
    
    x <- suit_fig11_df
    
    x$Month <- panel_grid$Month[i]
    x$Threshold <- panel_grid$Threshold[i]
    
    x
  }
)


suit_fig11_panels <- bind_rows(
  suit_panels
)


suit_fig11_panels$Month <- factor(
  suit_fig11_panels$Month,
  levels = c(
    "March",
    "April",
    "May"
  )
)

suit_fig11_panels$Threshold <- factor(
  suit_fig11_panels$Threshold,
  levels = c(
    "Tmin ≤ −1.0°C",
    "Tmin ≤ −1.5°C"
  )
)


# ------------------------------------------------------------------------------
# 7. FINAL PLOT
# ------------------------------------------------------------------------------

p_fig11_test_v2 <- ggplot() +
  
  # Continuous climatic suitability
  geom_tile(
    data = suit_fig11_panels,
    aes(
      x = x,
      y = y,
      fill = Suitability
    )
  ) +
  
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1.00
    ),
    name = "Climatic\nsuitability"
  ) +
  
  # Temperature-threshold contour
  geom_contour(
    data = contour_fig11_df,
    aes(
      x = x,
      y = y,
      z = z
    ),
    breaks = 0,
    colour = "red",
    linewidth = 0.40
  ) +
  
  # Country borders
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.18
  ) +
  
  # North arrow — repeated automatically in every facet
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.12, "in"),
    pad_y = grid::unit(0.12, "in"),
    height = grid::unit(0.34, "in"),
    width = grid::unit(0.27, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  # 2 rows × 3 columns
  facet_grid(
    Threshold ~ Month
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = c(34, 70),
    expand = FALSE
  ) +
  
  labs(
    
    title = expression(
      bold("Climatic suitability and spring temperature exposure for ") *
        bolditalic("Actinidia arguta")
    ),
    
    subtitle =
      "Red contour indicates the monthly climatological Tmin threshold",
    
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw(
    base_size = 13
  ) +
  
  theme(
    
    # --------------------------------------------------
    # MAIN TITLE
    # --------------------------------------------------
    
    plot.title = element_text(
      size = 18,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 12.5,
      hjust = 0.5,
      margin = margin(b = 9)
    ),
    
    
    # --------------------------------------------------
    # PANEL HEADINGS
    # --------------------------------------------------
    
    strip.background = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(
        t = 2,
        b = 5
      )
    ),
    
    strip.text.y = element_text(
      face = "bold",
      size = 11.5,
      margin = margin(
        l = 5,
        r = 5
      )
    ),
    
    
    # --------------------------------------------------
    # AXES
    # --------------------------------------------------
    
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text = element_text(
      size = 9
    ),
    
    
    # --------------------------------------------------
    # LEGEND
    # --------------------------------------------------
    
    legend.title = element_text(
      face = "bold",
      size = 11.5
    ),
    
    legend.text = element_text(
      size = 10.5
    ),
    
    
    # --------------------------------------------------
    # MAP APPEARANCE
    # --------------------------------------------------
    
    panel.grid = element_blank(),
    
    panel.spacing.x = grid::unit(1.4, "cm"),
    
    panel.spacing.y = grid::unit(
      0.55,
      "cm"
    ),
    
    plot.margin = margin(
      12,
      15,
      12,
      12
    )
  )


# ------------------------------------------------------------------------------
# 8. DISPLAY
# ------------------------------------------------------------------------------

print(
  p_fig11_test_v2
)


# ------------------------------------------------------------------------------
# 9. SAVE TEST ONLY
# ------------------------------------------------------------------------------

save_fig(
  p_fig11_test_v2,
  "Figure_11_Spring_temperature_overlay_2x3_v2.png",
  folder = "MAIN",
  width = 13,
  height = 9,
  dpi = 450
)

# ==============================================================================
# FIGURE 12 — GLM RESPONSE CURVES
# FINAL 4-PREDICTOR REV4.UNIFORM MODEL
# ==============================================================================

# Verify that the fitted main model contains exactly the four final predictors
final_vars_fig12 <- biomod2::get_formal_data(
  main_run$model,
  subinfo = "expl.var.names"
)

stopifnot(
  identical(
    final_vars_fig12,
    c("BIO4", "BIO5", "BIO15", "BIO18")
  )
)

# ------------------------------------------------------------------------------
# Select only the 20 GLM environmental cross-validation runs.
# The allRun model is excluded.
# ------------------------------------------------------------------------------

glm_models_fig12 <- grep(
  "_GLM$",
  biomod2::get_built_models(main_run$model),
  value = TRUE
)

glm_models_fig12 <- glm_models_fig12[
  !grepl("allRun", glm_models_fig12)
]

stopifnot(length(glm_models_fig12) == 20)

# ------------------------------------------------------------------------------
# Generate response curves from the FINAL 4-predictor model.
# Other predictors are fixed at their median values.
# ------------------------------------------------------------------------------

rc_glm_fig12 <- biomod2::bm_PlotResponseCurves(
  bm.out = main_run$model,
  models.chosen = glm_models_fig12,
  fixed.var = "median",
  show.variables = final_vars_fig12
)

df_glm_fig12 <- rc_glm_fig12$tab

# ------------------------------------------------------------------------------
# Mean response and 95% confidence interval across the 20 GLM ENV runs
# ------------------------------------------------------------------------------

response_summary_fig12 <- df_glm_fig12 %>%
  dplyr::group_by(expl.name, expl.val) %>%
  dplyr::summarise(
    mean_pred = mean(pred.val, na.rm = TRUE),
    sd_pred   = sd(pred.val, na.rm = TRUE),
    n         = dplyr::n(),
    se_pred   = sd_pred / sqrt(n),
    lower     = pmax(mean_pred - 1.96 * se_pred, 0),
    upper     = pmin(mean_pred + 1.96 * se_pred, 1),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    expl.name = factor(
      expl.name,
      levels = c("BIO4", "BIO5", "BIO15", "BIO18"),
      labels = c(
        "Temperature seasonality (BIO4)",
        "Max temperature of warmest month (BIO5)",
        "Precipitation seasonality (BIO15)",
        "Precipitation of warmest quarter (BIO18)"
      )
    )
  )

# ------------------------------------------------------------------------------
# Plot Figure 12
# ------------------------------------------------------------------------------

p_fig12 <- ggplot(
  response_summary_fig12,
  aes(x = expl.val, y = mean_pred)
) +
  
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
  
  facet_wrap(
    ~ expl.name,
    scales = "free_x",
    ncol = 2
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)
  ) +
  
  labs(
    title = "Response curves for selected environmental predictors",
    subtitle = "GLM environmental cross-validation runs: mean response with 95% confidence interval",
    x = NULL,
    y = "Predicted suitability"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "grey60"
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    panel.spacing = grid::unit(0.7, "cm")
  )

# ------------------------------------------------------------------------------
# Save FINAL Figure 12 to MAIN
# ------------------------------------------------------------------------------

save_fig(
  p_fig12,
  "Figure_12_Response_curves_GLM.png",
  folder = "MAIN",
  width = 10,
  height = 6.5,
  dpi = 450
)

print(p_fig12)

cat(
  "\nFigure 12 saved successfully:\n",
  file.path(
    output_root,
    "MAIN",
    "Figure_12_Response_curves_GLM.png"
  ),
  "\n"
)
# ==============================================================================
# SUPPLEMENTARY FIGURE S1 - DENSITY WEIGHTING SURFACE 
# ==============================================================================

bias_df <- as.data.frame(
  bias_weights,
  xy = TRUE,
  na.rm = TRUE
)

colnames(bias_df) <- c(
  "x",
  "y",
  "Weight"
)

bias_df$LogWeight <- log1p(
  bias_df$Weight
)

p_s2 <- ggplot() +
  
  geom_sf(
    data = study_region_sf,
    fill = "grey96",
    color = NA
  ) +
  
  geom_raster(
    data = bias_df,
    aes(
      x = x,
      y = y,
      fill = LogWeight
    )
  ) +
  
  geom_sf(
    data = study_region_sf,
    fill = NA,
    color = "grey55",
    linewidth = 0.15
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    name = "log(1 + occurrence-\ndensity weight)"
  ) +
  
  coord_sf(
    xlim = ASIA_XLIM,
    ylim = ASIA_YLIM,
    expand = FALSE
  ) +
  
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.12, "in"),
    pad_y = grid::unit(0.12, "in"),
    height = grid::unit(0.35, "in"),
    width = grid::unit(0.28, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = "Occurrence-density weighting surface for background sampling",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_publication(base_size = 13) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    axis.title = element_text(
      face = "plain"
    ),
    legend.title = element_text(
      face = "plain"
    )
  )

save_fig(
  p_s2,
  "Figure_S1_Occurrence_density_weighting_surface.png",
  folder = "SUPPLEMENTARY",
  width = 9,
  height = 7,
  dpi = 450
)

print(p_s2)

cat(
  "\nFigure S1 saved successfully:\n",
  file.path(
    output_root,
    "SUPPLEMENTARY",
    "Figure_S1_Occurrence_density_weighting_surface.png"
  ),
  "\n"
)

# ==============================================================================
# FIGURE S2 — QUANTITATIVE MARCH TEMPERATURE-EXPOSURE OVERLAP
# ==============================================================================
# This supplementary figure quantifies the spatial overlap between
# climatically suitable areas and the two March Tmin exposure scenarios.
# Rows = climatic suitability thresholds (>= 0.50 and >= 0.60)
# Columns = March Tmin thresholds (<= -1.0°C and <= -1.5°C)

s1_cases <- tibble::tribble(
  ~Suitability_threshold, ~Tmin_threshold,
  0.50,                  -1.0,
  0.50,                  -1.5,
  0.60,                  -1.0,
  0.60,                  -1.5
)

s1_list <- vector("list", nrow(s1_cases))

for (i in seq_len(nrow(s1_cases))) {
  
  s_thr <- s1_cases$Suitability_threshold[i]
  t_thr <- s1_cases$Tmin_threshold[i]
  
  # Suitable cells
  suitable_i <- suit_eur >= s_thr
  
  # March low-temperature exposure
  exposed_i <- march_eur <= t_thr
  
  # Classification:
  # 1 = climatically suitable but outside the exposure scenario
  # 2 = climatically suitable and within the exposure scenario
  class_i <- terra::ifel(
    suitable_i,
    terra::ifel(exposed_i, 2, 1),
    NA
  )
  
  class_df_i <- as.data.frame(
    class_i,
    xy = TRUE,
    na.rm = TRUE
  )
  
  colnames(class_df_i) <- c(
    "x",
    "y",
    "ExposureClass"
  )
  
  class_df_i$Suitability_threshold <- sprintf(
    "Suitability \u2265 %.2f",
    s_thr
  )
  
  class_df_i$Tmin_threshold <- sprintf(
    "March Tmin \u2264 %.1f\u00b0C",
    t_thr
  )
  
  s1_list[[i]] <- class_df_i
}

s1_df <- dplyr::bind_rows(s1_list)

# Factor order controls the 2 x 2 panel arrangement
s1_df$Suitability_threshold <- factor(
  s1_df$Suitability_threshold,
  levels = c(
    "Suitability \u2265 0.50",
    "Suitability \u2265 0.60"
  )
)

s1_df$Tmin_threshold <- factor(
  s1_df$Tmin_threshold,
  levels = c(
    "March Tmin \u2264 -1.0\u00b0C",
    "March Tmin \u2264 -1.5\u00b0C"
  )
)

s1_df$ExposureClass <- factor(
  s1_df$ExposureClass,
  levels = c(1, 2),
  labels = c(
    "Suitable \u2013 outside exposure",
    "Suitable \u2013 exposed"
  )
)

# --------------------------------------------------------------------------
# Percentages are taken directly from area_stats
# --------------------------------------------------------------------------

s1_labels <- area_stats %>%
  dplyr::filter(
    Month == "March",
    Suitability_threshold %in% c(0.50, 0.60),
    Tmin_exposure_threshold_C %in% c(-1.0, -1.5)
  ) %>%
  dplyr::mutate(
    Suitability_threshold = sprintf(
      "Suitability \u2265 %.2f",
      Suitability_threshold
    ),
    Tmin_threshold = sprintf(
      "March Tmin \u2264 %.1f\u00b0C",
      Tmin_exposure_threshold_C
    ),
    Label = sprintf(
      "%.1f%% of suitable area exposed",
      Percent_of_suitable_area_exposed
    )
  )

s1_labels$Suitability_threshold <- factor(
  s1_labels$Suitability_threshold,
  levels = c(
    "Suitability \u2265 0.50",
    "Suitability \u2265 0.60"
  )
)

s1_labels$Tmin_threshold <- factor(
  s1_labels$Tmin_threshold,
  levels = c(
    "March Tmin \u2264 -1.0\u00b0C",
    "March Tmin \u2264 -1.5\u00b0C"
  )
)

# Fixed label position in geographic coordinates
s1_labels$x <- -8.5
s1_labels$y <- 35.5

# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------

p_s1 <- ggplot() +
  
  geom_sf(
    data = europe,
    fill = "grey96",
    color = "grey75",
    linewidth = 0.12
  ) +
  
  geom_tile(
    data = s1_df,
    aes(
      x = x,
      y = y,
      fill = ExposureClass
    )
  ) +
  
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey55",
    linewidth = 0.12
  ) +
  
  facet_grid(
    rows = vars(Suitability_threshold),
    cols = vars(Tmin_threshold)
  ) +
  
  scale_fill_manual(
    values = c(
      "Suitable – outside exposure" = "grey72",
      "Suitable – exposed" = "#B2180B"
    ),
    name = NULL,
    drop = FALSE
  ) +
  
  geom_label(
    data = s1_labels,
    aes(
      x = x,
      y = y,
      label = Label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0,
    size = 3.7,
    fontface = "bold",
    fill = "white",
    label.size = 0.25
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = EU_YLIM,
    expand = FALSE
  ) +
  
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.08, "in"),
    pad_y = grid::unit(0.08, "in"),
    height = grid::unit(0.25, "in"),
    width = grid::unit(0.20, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = expression(
      bold("Climatic suitability and March temperature exposure for ") *
        bolditalic("Actinidia arguta")
    ),
    subtitle = "Quantitative overlap under two suitability and two Tmin thresholds",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_publication(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    strip.text = element_text(
      size = 11.5,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    panel.spacing.x = grid::unit(1.2, "cm"),
    
    legend.position = "bottom",
    legend.text = element_text(
      size = 11
    )
  )
save_fig(
  p_s1,
  "Figure_S2_Quantitative overlap between climatically suitable areas and March low-temperature exposure scenarios for A. arguta in Europe.png",
  folder = "SUPPLEMENTARY",
  width = 12,
  height = 9,
  dpi = 450
)

print(p_s1)

cat(
  "\nFigure S2 saved successfully:\n",
  file.path(
    output_root,
    "SUPPLEMENTARY",
    "Figure_S2_Quantitative overlap between climatically suitable areas and March low-temperature exposure scenarios for A. arguta in Europe.png"
  ),
  "\n"
)

print(
  s1_labels %>%
    dplyr::select(
      Suitability_threshold,
      Tmin_threshold,
      Suitable_area_km2,
      Suitable_area_exposed_km2,
      Percent_of_suitable_area_exposed
    )
)



# ==============================================================================
# FIGURE S3 — CLAMPING / EXTRAPOLATION
# BIAS-WEIGHTED BACKGROUND SENSITIVITY ANALYSIS
# ==============================================================================

# Locate the clamping mask produced by the bias-weighted BIOMOD2 projection
bias_clamp_files <- list.files(
  biomod_root,
  pattern = "ClampingMask\\.tif$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

bias_clamp_candidates <- bias_clamp_files[
  grepl(
    "MiniKiwi_REV4_BIASWEIGHT",
    bias_clamp_files,
    fixed = TRUE
  ) &
    !grepl(
      "UNIFORM",
      bias_clamp_files,
      fixed = TRUE
    )
]

if (length(bias_clamp_candidates) == 0) {
  stop("Bias-weighted clamping mask was not found.")
}

# Show exactly which file is being used
cat(
  "\nBias-weighted clamping source:\n",
  bias_clamp_candidates[1],
  "\n"
)

clamp_bias <- terra::rast(
  bias_clamp_candidates[1]
)

# Crop and mask to Europe
clamp_bias_eur <- terra::crop(
  clamp_bias,
  suit_eur
)

if (!terra::compareGeom(
  clamp_bias_eur,
  suit_eur,
  stopOnError = FALSE
)) {
  
  clamp_bias_eur <- terra::resample(
    clamp_bias_eur,
    suit_eur,
    method = "near"
  )
}

clamp_bias_eur <- terra::mask(
  clamp_bias_eur,
  europe_vect
)

# Export raster for reproducibility
terra::writeRaster(
  clamp_bias_eur,
  file.path(
    output_root,
    "RASTERS",
    "clamping_mask_REV4_biasweighted_Europe.tif"
  ),
  overwrite = TRUE
)

# Convert to data frame
clamp_bias_df <- as.data.frame(
  clamp_bias_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(clamp_bias_df) <- c(
  "x",
  "y",
  "N_outside_range"
)

# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------

p_s3 <- ggplot() +
  
  geom_sf(
    data = europe,
    fill = "grey96",
    color = NA
  ) +
  
  geom_raster(
    data = clamp_bias_df,
    aes(
      x = x,
      y = y,
      fill = N_outside_range
    )
  ) +
  
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey55",
    linewidth = 0.12
  ) +
  
  scale_fill_viridis_c(
    option = "cividis",
    limits = c(0, 2),
    breaks = c(0, 0.5, 1.0, 1.5, 2.0),
    name = "Predictors outside\ncalibration range"
  ) +
  
  coord_sf(
    xlim = EU_XLIM,
    ylim = EU_YLIM,
    expand = FALSE
  ) +
  
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.12, "in"),
    pad_y = grid::unit(0.12, "in"),
    height = grid::unit(0.35, "in"),
    width = grid::unit(0.28, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  
  labs(
    title = "Clamping / extrapolation diagnostic for the European projection",
    subtitle = "Bias-weighted background sensitivity analysis",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_publication(base_size = 13) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 12,
      face = "plain",
      hjust = 0
    ),
    axis.title = element_text(
      face = "plain"
    ),
    legend.title = element_text(
      face = "plain"
    )
  )

save_fig(
  p_s3,
  "Figure_S3_Bias_weighted_clamping_extrapolation.png",
  folder = "SUPPLEMENTARY",
  width = 9,
  height = 7,
  dpi = 450
)

print(p_s3)

cat(
  "\nFigure S3 saved successfully:\n",
  file.path(
    output_root,
    "SUPPLEMENTARY",
    "Figure_S3_Bias_weighted_clamping_extrapolation.png"
  ),
  "\n"
)
# ==============================================================================
# TABLE S1 — OCCURRENCE DATA FILTERING
# ==============================================================================

table_s1 <- filtering_summary %>%
  dplyr::rename(
    "Filtering step" = Step,
    "Number of records" = N
  )

save_xlsx(
  list(
    "Table S1" = table_s1
  ),
  "Table_S1_Occurrence_data_filtering.xlsx"
)

readr::write_csv(
  table_s1,
  file.path(
    output_root,
    "TABLES",
    "Table_S1_Occurrence_data_filtering.csv"
  )
)

cat(
  "\nTable S1 saved successfully:\n",
  file.path(
    output_root,
    "TABLES",
    "Table_S1_Occurrence_data_filtering.xlsx"
  ),
  "\n"
)

print(table_s1)
# ==============================================================================
# TABLE S2
# ==============================================================================

# Start with all candidate predictors:
# BIO1-BIO19 + spring_min + chilling_proxy
vif_work <- candidate_vals

vif_removed_list <- list()
step_id <- 1

# ------------------------------------------------------------------------------
# Iterative VIF screening
# Predictors with VIF > 10 are sequentially removed
# ------------------------------------------------------------------------------

repeat {
  
  current_vif <- vif_from_df(vif_work)
  
  max_row <- current_vif[
    which.max(current_vif$VIF),
    ,
    drop = FALSE
  ]
  
  max_vif <- max_row$VIF
  max_var <- max_row$Variable
  
  # Stop when all remaining predictors have VIF <= 10
  if (!is.finite(max_vif) || max_vif <= 10) {
    break
  }
  
  # Record VIF at the iteration in which the predictor was excluded
  vif_removed_list[[step_id]] <- tibble::tibble(
    Predictor = max_var,
    VIF = max_vif,
    Screening_step = step_id,
    VIF_screening_outcome = "Excluded"
  )
  
  # Remove predictor
  vif_work <- vif_work[
    ,
    setdiff(names(vif_work), max_var),
    drop = FALSE
  ]
  
  step_id <- step_id + 1
}

# ------------------------------------------------------------------------------
# Predictors retained after completion of VIF screening
# ------------------------------------------------------------------------------

vif_retained_final <- vif_from_df(vif_work) %>%
  dplyr::transmute(
    Predictor = Variable,
    VIF = VIF,
    Screening_step = NA_integer_,
    VIF_screening_outcome = "Retained"
  )

# ------------------------------------------------------------------------------
# Combine excluded and retained predictors
# ------------------------------------------------------------------------------

table_s2_vif <- dplyr::bind_rows(
  dplyr::bind_rows(vif_removed_list),
  vif_retained_final
)

# Final predictors used in the revised SDM
final_sdm_predictors <- c(
  "BIO4",
  "BIO5",
  "BIO15",
  "BIO18"
)

table_s2_vif <- table_s2_vif %>%
  dplyr::mutate(
    Final_SDM = ifelse(
      Predictor %in% final_sdm_predictors,
      "Included",
      "Not included"
    )
  )

# ------------------------------------------------------------------------------
# Restore original candidate-predictor order
# ------------------------------------------------------------------------------

original_order <- names(candidate_vals)

table_s2_vif$Predictor <- factor(
  table_s2_vif$Predictor,
  levels = original_order
)

table_s2_vif <- table_s2_vif %>%
  dplyr::arrange(Predictor) %>%
  dplyr::mutate(
    Predictor = as.character(Predictor),
    VIF = round(VIF, 2)
  ) %>%
  dplyr::rename(
    "Screening step" = Screening_step,
    "VIF screening outcome" = VIF_screening_outcome,
    "Final SDM" = Final_SDM
  )

# ------------------------------------------------------------------------------
# Verification against the actual usdm::vifstep result
# ------------------------------------------------------------------------------

retained_reconstructed <- names(vif_work)
retained_usdm <- vif_candidates_step@results$Variables

if (!setequal(retained_reconstructed, retained_usdm)) {
  stop(
    "ERROR: Reconstructed VIF screening does not match usdm::vifstep."
  )
}

cat("\nVIF screening verified successfully.\n")
cat(
  "Predictors retained after VIF screening:\n",
  paste(retained_reconstructed, collapse = ", "),
  "\n\n"
)

# ------------------------------------------------------------------------------
# Print final Supplementary Table S2
# ------------------------------------------------------------------------------

print(
  table_s2_vif,
  n = Inf
)

# ------------------------------------------------------------------------------
# Export Supplementary Table S2
# ------------------------------------------------------------------------------

save_csv(
  table_s2_vif,
  "Table_S2_VIF_screening.csv"
)

save_xlsx(
  list(
    "Table S2" = table_s2_vif
  ),
  "Table_S2_VIF_screening.xlsx"
)

cat(
  "\nSupplementary Table S2 saved successfully.\n"
)

# ==============================================================================
# TABLE S3 — BIAS-WEIGHTED BACKGROUND SENSITIVITY ANALYSIS
# ==============================================================================

table_s3 <- tss_bias_summary %>%
  dplyr::left_join(
    roc_bias_summary,
    by = "algo"
  ) %>%
  dplyr::transmute(
    Algorithm = algo,
    `Mean TSS` = round(Mean_TSS, 3),
    `SD TSS` = round(SD_TSS, 3),
    `Maximum TSS` = round(Max_TSS, 3),
    `Mean AUC` = round(Mean_AUC, 3),
    `SD AUC` = round(SD_AUC, 3),
    `Maximum AUC` = round(Max_AUC, 3)
  ) %>%
  dplyr::arrange(dplyr::desc(`Mean TSS`))

save_xlsx(
  list(
    "Table S3" = table_s3,
    "Threshold assessment" = bias_threshold_result
  ),
  "Table_S3_Bias_weighted_background_sensitivity.xlsx"
)

readr::write_csv(
  table_s3,
  file.path(
    output_root,
    "TABLES",
    "Table_S3_Bias_weighted_background_sensitivity.csv"
  )
)

cat(
  "\nTable S3 saved successfully:\n",
  file.path(
    output_root,
    "TABLES",
    "Table_S3_Bias_weighted_background_sensitivity.xlsx"
  ),
  "\n"
)

print(table_s3)
print(bias_threshold_result)
# ==============================================================================
# 14. COMBINED WORKBOOK + SESSION INFO
# ==============================================================================

save_xlsx(
  list(
    Occurrence_filtering = filtering_summary,
    Final_country_counts = pres_country,
    VIF_legacy5 = vif_legacy5_explicit,
    VIF_revised4 = vif_final4_explicit,
    Main_TSS = tss_summary_main,
    Main_AUC = roc_summary_main,
    Variable_importance = var_imp_main,
    Spring_area_exposure = area_stats,
    Spring_exposure_classes = class_area,
    Bias_threshold_result = bias_threshold_result,
    Bias_TSS = tss_bias_summary,
    Bias_AUC = roc_bias_summary
  ),
  "Reviewer_revision_key_results.xlsx"
)

capture.output(sessionInfo(), file = file.path(output_root, "TABLES", "R_sessionInfo.txt"))
cat("Spring exposure:        March, April, May climatological Tmin\n")

cat("\n=================================================================\n")
cat("REVISION ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Project root:          ", project_root, "\n")
cat("Outputs:               ", output_root, "\n")
cat("Final presence N:      ", n_final, "\n")
cat("Main predictors:       ", paste(final_keep, collapse = ", "), "\n")
cat("Main background:        uniform random (10,000)\n")
cat("Sensitivity background: occurrence-density weighted (10,000)\n")
cat("Temperature scaling:    WorldClim 2.1 °C, NO /10\n")
cat("=================================================================\n")
