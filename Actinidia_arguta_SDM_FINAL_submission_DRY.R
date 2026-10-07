# ==============================================================================
# Actinidia arguta Species Distribution Modelling — Comprehensive Reproducible Analysis
# ==============================================================================
#
# PURPOSE:
# Complete workflow for species distribution modeling incorporating reviewer feedback:
#
#   1. Occurrence data cleaning with transparent step-by-step counts.
#   2. Corrected WorldClim 2.1 temperature scaling (no division by 10).
#   3. VIF-based predictor selection: final 4-variable set (BIO4, BIO5, BIO15, BIO18).
#   4. European climatic suitability projections with environmental cross-validation.
#   5. Spring temperature-exposure assessment (March–May monthly Tmin).
#   6. Sensitivity analysis: occurrence-density-weighted background sampling.
#   7. Model diagnostics: inter-model disagreement, clamping/extrapolation masking.
#   8. Publication-ready figures and supplementary tables.
#
# USAGE:
# Place this script in the project root directory alongside occurrence.txt and data_env/ folder.
# In R/RStudio: Source this entire script. Outputs are written to REVISION_2026_EXPORTS/.
#
# ==============================================================================

rm(list = ls())
graphics.off()
cat("\014")

set.seed(123)
options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. PACKAGES AND SETTINGS
# ==============================================================================

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "sf", "terra", "geodata", "usdm",
  "biomod2", "ggplot2", "ggspatial", "scales", "viridis", "rnaturalearth",
  "rnaturalearthdata", "pROC", "spThin", "maxnet", "writexl", "reshape2", "ggnewscale"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, then rerun the script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

SEED <- 123
N_BACKGROUND <- 10000
N_DOMAIN_CORR <- 100000
SUIT_THRESHOLDS <- c(0.50, 0.60)
TEMP_THRESHOLDS <- c(-1.0, -1.5)

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

terra::terraOptions(
  tempdir = file.path(output_root, "TEMP"),
  progress = 1,
  memfrac = 0.60
)

occurrence_candidates <- c(
  file.path(project_root, "occurrence.txt"),
  file.path(project_root, "occurrence(1).txt")
)
occurrence_file <- occurrence_candidates[file.exists(occurrence_candidates)][1]
if (is.na(occurrence_file) || length(occurrence_file) == 0) {
  stop("occurrence.txt not found in project root directory.")
}

data_env_dir <- file.path(project_root, "data_env")
dir.create(data_env_dir, recursive = TRUE, showWarnings = FALSE)

biomod_root <- file.path(output_root, "BIOMOD_WORK")

# ==============================================================================
# 2. HELPER FUNCTIONS
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

  message("\n", resp_name)
  message("Presences: ", sum(resp_var2 == 1))
  message("Background: ", sum(resp_var2 == 0))

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
          "\nAttempting to recover ensemble projection from disk."
        )

        emproj_files <- list.files(
          file.path(biomod_root, resp_name),
          pattern = "\\.ensemble\\.projection\\.out$",
          recursive = TRUE,
          full.names = TRUE
        )

        if (length(emproj_files) == 0) {
          stop("Ensemble forecasting failed and no saved ensemble projection was found.")
        }

        emproj_info <- file.info(emproj_files)
        emproj_file <- rownames(emproj_info)[which.max(emproj_info$mtime)]

        loaded_object_name <- load(emproj_file)
        recovered <- get(loaded_object_name)

        if (!inherits(recovered, "BIOMOD.projection.out")) {
          stop("Recovered ensemble projection is not a BIOMOD.projection.out object.")
        }

        message("Recovered ensemble projection from: ", emproj_file)

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
# 3. OCCURRENCE DATA CLEANING
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
# 4. WORLDCLIM 2.1 DATA — CORRECTED TEMPERATURE SCALING
# ==============================================================================

bio_global  <- geodata::worldclim_global(var = "bio",  res = 2.5, path = data_env_dir)
tmin_global <- geodata::worldclim_global(var = "tmin", res = 2.5, path = data_env_dir)
tavg_global <- geodata::worldclim_global(var = "tavg", res = 2.5, path = data_env_dir)

bio_asia <- terra::crop(bio_global, extent_asia)
names(bio_asia) <- paste0("BIO", 1:19)

spring_min <- terra::app(terra::crop(tmin_global[[3:4]], extent_asia), fun = min)
names(spring_min) <- "spring_min"

chilling_proxy <- terra::app(
  terra::crop(tavg_global[[c(1, 2, 3, 11, 12)]], extent_asia),
  fun = mean
)
names(chilling_proxy) <- "chilling_proxy"

candidate_stack <- c(bio_asia, spring_min, chilling_proxy)

cell_template <- bio_asia[[1]]
thin_cells <- terra::cellFromXY(cell_template, as.matrix(pres_thinned[, c("x", "y")]))
keep_cell <- !is.na(thin_cells) & !duplicated(thin_cells)
pres_raster_filtered <- pres_thinned[keep_cell, , drop = FALSE]
n_raster_filtered <- nrow(pres_raster_filtered)

legacy5_keep <- c("spring_min", "BIO4", "BIO5", "BIO15", "BIO18")
final_keep <- c("BIO4", "BIO5", "BIO15", "BIO18")

legacy5_stack <- candidate_stack[[legacy5_keep]]
bio_selected <- candidate_stack[[final_keep]]

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
    "Records with valid coordinates",
    "Unique coordinate pairs",
    "Terrestrial records",
    "Records in East Asia calibration region",
    "After 5-km spatial thinning",
    "After one-record-per-cell filtering",
    "Final presences with complete predictor data"
  ),
  N = c(n_raw, n_coords, n_unique, n_land, n_region, n_thinned, n_raster_filtered, n_final)
)
save_csv(filtering_summary, "Occurrence_filtering_summary.csv")

pres_final_sf <- st_as_sf(pres_final, coords = c("x", "y"), crs = 4326, remove = FALSE)
pres_country <- st_join(
  pres_final_sf,
  world_sf[, c("admin", "iso_a2")],
  join = st_within,
  left = TRUE
) %>%
  st_drop_geometry() %>%
  mutate(
    Region = case_when(
      admin == "China" ~ "China",
      admin == "Japan" ~ "Japan",
      grepl("Korea", admin, ignore.case = TRUE) ~ "Korea",
      TRUE ~ "Other"
    )
  ) %>%
  count(Region, name = "N") %>%
  arrange(desc(N))
save_csv(pres_country, "Final_presence_distribution_by_region.csv")

raw_country_codes <- gbif_data %>% count(countryCode, sort = TRUE, name = "N")
save_csv(raw_country_codes, "Raw_GBIF_country_codes.csv")

# ==============================================================================
# 5. MULTICOLLINEARITY DIAGNOSTICS
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
         "Correlation_final4_occurrence_points.csv")
save_csv(as.data.frame(cor_legacy5_domain) %>% tibble::rownames_to_column("Variable"),
         "Correlation_legacy5_calibration_domain.csv")
save_csv(as.data.frame(cor_final4_domain) %>% tibble::rownames_to_column("Variable"),
         "Correlation_final4_calibration_domain.csv")
save_csv(vif_legacy5_explicit, "VIF_legacy5_explicit.csv")
save_csv(vif_final4_explicit, "VIF_final4_explicit.csv")
save_csv(tibble::as_tibble(vif_final4_step@results), "VIF_final4_vifstep.csv")
save_csv(tibble::as_tibble(vif_candidates_step@results), "VIF_candidate_predictors_vifstep.csv")

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

bg_uniform <- sample_background(
  sampling_raster = asia_land_mask,
  n = N_BACKGROUND,
  pres_xy = pres_final,
  template = bio_selected[[1]],
  method = "random",
  seed = SEED
)
save_csv(bg_uniform, "Background_uniform_10000.csv")

points_region_vect <- terra::vect(points_region)
bias_count <- terra::rasterize(
  points_region_vect,
  bio_selected[[1]],
  field = 1,
  fun = "sum",
  background = 0
)

bias_smooth <- terra::focal(
  bias_count,
  w = matrix(1, 25, 25),
  fun = "sum",
  na.policy = "omit",
  fillvalue = 0
)

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

tmin_spring_eu <- terra::crop(tmin_global[[3:5]], extent_eu)
names(tmin_spring_eu) <- c("MarchTmin_C", "AprilTmin_C", "MayTmin_C")

tmin_march_eu <- tmin_spring_eu[["MarchTmin_C"]]
tmin_april_eu <- tmin_spring_eu[["AprilTmin_C"]]
tmin_may_eu   <- tmin_spring_eu[["MayTmin_C"]]
names(tmin_march_eu) <- "MarchTmin_C"

# ==============================================================================
# 8. MAIN REVISED SDM — FOUR PREDICTORS + UNIFORM BACKGROUND
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

clamp_files <- list.files(
  biomod_root,
  pattern = "clamp.*\\.(tif|grd)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

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
# 11. BIAS-SENSITIVITY RUN — DENSITY-WEIGHTED BACKGROUND
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

evals_bias <- get_evaluations(bias_run$model)

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

# Figure 1: Occurrence records
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

# Figure 2: Candidate predictor correlations
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

n_cor <- nrow(candidate_bio_vals)

t_cor <- cor_bio_candidates * sqrt(
  (n_cor - 2) / (1 - cor_bio_candidates^2)
)

p_bio_candidates <- 2 * pt(
  -abs(t_cor),
  df = n_cor - 2
)

diag(p_bio_candidates) <- NA

save_csv(
  as.data.frame(p_bio_candidates) %>%
    tibble::rownames_to_column("Variable"),
  "Pvalues_candidate_predictor_correlations.csv"
)

strong_cor_mask <- abs(cor_bio_candidates) >= 0.75 &
  abs(cor_bio_candidates) < 0.999

strong_p_values <- p_bio_candidates[strong_cor_mask]

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

# Figure 3: Final predictor correlation matrix
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
      face = "plain"
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

# Prepare European map data
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

# Figure 4: GLM ROC curve
glm_obs <- biomod2::get_formal_data(
  main_run$model,
  subinfo = "resp.var"
)

glm_obs <- as.numeric(glm_obs)

glm_pred <- biomod2::get_predictions(
  main_run$model,
  algo = "GLM",
  run = paste0("RUN", 1:20),
  model.as.col = TRUE
)

glm_pred <- as.data.frame(glm_pred)

glm_cv <- as.data.frame(
  biomod2::get_calib_lines(main_run$model)
)

glm_cv <- glm_cv[
  ,
  paste0("_allData_RUN", 1:20),
  drop = FALSE
]

if (nrow(glm_pred) != length(glm_obs)) {
  stop(
    "Prediction rows (", nrow(glm_pred),
    ") do not match response length (",
    length(glm_obs), ")."
  )
}

fpr_grid <- seq(0, 1, length.out = 501)

roc_list <- vector("list", 20)
auc_values <- numeric(20)

for (i in 1:20) {

  validation_index <- which(glm_cv[[i]] %in% FALSE)

  obs_i <- glm_obs[validation_index]
  pred_i <- glm_pred[validation_index, i]

  keep_i <- is.finite(obs_i) & is.finite(pred_i)

  obs_i <- obs_i[keep_i]
  pred_i <- pred_i[keep_i]

  if (length(unique(obs_i)) < 2) {
    stop("RUN", i, " does not contain both response classes.")
  }

  roc_i <- pROC::roc(
    response = obs_i,
    predictor = pred_i,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  auc_values[i] <- as.numeric(pROC::auc(roc_i))

  roc_tmp <- data.frame(
    FPR = 1 - roc_i$specificities,
    TPR = roc_i$sensitivities
  )

  roc_tmp <- roc_tmp %>%
    dplyr::group_by(FPR) %>%
    dplyr::summarise(TPR = max(TPR, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(FPR)

  roc_tmp <- dplyr::bind_rows(
    data.frame(FPR = 0, TPR = 0),
    roc_tmp,
    data.frame(FPR = 1, TPR = 1)
  ) %>%
    dplyr::group_by(FPR) %>%
    dplyr::summarise(TPR = max(TPR, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(FPR)

  interp_tpr <- approx(
    x = roc_tmp$FPR,
    y = roc_tmp$TPR,
    xout = fpr_grid,
    method = "linear",
    ties = max,
    rule = 2
  )$y

  roc_list[[i]] <- data.frame(
    Run = paste0("RUN", i),
    FPR = fpr_grid,
    TPR = interp_tpr
  )
}

roc_all <- dplyr::bind_rows(roc_list)

roc_summary_glm <- roc_all %>%
  dplyr::group_by(FPR) %>%
  dplyr::summarise(
    Mean_TPR = mean(TPR, na.rm = TRUE),
    Lower_TPR = stats::quantile(TPR, probs = 0.025, na.rm = TRUE, names = FALSE),
    Upper_TPR = stats::quantile(TPR, probs = 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

glm_auc_biomod <- roc_main %>%
  dplyr::filter(algo == "GLM", is.finite(validation))

glm_auc_mean <- mean(glm_auc_biomod$validation, na.rm = TRUE)
glm_auc_sd <- sd(glm_auc_biomod$validation, na.rm = TRUE)

p_roc <- ggplot(
  roc_summary_glm,
  aes(x = FPR, y = Mean_TPR)
) +
  geom_ribbon(
    aes(ymin = Lower_TPR, ymax = Upper_TPR),
    fill = "#3182BD",
    alpha = 0.25
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey35"
  ) +
  geom_line(
    color = "#B30000",
    linewidth = 1.20
  ) +
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
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1)),
    expand = c(0, 0)
  ) +
  coord_equal() +
  labs(
    title = "ROC curve for the GLM model",
    subtitle = "Environmental cross-validation runs",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(base_size = 14) +
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
    plot.margin = margin(12, 15, 12, 12)
  )

save_fig(
  p_roc,
  "Figure_04_ROC_curve_GLM.png",
  folder = "MAIN",
  width = 7.5,
  height = 7.5,
  dpi = 450
)

# Figure 5: Model performance TSS
tss_fig5 <- tss_main %>%
  dplyr::filter(is.finite(validation), run != "allRun") %>%
  dplyr::transmute(
    Model = algo,
    Run   = run,
    TSS   = validation
  )

tss_fig5$Model <- factor(
  tss_fig5$Model,
  levels = c("RF", "MAXNET", "GLM", "GBM", "GAM")
)

tss_fig5_mean <- tss_fig5 %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(Mean_TSS = mean(TSS, na.rm = TRUE), .groups = "drop")

p_tss <- ggplot(
  tss_fig5,
  aes(x = TSS, y = Model)
) +
  geom_boxplot(
    fill = "grey65",
    color = "grey25",
    linewidth = 0.65,
    outlier.shape = 16,
    outlier.size = 2.2,
    outlier.colour = "grey45",
    width = 0.75
  ) +
  geom_vline(
    xintercept = 0.60,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.75
  ) +
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
    plot.margin = margin(10, 15, 10, 10)
  )

save_fig(
  p_tss,
  "Figure_05_Model_performance_TSS.png",
  folder = "MAIN",
  width = 10,
  height = 6.5,
  dpi = 450
)

# Figure 6: Variable importance
var_imp_fig6 <- var_imp_main %>%
  dplyr::filter(expl.var %in% c("BIO4", "BIO5", "BIO15", "BIO18"))

var_imp_fig6$expl.var <- factor(
  var_imp_fig6$expl.var,
  levels = c("BIO15", "BIO5", "BIO18", "BIO4")
)

p_varimp <- ggplot(
  var_imp_fig6,
  aes(x = expl.var, y = MeanImportance, fill = MeanImportance)
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

save_fig(
  p_varimp,
  "Figure_06_Variable_importance.png",
  folder = "MAIN",
  width = 9,
  height = 6,
  dpi = 450
)

# Figure 7: Climatic suitability across Europe
p_suit <- ggplot() +
  geom_raster(
    data = suit_eur_df,
    aes(x = x, y = y, fill = Suitability)
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
    plot.margin = margin(12, 15, 12, 12)
  )

save_fig(
  p_suit,
  "Figure_07_Climatic_suitability_Europe.png",
  folder = "MAIN",
  width = 9,
  height = 8,
  dpi = 450
)

# Figure 8: Inter-model disagreement
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

p_disagreement <- ggplot() +
  geom_sf(
    data = europe,
    fill = "grey96",
    color = NA
  ) +
  geom_raster(
    data = sd_eur_df,
    aes(x = x, y = y, fill = InterModel_SD)
  ) +
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.08
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Inter-model\nSD"
  ) +
  coord_sf(
    xlim = EU_XLIM,
    ylim = EU_YLIM,
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

save_fig(
  p_disagreement,
  "Figure_08_Intermodel_disagreement_Europe.png",
  folder = "MAIN",
  width = 10,
  height = 8,
  dpi = 450
)

# Figure 9: Suitability and extrapolation
clamp9 <- terra::crop(clamp_main, suit_eur)

if (!terra::compareGeom(clamp9, suit_eur, stopOnError = FALSE)) {
  clamp9 <- terra::resample(clamp9, suit_eur, method = "near")
}

clamp9 <- terra::mask(clamp9, europe_vect)

suit9_df <- as.data.frame(
  suit_eur,
  xy = TRUE,
  na.rm = TRUE
)

colnames(suit9_df) <- c("x", "y", "Suitability")

clamp9_overlay <- terra::ifel(clamp9 > 0, clamp9, NA)

clamp9_df <- as.data.frame(
  clamp9_overlay,
  xy = TRUE,
  na.rm = TRUE
)

colnames(clamp9_df) <- c("x", "y", "Clamping")

clamp9_df$Clamping <- factor(
  round(clamp9_df$Clamping),
  levels = c(1, 2),
  labels = c(
    "1 predictor outside range",
    "2 predictors outside range"
  )
)

p9 <- ggplot() +
  geom_raster(
    data = suit9_df,
    aes(x = x, y = y, fill = Suitability)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    name = "Climatic\nsuitability"
  ) +
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey55",
    linewidth = 0.12
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = clamp9_df,
    aes(x = x, y = y, fill = Clamping),
    alpha = 0.45
  ) +
  scale_fill_manual(
    values = c(
      "1 predictor outside range" = "#F46D43",
      "2 predictors outside range" = "#B2182B"
    ),
    name = "Environmental\nextrapolation",
    drop = FALSE
  ) +
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
    subtitle = "Clamping indicates predictors outside the East Asian calibration range",
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
      hjust = 0
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )

save_fig(
  p9,
  "Figure_09_Climatic_suitability_environmental_extrapolation.png",
  folder = "MAIN",
  width = 12,
  height = 8,
  dpi = 450
)

# Figure 10: Spring Tmin exposure
p_tmin <- ggplot() +
  geom_raster(
    data = spring_long,
    aes(x = x, y = y, fill = Tmin_C)
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = expression(
      "Mean minimum\ntemperature (" * degree * "C)"
    )
  ) +
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.20
  ) +
  facet_wrap(~ Month, ncol = 3) +
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
    panel.spacing.x = grid::unit(1.2, "cm"),
    plot.margin = margin(12, 15, 12, 12)
  )

save_fig(
  p_tmin,
  "Figure_10_Spring_Tmin_exposure_Europe.png",
  folder = "MAIN",
  width = 12,
  height = 5.8,
  dpi = 450
)

# Figure 11: Suitability × Spring temperature exposure
march_fig11 <- march_eur

april_fig11 <- terra::crop(tmin_april_eu, suit_eur)

if (!terra::compareGeom(april_fig11, suit_eur, stopOnError = FALSE)) {
  april_fig11 <- terra::resample(april_fig11, suit_eur, method = "bilinear")
}

april_fig11 <- terra::mask(april_fig11, europe_vect)

may_fig11 <- terra::crop(tmin_may_eu, suit_eur)

if (!terra::compareGeom(may_fig11, suit_eur, stopOnError = FALSE)) {
  may_fig11 <- terra::resample(may_fig11, suit_eur, method = "bilinear")
}

may_fig11 <- terra::mask(may_fig11, europe_vect)

suit_fig11_df <- as.data.frame(
  suit_eur,
  xy = TRUE,
  na.rm = TRUE
)

names(suit_fig11_df)[1:3] <- c("x", "y", "Suitability")

make_contour_df <- function(r, threshold, month_name, threshold_name) {
  z <- r - threshold
  df <- as.data.frame(z, xy = TRUE, na.rm = TRUE)
  names(df)[1:3] <- c("x", "y", "z")
  df$Month <- month_name
  df$Threshold <- threshold_name
  return(df)
}

contour_march_10 <- make_contour_df(march_fig11, -1.0, "March", "Tmin ≤ −1.0°C")
contour_april_10 <- make_contour_df(april_fig11, -1.0, "April", "Tmin ≤ −1.0°C")
contour_may_10 <- make_contour_df(may_fig11, -1.0, "May", "Tmin ≤ −1.0°C")
contour_march_15 <- make_contour_df(march_fig11, -1.5, "March", "Tmin ≤ −1.5°C")
contour_april_15 <- make_contour_df(april_fig11, -1.5, "April", "Tmin ≤ −1.5°C")
contour_may_15 <- make_contour_df(may_fig11, -1.5, "May", "Tmin ≤ −1.5°C")

contour_fig11_df <- bind_rows(
  contour_march_10, contour_april_10, contour_may_10,
  contour_march_15, contour_april_15, contour_may_15
)

contour_fig11_df$Month <- factor(
  contour_fig11_df$Month,
  levels = c("March", "April", "May")
)

contour_fig11_df$Threshold <- factor(
  contour_fig11_df$Threshold,
  levels = c("Tmin ≤ −1.0°C", "Tmin ≤ −1.5°C")
)

panel_grid <- expand.grid(
  Month = c("March", "April", "May"),
  Threshold = c("Tmin ≤ −1.0°C", "Tmin ≤ −1.5°C"),
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

suit_fig11_panels <- bind_rows(suit_panels)

suit_fig11_panels$Month <- factor(
  suit_fig11_panels$Month,
  levels = c("March", "April", "May")
)

suit_fig11_panels$Threshold <- factor(
  suit_fig11_panels$Threshold,
  levels = c("Tmin ≤ −1.0°C", "Tmin ≤ −1.5°C")
)

p_fig11 <- ggplot() +
  geom_tile(
    data = suit_fig11_panels,
    aes(x = x, y = y, fill = Suitability)
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    name = "Climatic\nsuitability"
  ) +
  geom_contour(
    data = contour_fig11_df,
    aes(x = x, y = y, z = z),
    breaks = 0,
    colour = "red",
    linewidth = 0.40
  ) +
  geom_sf(
    data = europe,
    fill = NA,
    color = "grey65",
    linewidth = 0.18
  ) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(0.12, "in"),
    pad_y = grid::unit(0.12, "in"),
    height = grid::unit(0.34, "in"),
    width = grid::unit(0.27, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  facet_grid(Threshold ~ Month) +
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
    subtitle = "Red contour indicates the monthly climatological Tmin threshold",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 13) +
  theme(
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
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 2, b = 5)
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 11.5,
      margin = margin(l = 5, r = 5)
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 9
    ),
    legend.title = element_text(
      face = "bold",
      size = 11.5
    ),
    legend.text = element_text(
      size = 10.5
    ),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.55, "cm"),
    panel.spacing.y = grid::unit(0.55, "cm"),
    plot.margin = margin(12, 15, 12, 12)
  )

save_fig(
  p_fig11,
  "Figure_11_Spring_temperature_overlay_2x3.png",
  folder = "MAIN",
  width = 13,
  height = 9,
  dpi = 450
)

# Figure 12: GLM Response curves
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

glm_models_fig12 <- grep(
  "_GLM$",
  biomod2::get_built_models(main_run$model),
  value = TRUE
)

glm_models_fig12 <- glm_models_fig12[
  !grepl("allRun", glm_models_fig12)
]

stopifnot(length(glm_models_fig12) == 20)

rc_glm_fig12 <- biomod2::bm_PlotResponseCurves(
  bm.out = main_run$model,
  models.chosen = glm_models_fig12,
  fixed.var = "median",
  show.variables = final_vars_fig12
)

df_glm_fig12 <- rc_glm_fig12$tab

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

save_fig(
  p_fig12,
  "Figure_12_Response_curves_GLM.png",
  folder = "MAIN",
  width = 10,
  height = 6.5,
  dpi = 450
)

# ==============================================================================
# 14. SUPPLEMENTARY TABLES AND SESSION INFO
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

cat("\n")
cat("=================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("=================================================================\n")
cat("Project root:          ", project_root, "\n")
cat("Outputs:               ", output_root, "\n")
cat("Final presence records:", n_final, "\n")
cat("Main predictors:       ", paste(final_keep, collapse = ", "), "\n")
cat("Background sampling:   uniform random (10,000 cells)\n")
cat("Sensitivity analysis:  occurrence-density weighted (10,000 cells)\n")
cat("WorldClim version:     2.1 (corrected °C scaling)\n")
cat("=================================================================\n")
