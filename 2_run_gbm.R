# =============================================================
# Predicting undiscovered non‑human‑primate hosts of
# Semliki Forest complex Alphaviruses
# Author: Michael Celone   •   Updated: 2025‑06‑17
# =============================================================
#  ──  Pipeline overview  ─────────────────────────────────────
#    1) Parallel grid–search to identify optimal hyper‑parameters
#    2) Bootstrap‑wrapped GBM modelling with full predictor set
#    3) Variable‑importance filtering (> 1 % mean rel.inf)
#    4) Reduced‑predictor GBM bootstrap run
#    5) Label‑permuted null model → corrected AUC
# =============================================================

# ── 0. Setup ─────────────────────────────────────────────────

# Install missing pkgs on the fly (optional)
req_pkgs <- c(
  "gbm", "pROC", "caret", "pdp", "dplyr", "tidyr", "foreach", "progressr",
  "doParallel", "here", "Metrics", "readr", "arrow" # arrow → Parquet
)
new_pkgs <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

# Load libraries -----------------------------------------------------------
lapply(req_pkgs, library, character.only = TRUE)
options(dplyr.summarise.inform = FALSE)
handlers(global = TRUE)
handlers("txtprogressbar")

# ── 1. Paths & global parameters ─────────────────────────────────────────
project_dir  <- here::here()                                   # root dir
output_dir   <- file.path(project_dir, sprintf("NHP_model_run_%s", Sys.Date()))
if (!dir.exists(output_dir)) dir.create(output_dir)

data_path    <- file.path(project_dir, "data", "prim_data_final.csv")

# Seed for reproducibility (global)
SEED <- 123
set.seed(SEED)

# ── 2. Read data & define metadata  ──────────────────────────────────────
prim_data <- readr::read_csv(data_path, show_col_types = FALSE)

label   <- "vir_pos"                         # response column
species <- "species"                         # species column

# Exclude ID cols here explicitly (safer than numeric positions)
id_cols <- c(species, label)
vars    <- setdiff(names(prim_data), id_cols) # predictor names

# ── 3. Parallel back‑end --------------------------------------------------
N_CORES <- max(1, parallel::detectCores() - 1) # leave 1 core free
cl <- parallel::makeCluster(N_CORES)
doParallel::registerDoParallel(cl)
# (stopCluster(cl) is called in on.exit)
on.exit(parallel::stopCluster(cl), add = TRUE)

# ── 4. Helper functions  ─────────────────────────────────────────────────

# 4.1  Grid‑search GBM -----------------------------------------------------
run_grid_search <- function(train_df, grid, n_trees = 50000, seed = SEED) {
  foreach(i = seq_len(nrow(grid)), .combine = dplyr::bind_rows, .packages = "gbm") %dopar% {
    set.seed(seed + i)  # vary seed slightly across workers
    fit <- gbm(
      formula            = vir_pos ~ . -species -wos_hits,
      distribution       = "bernoulli",
      data               = train_df,
      n.trees            = n_trees,
      interaction.depth  = grid$interaction.depth[i],
      shrinkage          = grid$shrinkage[i],
      n.minobsinnode     = grid$n.minobsinnode[i],
      bag.fraction       = 0.5,
      train.fraction     = 0.75,
      verbose            = FALSE
    )
    tibble::tibble(
      shrinkage          = grid$shrinkage[i],
      interaction.depth  = grid$interaction.depth[i],
      n.minobsinnode     = grid$n.minobsinnode[i],
      optimal_trees      = which.min(fit$valid.error),
      min_RMSE           = sqrt(min(fit$valid.error))
    )
  }
}

# 4.2  Bootstrap GBM runner -----------------------------------------------
run_bootstrap_gbm <- function(df, params, label_col, species_col,
                              n_runs = 10, n_trees = 50000,
                              cv_folds = 5, seed = SEED) {
  foreach(m = seq_len(n_runs), .packages = c("gbm", "caret", "dplyr")) %dopar% {
    set.seed(seed + m)
    idx   <- caret::createDataPartition(df[[label_col]], p = 0.8)[[1]]
    train <- df[idx, ]
    test  <- df[-idx, ]
    
    fit <- gbm(
      formula            = as.formula(paste(label_col, "~ . -", species_col, "-wos_hits")),
      distribution       = "bernoulli",
      data               = train,
      n.trees            = n_trees,
      interaction.depth  = params$interaction.depth,
      shrinkage          = params$shrinkage,
      n.minobsinnode     = params$n.minobsinnode,
      cv.folds           = cv_folds,
      bag.fraction       = 0.5,
      verbose            = FALSE
    )
    
    best_iter <- gbm.perf(fit, method = "cv", plot.it = FALSE)
    fit$var.levels <- lapply(fit$var.levels, function(x) replace(x, is.infinite(x), 0))
    
    pred_train <- predict(fit, newdata = train, n.trees = best_iter, type = "response")
    pred_test  <- predict(fit, newdata = test,  n.trees = best_iter, type = "response")
    
    auc_train <- pROC::auc(train[[label_col]], pred_train)[1]
    auc_test  <- pROC::auc(test[[label_col]],  pred_test)[1]
    
    rmse_train <- Metrics::rmse(train[[label_col]], pred_train)
    rmse_test  <- Metrics::rmse(test[[label_col]],  pred_test)
    
    vi <- tibble::tibble(var = fit$var.names,
                         rel.inf = 0) |>
      dplyr::left_join(as.data.frame(summary(fit, plotit = FALSE)), by = c("var" = "var")) |>
      dplyr::mutate(rel.inf = dplyr::coalesce(rel.inf.y, rel.inf.x)) |>
      dplyr::select(var, rel.inf)
    
    list(
      fit          = fit,
      best_iter    = best_iter,
      metrics      = tibble::tibble(run = m, auc_train, auc_test, rmse_train, rmse_test),
      var_imp      = vi,
      pred_all     = tibble::tibble(
        species        = df[[species_col]],
        prediction     = predict(fit, newdata = df, n.trees = best_iter, type = "response"),
        bootstrap_run  = m,
        observed       = df[[label_col]]
      )
    )
  }
}


# 4.3  Helper to aggregate bootstrap results ------------------------------
aggr_bootstrap <- function(boot_list) {
  metrics   <- dplyr::bind_rows(lapply(boot_list, `[[`, "metrics"))
  var_imp   <- dplyr::bind_rows(lapply(boot_list, `[[`, "var_imp")) |>
    dplyr::group_by(var) |>
    dplyr::summarise(mean_rel.inf = mean(rel.inf, na.rm = TRUE), .groups = "drop")
  preds_all <- dplyr::bind_rows(lapply(boot_list, `[[`, "pred_all"))
  list(metrics = metrics, var_imp = var_imp, preds_all = preds_all)
}

# 4.4  Null‑model runner ---------------------------------------------------
run_null_model <- function(df, params, n_runs = 10, n_trees = 50000, seed = SEED,
                           label_col, species_col) {
  foreach(m = seq_len(n_runs), .combine = c, .packages = c("gbm", "caret", "pROC")) %dopar% {
    set.seed(seed + m)
    idx   <- caret::createDataPartition(df[[label]], p = 0.8)[[1]]
    train <- df[idx, ]
    test  <- df[-idx, ]
    # Shuffle labels
    train[[label]] <- sample(train[[label]])
    test[[label]]  <- sample(test[[label]])
    
    fit <- gbm(
      formula            = vir_pos ~ . -species -wos_hits,
      distribution       = "bernoulli",
      data               = train,
      n.trees            = n_trees,
      interaction.depth  = params$interaction.depth,
      shrinkage          = params$shrinkage,
      n.minobsinnode     = params$n.minobsinnode,
      cv.folds           = 5,
      bag.fraction       = 0.5,
      verbose            = FALSE
    )
    best_iter <- gbm.perf(fit, method = "cv", plot.it = FALSE)
    pred_test <- predict(fit, newdata = test, n.trees = best_iter, type = "response")
    as.numeric(pROC::auc(test[[label]], pred_test)[1])
  }
}

# ── 5. Grid search  ─────────────────────────────────────────────────────-

grid <- expand.grid(
  shrinkage          = c(1e-4, 1e-3, 1e-2),
  interaction.depth  = 2:4,
  n.minobsinnode     = 2:5
)

# Train/test split for tuning only once -----------------------------------
tune_idx <- caret::createDataPartition(prim_data[[label]], p = 0.8)[[1]]
train_tune <- prim_data[tune_idx, ]

grid_res  <- run_grid_search(train_tune, grid)
readr::write_csv(grid_res, file.path(output_dir, "optimal_parameters.csv"))

best_row <- grid_res |> dplyr::slice_min(min_RMSE, n = 1)
params   <- best_row |> dplyr::select(shrinkage, interaction.depth, n.minobsinnode)

# ── 6. Full‑predictor bootstrap models ─────────────────────────────────––
cat(sprintf("\nRunning full‑predictor models with η=%.4f, depth=%d, n.minobs=%d\n",
            params$shrinkage, params$interaction.depth, params$n.minobsinnode))

boot_full <- run_bootstrap_gbm(prim_data, params, n_runs = 5, label_col = label,
                               species_col = species)
full_aggr <- aggr_bootstrap(boot_full)

# Save outputs -------------------------------------------------------------
write_results <- function(df, name) readr::write_csv(df, file.path(output_dir, name))

write_results(full_aggr$metrics,   "model_perf_all_vars.csv")
write_results(full_aggr$var_imp,   "var_importance_all_vars.csv")
write_results(full_aggr$preds_all,  "predictions_all_vars.csv")

write_results(grid_res,             "grid_results.csv")
write_results(best_row,             "best_hyperparams.csv")

# Also save as RDS for faster load
saveRDS(list(models = boot_full, summary = full_aggr), file.path(output_dir, "full_model_bootstrap.rds"))

# ── 7. Variable‑reduced model (>1 % mean rel.inf)  ───────────────────────
keep_vars <- full_aggr$var_imp |> dplyr::filter(mean_rel.inf > 1) |> dplyr::pull(var)
cat(sprintf("\nRetaining %d predictors with >1%% importance (of %d total)\n",
            length(keep_vars), length(vars)))

prim_data_red <- prim_data |> dplyr::select(all_of(c(id_cols, keep_vars)))

boot_red  <- run_bootstrap_gbm(prim_data_red, params, n_runs = 10, label_col = label,
                               species_col = species)
red_aggr  <- aggr_bootstrap(boot_red)

write_results(red_aggr$metrics,   "model_perf_reduced.csv")
write_results(red_aggr$var_imp,   "var_importance_reduced.csv")
write_results(red_aggr$preds_all, "predictions_reduced.csv")

saveRDS(list(models = boot_red, summary = red_aggr), file.path(output_dir, "reduced_model_bootstrap.rds"))

# ── 8. Null model & corrected AUC ─────────────────────────────────────────
cat("\nRunning null model …\n")
null_aucs <- run_null_model(prim_data, params, n_runs = 10, label_col = label,
                            species_col = species)

corrected_auc <- mean(red_aggr$metrics$auc_test) - (mean(null_aucs) - 0.5)
cat(sprintf("\nCorrected AUC  =  %.3f\n", corrected_auc))

readr::write_csv(tibble::tibble(corrected_auc = corrected_auc),
                 file.path(output_dir, "corrected_auc.csv"))

# ── 9. Done ───────────────────────────────────────────────────────────────
cat(sprintf("\nOutputs saved to: %s\n", output_dir))

# Clean up
parallel::stopCluster(cl)
