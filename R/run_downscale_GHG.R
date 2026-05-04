# -----------------------------------------------------------------------------
# fdr_run_downscaling()
#
# Purpose:
# - Estimate land-use transition "preferences" from observed historical changes
#   (calibration period) using multinomial logit (MNL).
# - Use these estimated coefficients ("betas") as priors in downscalr::downscale()
#   to spatially distribute national (FABLE) land-use change targets.
#
# Key idea:
# - We fit ONE MNL per origin land-use (lu.from). Example:
#     cropland -> {cropland, pasture, urban, ...}
# - For each origin, the dependent variable Y is the *share* of the origin land
#   that transitions to each destination land-use in each grid cell.
# - Explanatory variables X are the covariates (priors) at the grid-cell level.
#
# Why so many checks?
# - Spatial calibration data often include:
#     * missing covariates (NAs),
#     * destinations that are always zero (no information),
#     * near-zero destinations (quasi-separation / numerical issues),
#     * collinearity in covariates (X) => singular matrix in MNL internals,
#     * cells where nothing happens for an origin (row sums = 0),
#   which can make MNL fail or produce unstable priors.
# -----------------------------------------------------------------------------

#' Fit MNL models and run DownscalR downscaling
#'
#' @param targets FABLE transition targets (long: times, lu.from, lu.to, value)
#' @param country_luc Calibration transitions (output of lc_build_country_luc)
#' @param priors Output of fdr_build_priors()
#' @param mnl_niter Number of MCMC iterations
#' @param mnl_nburn Burn-in iterations
#'
#' @return List with:
#'   - downscaled_LUC
#'   - X_long
#'   - pred_coeff_long
#'   - country_start_areas
#'
#' @export

fdr_run_downscaling_GHG <- function(
    targets,
    country_luc,
    priors,
    mnl_niter = 100,
    mnl_nburn = 50
) {

  if (!requireNamespace("downscalr", quietly = TRUE)) {
    stop("Package 'downscalr' is required.")
  }

  # ---------------------------------------------------------------------------
  # Helper: stabilize X (covariates) so MNL is numerically feasible
  #
  # Why:
  # - MNL inside downscalr involves solving linear systems.
  # - Collinearity / near-collinearity (high correlation among covariates)
  #   can make matrices close to singular → solve() fails.
  # - Near-constant columns are also problematic (no information, can create
  #   numerical degeneracy after scaling).
  #
  # What we do:
  # 1) Drop near-constant columns (sd ~ 0)
  # 2) Optionally scale X (centering + dividing by sd) to improve conditioning
  # 3) Drop (near) linear dependencies using QR pivoting (keeps full-rank set)
  #
  # Returns:
  # - list(X = cleaned matrix, keep_cols = names of columns kept)
  # ---------------------------------------------------------------------------
  fdr_make_full_rank_X <- function(X, tol_sd = 1e-10, tol_qr = 1e-6, scale_X = TRUE) {
    stopifnot(is.matrix(X))

    # (1) Drop near-constant columns
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep_sd <- which(is.finite(sds) & sds > tol_sd)
    X <- X[, keep_sd, drop = FALSE]

    if (ncol(X) < 2) stop("X has <2 usable columns after sd filtering.")

    # (2) Scale to improve numerical conditioning (optional but usually helpful)
    if (scale_X) {
      mu  <- colMeans(X, na.rm = TRUE)
      sdv <- apply(X, 2, sd, na.rm = TRUE)
      sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
      X <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    }

    # (3) QR pivoting to drop collinear columns and keep full-rank subset
    q <- qr(X, tol = tol_qr, LAPACK = TRUE)
    keep_qr <- q$pivot[seq_len(q$rank)]
    X <- X[, keep_qr, drop = FALSE]

    list(X = X, keep_cols = colnames(X))
  }

  # ---------------------------------------------------------------------------
  # 1) Prepare X (covariate matrix) in wide format for MNL
  #
  # priors$xmat is long format: ns, ks, value
  # MNL expects X as a numeric matrix: rows=ns (cells), cols=covariates (ks)
  # ---------------------------------------------------------------------------
  X_long    <- priors$xmat
  lu_levels <- priors$lu_levels

  X_wide <- X_long %>%
    tidyr::pivot_wider(
      names_from  = "ks",
      values_from = "value",
      values_fill = 0
    ) %>%
    dplyr::arrange(ns)

  ns_ids <- X_wide$ns

  X_mat <- as.matrix(X_wide[, -1])
  rownames(X_mat) <- ns_ids

  # Global cleaning ONCE:
  # Why once globally?
  # - Keeps the same covariate set across origins (easier to reason about),
  # - avoids origin-specific covariate subsets that complicate later use of betas.
  clean <- fdr_make_full_rank_X(
    X_mat,
    tol_sd  = 1e-10,
    tol_qr  = 1e-6,
    scale_X = TRUE
  )
  X_mat   <- clean$X
  keep_ks <- clean$keep_cols

  # Keep long xmat consistent with the cleaned wide matrix
  # (DownscalR downscale() later expects the ks in xmat to match betas)
  X_long <- X_long %>%
    dplyr::filter(ks %in% keep_ks)

  # ---------------------------------------------------------------------------
  # 2) Prepare calibration transitions (Yraw)
  #
  # country_luc is long transitions: ns, lu.from, lu.to, Ts, value (shares)
  # We restrict to ns that exist in X_mat so we can align rows.
  # ---------------------------------------------------------------------------
  Yraw <- country_luc %>%
    dplyr::mutate(ns = as.character(ns)) %>%
    dplyr::filter(ns %in% rownames(X_mat))

  origins <- sort(unique(Yraw$lu.from))

  pred_coeff_long     <- tibble::tibble()
  country_start_areas <- lu_levels %>%
    dplyr::mutate(times = min(targets$times))


  # ---------------------------------------------------------------------------
  # 3) Fit ONE MNL per origin land-use
  #
  # For each lu_from:
  # - Build Y matrix: rows=ns, cols=destinations (lu.to), entries=shares
  # - Filter destinations with no/too little information
  # - Filter rows with no transition mass
  # - Fit MNL with ridge fallback (A0 grid)
  # - Store posterior mean coefficients as (lu.from, lu.to, ks, value)
  # ---------------------------------------------------------------------------
  for (lu_from in origins) {

    # --- Build destination share table for this origin ---
    # Why:
    # - MNL models a multinomial distribution over destinations
    # - We want a consistent Y matrix per origin: ns x destinations
    Y_df <- Yraw %>%
      dplyr::filter(lu.from == lu_from, Ts == 2015) %>%
      dplyr::select(ns, lu.to, value) %>%
      dplyr::mutate(
        ns    = as.character(ns),
        lu.to = as.character(lu.to),
        value = as.numeric(value)
      ) %>%
      tidyr::pivot_wider(
        names_from  = lu.to,
        values_from = value,
        values_fill = 0
      ) %>%
      dplyr::arrange(ns)

    # Convert to numeric matrix (drop ns key column)
    Y <- as.matrix(Y_df[, setdiff(names(Y_df), "ns"), drop = FALSE])
    storage.mode(Y) <- "double"
    Y[is.na(Y)] <- 0

    # -----------------------------------------------------------------------
    # Drop destination columns with no information or extremely rare mass
    #
    # Why:
    # - If a destination is always 0 in calibration, it provides no signal and
    #   can destabilize MNL routines.
    # - Extremely rare destinations can cause quasi-separation / numerical
    #   singularities (common in multinomial settings).
    #
    # Example in our case:
    # - "forest" may be structurally zero because we disallow transitions to forest.
    # -----------------------------------------------------------------------
    col_totals <- colSums(Y)

    keep_dest <- which(col_totals > 0)

    min_share <- 1e-4  # tune if needed (e.g., 1e-3 for more aggressive filtering)
    keep_dest <- keep_dest[col_totals[keep_dest] / sum(col_totals) > min_share]

    Y <- Y[, keep_dest, drop = FALSE]

    # Need at least 2 destinations for a multinomial model
    if (ncol(Y) < 2) {
      message("Skipping origin ", lu_from,
              " (not enough informative destinations after filtering) - DownscalR will later use flat/zero priors for this transition")
      next
    }

    # -----------------------------------------------------------------------
    # Align rows (ns) between Y and X
    #
    # Why:
    # - X_mat is built from priors; Y_df is built from calibration transitions.
    # - Some ns may be present in Y_df but removed earlier from priors (NAs).
    # - We keep only overlapping ns and ensure SAME ORDER.
    # -----------------------------------------------------------------------
    idx <- match(Y_df$ns, rownames(X_mat))
    keep_rows <- !is.na(idx)

    if (!any(keep_rows)) {
      message("Skipping origin ", lu_from, " (no overlapping ns with X)")
      next
    }

    Y <- Y[keep_rows, , drop = FALSE]
    X_use <- X_mat[idx[keep_rows], , drop = FALSE]

    stopifnot(nrow(Y) == nrow(X_use))

    # -----------------------------------------------------------------------
    # Drop rows where nothing happens (row sum = 0)
    #
    # Why:
    # - After dropping rare destinations, some cells may end up with all zeros.
    # - These rows contain no information for estimation and break normalization.
    # -----------------------------------------------------------------------
    rs <- rowSums(Y)
    good_rows <- rs > 0

    if (sum(good_rows) < 2) {
      message("Skipping origin ", lu_from,
              " (too few non-zero rows after filtering)")
      next
    }

    Y     <- Y[good_rows, , drop = FALSE]
    X_use <- X_use[good_rows, , drop = FALSE]

    # -----------------------------------------------------------------------
    # Renormalize rows to sum to 1 (shares)
    #
    # Why:
    # - Some workflows store Y as areas or partial shares.
    # - MNL expects each row to be a probability vector over destinations.
    # - Even if inputs are shares already, filtering columns/rows changes totals.
    # -----------------------------------------------------------------------
    Y <- Y / rowSums(Y)

    # Final sanity checks before estimation
    stopifnot(is.matrix(X_use), is.numeric(X_use))
    stopifnot(is.matrix(Y),     is.numeric(Y))
    stopifnot(!anyNA(X_use), !anyNA(Y))
    stopifnot(nrow(X_use) == nrow(Y))

    # Diagnostic printed to understand why n varies across origins
    message(
      "Origin=", lu_from,
      " | n=", nrow(X_use),                  # number of cells kept after filtering
      " p=", ncol(X_use),                    # number of covariates used
      " kappa=", signif(kappa(X_use), 3),    # conditioning diagnostic
      " | dest=", ncol(Y)                    # number of destination classes
    )

    # -----------------------------------------------------------------------
    # Fit MNL with ridge fallback
    #
    # Why:
    # - Even after cleaning X, some origins remain hard (e.g., cropland),
    #   causing singular systems inside mnlogit.
    # - Ridge (A0) stabilizes the posterior by shrinking coefficients.
    # - We try increasing ridge until it works.
    # -----------------------------------------------------------------------
    ridge_grid <- c(100, 500, 1000, 5000, 10000, 50000)

    fit <- NULL
    for (A0 in ridge_grid) {
      fit_try <- try(
        downscalr::mnlogit(
          X        = X_use,
          Y        = Y,
          baseline = which.max(colSums(Y)),
          niter    = mnl_niter,
          nburn    = mnl_nburn,
          A0       = A0
        ),
        silent = TRUE
      )

      if (!inherits(fit_try, "try-error")) {
        fit <- fit_try
        break
      } else {
        msg <- conditionMessage(attr(fit_try, "condition"))
        message("mnlogit failed for ", lu_from, " with A0=", A0, " : ", msg)
      }
    }

    if (is.null(fit)) {
      message("Skipping origin ", lu_from,
              " (mnlogit failed even with ridge; priors will be missing for this origin)")
      next
    }
    # -----------------------------------------------------------------------
    # Extract posterior mean betas
    #
    # Why:
    # - downscalr::downscale() needs betas in long format:
    #     lu.from, lu.to, ks, value
    # - We use posterior means as stable priors for the downscaling step.
    # -----------------------------------------------------------------------
    post <- fit$postb
    Bmean <- apply(post, c(1, 2), mean)

    coef_long <- as.data.frame(Bmean) %>%
      tibble::rownames_to_column("ks") %>%
      tidyr::pivot_longer(-ks, names_to = "lu.to", values_to = "value") %>%
      dplyr::mutate(lu.from = lu_from, .before = 1)

    pred_coeff_long <- dplyr::bind_rows(pred_coeff_long, coef_long)

  }
  message("n = number of spatial cells that have usable transition information for this origin land-use class after all filtering steps.")

  # --- Ensure we have betas for every target pair required by downscale() ---

  ks_all <- sort(unique(as.character(X_long$ks)))

  target_pairs <- targets %>%
    dplyr::transmute(
      lu.from = as.character(lu.from),
      lu.to   = as.character(lu.to)
    ) %>%
    dplyr::distinct()

  existing_beta_keys <- pred_coeff_long %>%
    dplyr::transmute(
      lu.from = as.character(lu.from),
      lu.to   = as.character(lu.to),
      ks      = as.character(ks)
    ) %>%
    dplyr::distinct()

  required_beta_keys <- tidyr::crossing(
    target_pairs,
    ks = ks_all
  )

  missing_beta_keys <- dplyr::anti_join(
    required_beta_keys,
    existing_beta_keys,
    by = c("lu.from", "lu.to", "ks")
  )

  if (nrow(missing_beta_keys) > 0) {
    message(
      "Completing missing beta coefficients with 0 for ",
      nrow(missing_beta_keys), " (lu.from, lu.to, ks) combinations."
    )

    filler <- missing_beta_keys %>%
      dplyr::mutate(value = 0)

    pred_coeff_long <- dplyr::bind_rows(pred_coeff_long, filler)
  }

  pred_coeff_long <- pred_coeff_long %>%
    dplyr::mutate(
      lu.from = as.character(lu.from),
      lu.to   = as.character(lu.to),
      ks      = as.character(ks),
      value   = as.numeric(value)
    ) %>%
    dplyr::summarise(
      value = dplyr::first(value),
      .by = c(lu.from, lu.to, ks)
    )

  # ---------------------------------------------------------------------------
  # 4) Run DownscalR spatial allocation
  #
  # Why:
  # - targets: national (or large-scale) transition constraints from FABLE
  # - start.areas: baseline per-cell land use
  # - xmat + betas: priors controlling where transitions happen
  # ---------------------------------------------------------------------------
  results <- downscalr::downscale(
    targets     = targets,
    start.areas = country_start_areas,
    xmat        = X_long,
    betas       = pred_coeff_long
  )


  # ---------------------------------------------------------------------------
  # 5) Computing GHG due to Land use change by ecoregion
  # ---------------------------------------------------------------------------
  library(readxl)
  library(readr)
  EF_Pools_transition_Ecoregion <- read_csv("Data/EF_Pools_transition_Ecoregion.csv")

  results$out.res <- results$out.res %>%
    left_join(EF_Pools_transition_Ecoregion %>%
                select(id_c, from, to, ef_biomass),
              by = c("ns" ="id_c", "lu.from" = "from", "lu.to" = "to")
    ) %>%
    mutate(GHG_biomass = ef_biomass * value)

  # ---------------------------------------------------------------------------
  # 6) Return clean outputs to the workflow
  # ---------------------------------------------------------------------------

   list(
    out.res            = results$out.res,
    downscaled_LUC      = results$out.res,
    X_long              = X_long,
    pred_coeff_long     = pred_coeff_long,
    country_start_areas = country_start_areas
  )





}
