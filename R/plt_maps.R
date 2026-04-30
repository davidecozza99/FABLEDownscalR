# -----------------------------------------------------------------------------
# fdr_to_ns_int
#' Convert a results table keyed by id_c into ns_int for raster plotting
#'
#' @param df long results with column ns (character id_c)
#' @param ns_map data.frame with columns id_c and ns_int
#' @return same df but ns replaced by ns_int
#' @export
# -----------------------------------------------------------------------------
fdr_to_ns_int <- function(df, ns_map) {
  chk_required_cols(df, c("ns"))
  chk_required_cols(ns_map, c("id_c", "ns_int"))

  df %>%
    dplyr::mutate(ns = as.character(ns)) %>%
    dplyr::left_join(
      ns_map %>% dplyr::mutate(id_c = as.character(id_c)),
      by = c("ns" = "id_c")
    ) %>%
    dplyr::mutate(ns = as.integer(ns_int)) %>%
    dplyr::select(-ns_int)
}

# -----------------------------------------------------------------------------
#' Plot downscaled results on an ID raster (facet by time and origins)
#'
#' @param out_res results table (typically fdr_run_downscaling()$out.res or $downscaled_LUC)
#' must have: ns, lu.to, times, value
#' @param rasterized_layer SpatRaster with ns_int values (from fdr_build_id_maps())
#' @param ns_map id_c -> ns_int mapping (from fdr_build_ns_map())
#' @param limits numeric vector length 2 for consistent color scales across facets
#' @param palette distiller palette name
#' @param na_color fill color for NA pixels
#' @export
# -----------------------------------------------------------------------------
theme_fdr_map <- function(base_size = 11) {

  ggplot2::theme_minimal(base_size = base_size) +

    ggplot2::theme(

      # -----------------------
      # Background
      # -----------------------
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),

      # -----------------------
      # MATRIX GRID
      # -----------------------
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "grey90",   # light grid lines
        linewidth = 0.4
      ),

      panel.spacing = ggplot2::unit(0, "lines"),

      # -----------------------
      # Facet strips
      # -----------------------
      strip.background = ggplot2::element_blank(),

      strip.text.x = ggplot2::element_text(
        face = "bold",
        color = "black",
        size = 11
      ),

      strip.text.y = ggplot2::element_text(
        face = "bold",
        color = "black",
        size = 11,
        angle = 270
      ),

      strip.placement = "outside",

      # -----------------------
      # Map style
      # -----------------------
      panel.grid = ggplot2::element_blank(),

      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),

      # -----------------------
      # Legend
      # -----------------------
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )
}

# LAND USE

fdr_plot_downscaled_LU <- function(
    out_res,
    rasterized_layer,
    ns_map,
    year = NULL,
    LU = NULL,
    limits = NULL,
    na_color = "grey90",
    add_border = TRUE
) {

  chk_required_cols(out_res, c("ns", "lu.to", "times", "value"))

  out_int <- fdr_to_ns_int(out_res, ns_map)

  df_pix <- terra::as.data.frame(rasterized_layer, xy = TRUE, na.rm = FALSE)
  names(df_pix)[3] <- "ns"
  df_pix <- dplyr::filter(df_pix, !is.na(ns))

  inputs <- out_int %>%
    dplyr::group_by(ns, lu.to, times) %>%
    dplyr::summarise(value = sum(value), .groups = "drop")

  if (!is.null(LU)) {
    inputs <- inputs %>% dplyr::filter(lu.to %in% LU)
  }

  if (!is.null(year)) {
    inputs <- inputs %>% dplyr::filter(times == year)
  }

  plot_df <- df_pix %>%
    dplyr::left_join(inputs, by = "ns") %>%
    dplyr::filter(!is.na(lu.to), !is.na(times))

  lu_order <- c("newforest", "cropland", "otherland", "forest", "pasture")
  plot_df$lu.to <- factor(plot_df$lu.to, levels = lu_order)

  lu_present <- na.omit(unique(as.character(plot_df$lu.to)))

  if (is.null(limits)) {
    limits <- range(plot_df$value, na.rm = TRUE)
  }

  lu_colors <- list(
    cropland = "#B8860B",
    forest = "#006400",
    newforest = "#90EE90",
    otherland = "#6A0DAD",
    pasture = "#B22222"
  )

  lu_labels <- c(
    cropland = "Cropland",
    forest = "Forest",
    newforest = "New forest",
    otherland = "Other land",
    pasture = "Pasture"
  )

  library(ggnewscale)

  p <- ggplot2::ggplot()

  for (i in seq_along(lu_present)) {

    lu <- lu_present[i]

    p <- p +
      ggplot2::geom_raster(
        data = dplyr::filter(plot_df, lu.to == lu),
        ggplot2::aes(x = x, y = y, fill = value)
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = lu_colors[[lu]],
        limits = limits,
        na.value = na_color,
        name = paste0(lu_labels[[lu]], " (1000 ha)")
      )

    if (i < length(lu_present)) {
      p <- p + ggnewscale::new_scale_fill()
    }
  }

  p <- p +
    ggplot2::coord_equal(expand = FALSE) +
    theme_fdr_map() +
    ggplot2::facet_grid(
      times ~ lu.to,
      labeller = ggplot2::labeller(lu.to = lu_labels)
    )

  if (add_border) {

    r <- rasterized_layer
    r[!is.na(r)] <- 1

    country_border <- terra::as.polygons(r, dissolve = TRUE)
    country_border <- sf::st_as_sf(country_border)

    p <- p +
      ggplot2::geom_sf(
        data = country_border,
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
  }

  return(p)
}




# LAND USE CHANGE

fdr_plot_downscaled_LUC <- function(
    out_res,
    rasterized_layer,
    ns_map,
    year = NULL,
    LU = NULL,
    limits = NULL,
    na_color = "grey90",
    add_border = TRUE
) {

  chk_required_cols(out_res, c("ns", "lu.to", "times", "value"))

  out_int <- fdr_to_ns_int(out_res, ns_map)

  df_pix <- terra::as.data.frame(rasterized_layer, xy = TRUE, na.rm = FALSE)
  names(df_pix)[3] <- "ns"
  df_pix <- dplyr::filter(df_pix, !is.na(ns))

  inputs <- out_int %>%
    dplyr::filter(lu.from != lu.to) %>%
    # Gains by destination class
    dplyr::group_by(lu.to, ns, times) %>%
    dplyr::summarise(gain = sum(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Type = "gain", lu = lu.to) %>%

    # Bind losses by origin class (negative sign for map interpretation)
    dplyr::bind_rows(
      out_int %>%
        dplyr::filter(lu.from != lu.to) %>%
        dplyr::group_by(lu.from, ns, times) %>%
        dplyr::summarise(loss = -sum(value, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(Type = "loss", lu = lu.from)
    ) %>%

    # Collapse gains and losses into one signed number per cell and land-use
    dplyr::group_by(lu, times, ns) %>%
    dplyr::summarise(value = sum(gain, loss, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(lu.to = lu)

  if (!is.null(LU)) {
    inputs <- inputs %>% dplyr::filter(lu.to %in% LU)
  }

  if (!is.null(year)) {
    inputs <- inputs %>% dplyr::filter(times == year)
  }

  plot_df <- df_pix %>%
    dplyr::left_join(inputs, by = "ns") %>%
    dplyr::filter(!is.na(lu.to), !is.na(times))

  lu_order <- c("cropland", "newforest", "otherland", "pasture", "forest", "urban")
  plot_df$lu.to <- factor(plot_df$lu.to, levels = lu_order)

  lu_labels <- c(
    cropland = "Cropland",
    forest = "Forest",
    newforest = "New forest",
    otherland = "Other land",
    pasture = "Pasture",
    urban = "Urban"
  )

  # ----------------------------
  # GLOBAL LIMITS
  # ----------------------------
  if (is.null(limits)) {
    max_abs <- max(abs(plot_df$value), na.rm = TRUE)
    limits <- c(-max_abs, max_abs)
  }

  # ----------------------------
  # PLOT
  # ----------------------------
  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_raster(
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b",      # dark red
      mid = "white",
      high = "#1a7f37",     # dark green
      midpoint = 0,
      limits = limits,
      na.value = na_color,
      name = "1000 ha"
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    theme_fdr_map() +
    ggplot2::facet_grid(
      times ~ lu.to,
      labeller = ggplot2::labeller(lu.to = lu_labels)
    )

  # ----------------------------
  # BORDER (same logic as your LU function)
  # ----------------------------
  if (add_border) {

    r <- rasterized_layer
    r[!is.na(r)] <- 1

    country_border <- terra::as.polygons(r, dissolve = TRUE)
    country_border <- sf::st_as_sf(country_border)

    p <- p +
      ggplot2::geom_sf(
        data = country_border,
        fill = NA,
        color = "black",
        linewidth = 0.5
      )
  }

  return(p)
}
