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
#'               must have: ns, lu.to, times, value
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
      # LIGHT GRID BETWEEN PANELS (maps structure)
      # -----------------------
      panel.grid.major = ggplot2::element_line(
        color = "grey95",
        linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),

      # no heavy frame
      panel.border = ggplot2::element_blank(),

      # spacing between LU × years
      panel.spacing = ggplot2::unit(0.8, "lines"),

      # -----------------------
      # FACET STRIPS
      # -----------------------
      strip.background = ggplot2::element_blank(),

      strip.text.x = ggplot2::element_text(
        face = "bold",
        color = "#2C3E50",
        size = 11
      ),

      strip.text.y = ggplot2::element_text(
        face = "bold",
        color = "#1B4F72",
        size = 11,
        angle = 90
      ),

      strip.placement = "outside",

      # -----------------------
      # AXES (clean map style)
      # -----------------------
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),

      # -----------------------
      # LEGEND
      # -----------------------
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )
}




fdr_plot_downscaled_maps <- function(
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

  # -----------------------------
  # Convert ids to raster index
  # -----------------------------
  out_int <- fdr_to_ns_int(out_res, ns_map)

  df_pix <- terra::as.data.frame(rasterized_layer, xy = TRUE, na.rm = FALSE)
  names(df_pix)[3] <- "ns"
  df_pix <- dplyr::filter(df_pix, !is.na(ns))

  # -----------------------------
  # Aggregate + filter
  # -----------------------------
  inputs <- out_int %>%
    dplyr::group_by(ns, lu.to, times) %>%
    dplyr::summarise(value = sum(value), .groups = "drop")

  if (!is.null(LU)) {
    inputs <- inputs %>% dplyr::filter(lu.to == LU)
  }
  if (!is.null(year)) {
    inputs <- inputs %>% dplyr::filter(times == year)
  }

  plot_df <- df_pix %>%
    dplyr::left_join(inputs, by = "ns") %>%
    dplyr::filter(!is.na(lu.to), !is.na(times))

  lu_order <- c("cropland",  "newforest","otherland",  "pasture", "forest")

  plot_df$lu.to <- factor(plot_df$lu.to, levels = lu_order)

  if (is.null(limits)) {
    limits <- range(plot_df$value, na.rm = TRUE)
  }

  # -----------------------------
  # Build plot with ggnewscale
  # -----------------------------
  library(ggnewscale)

  p <- ggplot2::ggplot()

  # ---- CROPLAND
  p <- p +
    ggplot2::geom_raster(
      data = dplyr::filter(plot_df, lu.to == "cropland"),
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#B8860B",
      limits = limits,
      name = "Cropland (1000 ha)"
    )

  p <- p + ggnewscale::new_scale_fill()

  # ---- FOREST
  p <- p +
    ggplot2::geom_raster(
      data = dplyr::filter(plot_df, lu.to == "forest"),
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#006400",
      limits = limits,
      name = "Forest (1000 ha)"
    )

  p <- p + ggnewscale::new_scale_fill()

  # ---- NEW FOREST
  p <- p +
    ggplot2::geom_raster(
      data = dplyr::filter(plot_df, lu.to == "newforest"),
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#90EE90",
      limits = limits,
      name = "New forest (1000 ha)"
    )

  p <- p + ggnewscale::new_scale_fill()

  # ---- OTHER LAND
  p <- p +
    ggplot2::geom_raster(
      data = dplyr::filter(plot_df, lu.to == "otherland"),
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#6A0DAD",
      limits = limits,
      name = "Other land (1000 ha)"
    )

  p <- p + ggnewscale::new_scale_fill()

  # ---- PASTURE
  p <- p +
    ggplot2::geom_raster(
      data = dplyr::filter(plot_df, lu.to == "pasture"),
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#B22222",
      limits = limits,
      name = "Pasture (1000 ha)"
    )

  # -----------------------------
  # Layout
  # -----------------------------
  p <- p +
    ggplot2::coord_equal(expand = FALSE) +
    theme_fdr_map()+
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::facet_grid(
      times ~ lu.to,
      labeller = ggplot2::labeller(lu.to = c(
        cropland = "Cropland",
        forest = "Forest",
        newforest = "New forest",
        otherland = "Other land",
        pasture = "Pasture"
      ))
    )

  # -----------------------------
  # Border
  # -----------------------------
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
