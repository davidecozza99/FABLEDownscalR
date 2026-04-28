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
fdr_plot_downscaled_maps <- function(
    out_res,
    rasterized_layer,
    ns_map,
    year=NULL, LU=NULL,
    limits = NULL,
    palette = "Greens",
    na_color = "grey90"
) {
  chk_required_cols(out_res, c("ns","lu.to","times","value"))

  # Convert id_c -> ns_int (so we can paint the raster)
  out_int <- fdr_to_ns_int(out_res, ns_map)

  # Pixel canvas from raster
  df_pix <- terra::as.data.frame(rasterized_layer, xy = TRUE, na.rm = FALSE)
  names(df_pix)[3] <- "ns"
  df_pix <- dplyr::filter(df_pix, !is.na(ns))

  ns = lu.to = times = value = x = y= NULL
  if(is.null(year) & is.null(LU)){
    inputs <- out_int %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep")
  } else if(!(is.null(year) | is.null(LU))){
    inputs <- out_int %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU & times==year)
  } else if(is.null(year)){
    inputs <- out_int %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU)
  } else {
    inputs <- out_int %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(times==year)
  }

  # Join pixels to values (many-to-many because we facet by lu.to and times)
  plot_df <- dplyr::left_join(df_pix, inputs, by = "ns", relationship = "many-to-many") %>%
    dplyr::filter(!is.na(lu.to), !is.na(times))

  # Stable global limits if not provided
  if (is.null(limits)) {
    limits <- range(plot_df$value, na.rm = TRUE)
  }

  ggplot2::ggplot(plot_df) +
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value, group = lu.to)) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "darkgreen",
      limits = limits,
      na.value = na_color
    ) +
    #ggplot2::labs(title = "")
    ggplot2::coord_equal(expand = FALSE) +
    ggthemes::theme_map() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::facet_grid(times ~ lu.to)
}
