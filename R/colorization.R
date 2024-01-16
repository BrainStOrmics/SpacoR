#' @rdname colorization
#' @title Colorize cell clusters based on spatial distribution
#' @name colorize
#' @description
#' Colorize cell clusters based on spatial distribution, so that spatially interlaced and spatially neighboring clusters are assigned with more perceptually different colors.
#' Spaco2 provides 3 basic color mapping mode:
#' 1. Optimize the mapping of a pre-defined color palette.
#' 2. Extract colors from image.
#' 3. Automatically generate colors within colorspace.
#' @param cell_coordinates a list like object containing spatial coordinates for each cell.
#' @param cell_labels a list like object containing cluster labels for each cell.
#' @param colorblind_type Optional parameter.
#' @param palette a list of colors (in hex). If given, `image_palette`will be ignored. See `Mode 1` above. Defaults to None.
#' @param image_palette an image in numpy array format. Should be a typical RGB image of shape (x, y, 3).
#' Ignored if `palette` is given. See `Mode 2`above. Defaults to None.
#' @param manual_mapping a data structure for manual color mapping including cluster names and manually assigned colors (in hex).
#' @param neighbor_weight Weight for calculating cell neighborhood.Defaults to 0.5.
#' @param radius radius used to calculate cell neighborhood.Defaults to 90.
#' @param n_neighbors k for KNN neighbor detection.Defaults to 16.
#' @param neighbor_args arguments passed to `spatial_distance` function.
#' @param mapping_args arguments passed to `map_graph` function.
#' @param embed_args arguments passed to `embed_graph` function.
#' @return Optimized color mapping for clusters, including cluster names and corresponding hex
#' @export
colorize <- function(
    cell_coordinates,
    cell_labels,
    colorblind_type = c("none", "protanopia", "deuteranopia", "tritanopia", "general"),
    palette = NULL, #NULL
    image_palette = NULL, #NULL
    manual_mapping = NULL, #NULL
    neighbor_weight = 0.5, # TODO: confirm default value
    radius = 90, # TODO: confirm default value
    n_neighbors = 16,
    ...
) {
  if (!is.null(manual_mapping)) {
    if (is.null(palette) && is.null(image_palette)) {
      print("Palette not provided, ignoring manual mapping to avoid color duplication with auto-generated ones.")
      manual_mapping <- list()
    } else {
      print(paste("Using manual color mapping for", names(manual_mapping), "."))
      print("Note that manual colors should be different with any color in `palette` or `image_palette`.")
      # Exclude cells with manually given color
      cell_coordinates <- cell_coordinates[!cell_labels %in% names(manual_mapping)]
      cell_labels <- cell_labels[!cell_labels %in% names(manual_mapping)]
    }
  } else {
    manual_mapping <- list()
  }
  if (!is.null(palette)) {
    # print(length(unique(cell_labels)))
    # print(length(palette))
    # stopifnot(length(unique(cell_labels)) <= length(palette),
    #           paste("Palette not sufficient for", length(unique(cell_labels)), "cell types."))
  }
  # Construct cluster spatial distance matrix based on cell neighborhood
  print("Calculating cluster distance graph...")
  cluster_distance_matrix <- spatial_distance(
    cell_coordinates = cell_coordinates,
    cell_labels = cell_labels,
    neighbor_weight = neighbor_weight,
    radius = radius,
    n_neighbors = n_neighbors
    # neighbor_args = list(...)
  )
  # Calculate color mapping
  color_mapping <- assign_color(
    cluster_distance_matrix = cluster_distance_matrix,
    palette = palette,
    image_palette = image_palette,
    colorblind_type = colorblind_type
    # mapping_args = list(...),
    # embed_args = list(...)
  )
  # Restore manual colors, reorder color mapping by cluster names.
  # color_mapping <- c(mapping_args, embed_args)
  # color_mapping <- color_mapping[order(names(color_mapping))]
  return(color_mapping)
}



