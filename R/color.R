#' @rdname color
#' @title Core color mapping function for Spaco2.
#' @name assign_color
#' @description
#' Assign Colors to Clusters Based on Distance Matrix.
#' This function assigns colors to clusters based on a given distance matrix.
#' It supports colorblind-friendly options and can automatically generate a color palette,
#' use a predefined palette, or extract colors from an image.
#' Spaco2 provides 3 basic color mapping mode in this function:
#' 1. Optimize the mapping of a pre-defined color palette.
#' 2. Extract colors from image.
#' 3. Automatically generate colors within colorspace.
#' @param cluster_distance_matrix A matrix representing the distances between clusters.A DataFrame with unique cluster names as `index` and `columns`, which contains a distance adjacent matrix for clusters, representing the dissimilarity between clusters.
#' @param colorblind_type A character vector specifying the type of colorblindness
#'   to accommodate. Options include "none", "protanopia", "deuteranopia", "tritanopia",
#'   and "general". Default is "none".
#' @param palette An optional vector of color values (in hex format). If provided,
#'   this palette will be used and `image_palette` will be ignored.Defaults to None.
#' @param image_palette An optional image (in a format compatible with R) used to
#'   extract a color palette. Ignored if `palette` is provided.Defaults to None.
#' @param mapping_args A list of additional arguments to pass to the `map_graph` function.
#' @param embed_args A list of additional arguments to pass to the `embed_graph` function.
#' @return A named vector where names are cluster identifiers and values are the assigned
#'   hex color codes.
#' @export
assign_color <- function(
    cluster_distance_matrix = data.frame(),
    colorblind_type = c("none", "protanopia", "deuteranopia", "tritanopia", "general"),
    palette = NULL, #NULL
    image_palette = NULL, #NULL
    ...
) {
  # Auto-generate a palette if not provided
  if (is.null(palette)) {
    print("palette not provided.")
    if (is.null(image_palette)) {
      # Mode 3
      print("Auto-generating colors from CIE Lab colorspace...")
      color_mapping <- embed_graph(
        cluster_distance = cluster_distance_matrix,
        # embed_args = list(...)
        )
      color_mapping <- color_mapping[order(names(color_mapping))]
      return(color_mapping)
    } else {
      # Mode 2
      print(paste("Using", image_palette, "..."))
      print("Drawing appropriate colors from provided image...")
      palette <- extract_palette(reference_image = image_palette, n_colors = length(cluster_distance_matrix), colorblind_type=colorblind_type)
    }
  }
  # Construct color perceptual distance matrix
  print("Calculating color distance graph...")
  color_distance_matrix <- perceptual_distance(colors = palette, colorblind_type = colorblind_type) + 1e-5
  # Map clusters and colors via graph
  print("Optimizing color mapping...")
  color_mapping <- map_graph(
    cluster_distance = cluster_distance_matrix,
    color_distance = color_distance_matrix
    # mapping_args = list(...)
  )
  color_mapping <- setNames(paste0("#", color_mapping),names(color_mapping))
  return(color_mapping)
}
