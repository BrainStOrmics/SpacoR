#' @rdname distance
#' @title calculate spatial interlacement distance graph for cell clusters
#' @name spatial_distance
#' @description
#' Function to calculate spatial interlacement distance graph for cell clusters, where
#' we define the interlacement distance as the number of neighboring cells between two cluster
#' @param cell_coordinates a list like object containing spatial coordinates for each cell
#' @param cell_labels a list like object containing cluster labels for each cell
#' @param neighbor_weight cell weight to calculate cell neighborhood. Defaults to 0.5
#' @param radius radius used to calculate cell neighborhood. Defaults to 90
#' @param n_neighbors k for KNN neighbor detection. Defaults to 16
#' @param n_cells nly calculate neighborhood with more than `n_cells`. Defaults to 3
#' @return a Data.Frame with unique cluster names as `index` and `columns`,
#' which contains interlacement distance between clusters
#' @export
#' @importFrom FNN get.knnx

spatial_distance <- function( #TODO:optimize neighbor calculation
    cell_coordinates,
    cell_labels,
    neighbor_weight = 0.5,
    radius = 90,
    n_neighbors = 16,
    n_cells = 3
) {
  unique_clusters <- unique(cell_labels)
  cluster_index <- as.list(seq_along(unique_clusters))
  names(cluster_index) <- unique_clusters

  # Calculate neighborhoods for all cells
  print("Calculating cell neighborhood...")
  cell_coordinates <- t(cell_coordinates)
  knn_result <- get.knnx(data = cell_coordinates, query = cell_coordinates, k = n_neighbors)
  neighbor_index_knn <- knn_result$nn.index
  neighbor_distance_knn <- knn_result$nn.dist

  # Intersection between knn neighbors and radius neighbors
  print("Filtering out neighborhood outliers...")
  neighbor_index_filtered <- list()
  neighbor_distance_filtered <- list()

  for (i in 1:dim(cell_coordinates)[1]) {
    # filter by radius is equalized to intersection
    neighbor_index_filtered[[i]] <- neighbor_index_knn[i,][neighbor_distance_knn[i,] <= radius]
    neighbor_distance_filtered[[i]] <- neighbor_distance_knn[i,][neighbor_distance_knn[i,] <= radius]
    # filter banished cell
    if (sum(cell_labels[neighbor_index_filtered[[i]]] == cell_labels[i]) < n_cells) {
      # keep an empty network with only cell i itself
      neighbor_index_filtered[[i]] <- i
      neighbor_distance_filtered[[i]] <- 0
    }
    # neighbor_index_filtered <- c(neighbor_index_filtered, neighbor_index_filtered[[i]])
    # neighbor_distance_filtered <- c(neighbor_distance_filtered, neighbor_distance_filtered[[i]])
  }
  neighbor_index_filtered <- as.list(neighbor_index_filtered)
  neighbor_distance_filtered <- as.list(neighbor_distance_filtered)

  # Calculate score matrix
  print("Calculating cluster interlacement score...")
  score_matrix <- matrix(0, nrow = length(unique_clusters), ncol = length(unique_clusters))
  rownames(score_matrix) <- unique_clusters
  colnames(score_matrix) <- unique_clusters

  # print("Dimnames of score matrix:")
  # print(dimnames(score_matrix))

  for (cell_i in seq_along(neighbor_index_filtered)) {
    size_n_i <- length(neighbor_index_filtered[[cell_i]])
    if (size_n_i == 0) {
      next
    }

    cell_cluster_i <- cell_labels[cell_i]
    # cell_cluster_i <- cluster_index[cell_labels[cell_i]]
    for (j in seq(2, size_n_i)) {
      # cell_cluster_j <- cluster_index[cell_labels[neighbor_index_filtered[[cell_i]][j]]]
      cell_cluster_j <- cell_labels[neighbor_index_filtered[[cell_i]][j]]
      # 尝试索引score_matrix之前，打印cell_cluster_i和cell_cluster_j的值
      # print(paste("cell_cluster_i:", cell_cluster_i, "cell_cluster_j:", cell_cluster_j))
      inversed_euclidean <- 1 / neighbor_distance_filtered[[cell_i]][j]
      if (!(cell_cluster_i %in% rownames(score_matrix) && cell_cluster_j %in% colnames(score_matrix))) {
        # print(paste("Invalid indices:", cell_cluster_i, cell_cluster_j))
      } else {
        score_matrix[cell_cluster_i, cell_cluster_j] <- score_matrix[cell_cluster_i, cell_cluster_j] + (inversed_euclidean / size_n_i)
      # score_matrix[cell_cluster_i, cell_cluster_j] <- score_matrix[cell_cluster_i, cell_cluster_j] + (inversed_euclidean / size_n_i)
    }
  }}
  # Keep maximum between score_matrix[x][y] and score_matrix[y][x], set diagonal to zero
  for (x in seq_along(unique_clusters)) {
    for (y in seq_along(unique_clusters)) {
      if (x == y) {
        score_matrix[x, y] <- 0
      } else {
        score_matrix[x, y] <- max(score_matrix[x, y], score_matrix[y, x])
      }
    }
  }

  cluster_interlace_matrix <- score_matrix

  # print("Constructing cluster interlacement graph...")
  cluster_interlace_matrix <- as.data.frame(cluster_interlace_matrix)
  rownames(cluster_interlace_matrix) <- unique_clusters
  colnames(cluster_interlace_matrix) <- unique_clusters

  return(cluster_interlace_matrix)
}


#' @rdname distance
#' @title calculate color perceptual difference matrix
#' @name perceptual_distance
#' @description
#' See 'color_difference_rgb' for details
#' @param colors a list of colors (in hex)
#' @param colorblind_type optional parameter
#' @return a data.frame with unique colors (in hex) as index and columns,which contains perceptual distance between colors
#' @export
perceptual_distance <- function(
  colors,
  colorblind_type = c("none","protanopia", "deuteranopia", "tritanopia", "general")
){
  colors <- gsub("#", "", colors)
  # 验证颜色代码格式是否符合要求
  for (color in colors) {
    if (nchar(color) != 6) {
      stop(paste("Invalid color code:", color, "Color code must be a 6-character hex string."))
    }
  }

  difference_matrix <- matrix(0, nrow = length(colors), ncol = length(colors))

  if (colorblind_type != "none") {
    print(paste("Calculating color perceptual distance under", colorblind_type, "..."))

    if (colorblind_type == "general") {
      colors <- gsub("#", "", colors)
      for (i in seq_along(colors)) {
        for (j in seq_along(colors)) {
          difference_matrix[i, j] <- difference_matrix[i, j] + color_difference_rgb(colors[i], colors[j]) / 4
        }
      }

      for (cb_t in c("protanopia", "deuteranopia", "tritanopia")) {
        color_cvd <- simulate_cvd(colors, colorblind_type = cb_t)
        # Calculate difference between cvd colors
        for (i in seq_along(color_cvd)) {
          for (j in seq_along(color_cvd)) {
            difference_matrix[i, j] <- difference_matrix[i, j] + color_difference_rgb(color_cvd[i], color_cvd[j]) / 4
          }
        }
      }
    } else {
      color_cvd <- simulate_cvd(colors, colorblind_type = colorblind_type)
    }
  } else {
    # Calculate difference between colors
    # print("Calculating color perceptual distance...")
    colors <- gsub("#", "", colors)
    for (i in seq_along(colors)) {
      # print(colors)
      for (j in seq_along(colors)) {
        difference_matrix[i, j] <- color_difference_rgb(colors[i], colors[j])
        if (colorblind_type == "none") {
          # print("Calculating color perceptual distance...")
          for (i in seq_along(colors)) {
            for (j in seq_along(colors)) {
              # 在调用 color_difference_rgb 之前打印颜色代码
              # print(paste("Comparing colors:", colors[i], "and", colors[j]))
              difference_matrix[i, j] <- color_difference_rgb(colors[i], colors[j])
            }
          }
        }
      }
    }
  }

  min_difference <- 1e9
  for (i in seq_along(colors)) {
    for (j in seq_along(colors)) {
      if (i == j) {
        next
      }
      min_difference <- min(min_difference, difference_matrix[i, j])
    }
  }
  print("Constructing color distance graph...")
  difference_matrix = as.data.frame(difference_matrix)
  rownames(difference_matrix) <- colors
  colnames(difference_matrix) <- colors
  print(paste("Difference of the most similar pair in the palette is", min_difference))

  return(difference_matrix)
}
