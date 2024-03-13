#' @rdname mapping
#' @title map the vertices between two graph
#' @name map_graph
#' @param cluster_distance a data.frame with unique cluster names as index and columns, which
#' contains a distance adjacent matrix for clusters, representing the dissimilarity between clusters
#' @param color_distance a data.frame with unique colors (in hex) as index and columns, which
#' contains a distance adjacent matrix for colors, representing the perceptual difference between colors.
#' @param random_seed random seed for heuristic solver
#' @param distance_metric metric used for matrix mapping. Defaults to "manhattan"
#' @param random_max_iter optional parameter
#' @param verbose output info
#' @return optimized color mapping for clusters including cluster names and hex colors
#' @export
#' @importFrom stats setNames
map_graph <- function(
    cluster_distance,
    color_distance,
    random_seed = 123,
    distance_metric = "mul_1",
    random_max_iter = 5000,
    verbose = FALSE
) {
  set.seed(random_seed)

  color_shuffle_index <- 1:nrow(color_distance)
  shuffle_distance <- 1e9

  for (iseed in 1:random_max_iter) {
    color_shuffle_index_tmp <- color_shuffle_index
    set.seed(iseed)
    # color_shuffle_index_tmp_test <- unlist(color_shuffle_index_tmp)
    color_shuffle_index_tmp <- sample(color_shuffle_index_tmp)
    color_distance_shuffle <- t(color_distance[color_shuffle_index_tmp, ])[color_shuffle_index_tmp, ]
    shuffle_distance_tmp <- matrix_distance(
      matrix_x = as.matrix(cluster_distance),
      matrix_y = color_distance_shuffle,
      metric = distance_metric
    )
    if (shuffle_distance_tmp < shuffle_distance) {
      shuffle_distance <- shuffle_distance_tmp
      color_shuffle_index <- color_shuffle_index_tmp
    }
  }
  for (siter in 1:(random_max_iter*2)) {
    color_shuffle_index_tmp <- color_shuffle_index
    idx <- sample(1:length(color_shuffle_index_tmp), 2)
    tmp <- color_shuffle_index_tmp[idx[2]]
    color_shuffle_index_tmp[idx[2]] <- color_shuffle_index_tmp[idx[1]]
    color_shuffle_index_tmp[idx[1]] <- tmp
    color_distance_shuffle <- t(color_distance[color_shuffle_index_tmp, ])[color_shuffle_index_tmp, ]
    shuffle_distance_tmp <- matrix_distance(
      matrix_x = as.matrix(cluster_distance),
      matrix_y = color_distance_shuffle,
      metric = distance_metric
    )
    if (shuffle_distance_tmp < shuffle_distance) {
      shuffle_distance <- shuffle_distance_tmp
      color_shuffle_index <- color_shuffle_index_tmp
    }
  }
  # Return color mapping dictionary, sorted by keys
  color_mapping <- setNames(rownames(color_distance)[color_shuffle_index], rownames(cluster_distance))
  # color_mapping <- setNames(color_distance[color_shuffle_index], rownames(cluster_distance))
  color_mapping <- color_mapping[order(names(color_mapping))]

  return(color_mapping)
}

#' @rdname mapping
#' @title embed the cluster distance graph into chosen colorspace
#' @name embed_graph
#' @description
#' Function to embed the cluster distance graph into chosen colorspace, while keeping distance
#' relationship. Currently only supports CIE Lab space. Proper colors are selected within whole
#' colorspace based on the embedding of each cluster
#' @param cluster_distance a data.frame representing the dissimilarity between clusters
#' @param color_distance a data.frame representing the perceptual difference between colors
#' @param random_seed Integer for random seed in heuristic solver
#' @param distance_metric Metric used for matrix mapping, default is manhattan
#' @param verbose Boolean flag for outputting info, default is FALSE
#' @return A list where keys are cluster names and values are hex colors, representing the optimized color mapping
#' @export
embed_graph <- function(
    cluster_distance,
    transformation = "umap",
    l_range = c(30, 80),
    log_colors = FALSE,
    trim_fraction = 0.0125
) {
  # Embed clusters into 3-dimensional space
  print("Calculating cluster embedding...")

  if (transformation == "mds") {
    mds_result <- isoMDS(dist(cluster_distance), k = 3)
    embedding <- mds_result$points
    set.seed(123)
  } else if (transformation == "umap") {
    embedding <- umap::umap(cluster_distance, n_components = 3)
    set.seed(123)
  }

  # Rescale embedding to CIE Lab colorspace
  print("Rescaling embedding to CIE Lab colorspace...")
  embedding <- embedding$layout

  embedding <- sweep(embedding, 2, quantile(embedding, probs = trim_fraction, na.rm = TRUE), FUN = "-")
  embedding[embedding < 0] <- 0
  embedding <- sweep(embedding, 2, quantile(embedding, probs = 1 - trim_fraction, na.rm = TRUE), FUN = "/")
  embedding[embedding > 1] <- 1

  if (log_colors) {
    embedding <- log10(embedding + max(quantile(embedding, probs = 0.05, na.rm = TRUE), 1e-3))
    embedding <- sweep(embedding, 2, apply(embedding, 2, min), FUN = "-")
    embedding <- sweep(embedding, 2, apply(embedding, 2, max), FUN = "/")
  }

  embedding[, 1] <- embedding[, 1] * (l_range[2] - l_range[1]) + l_range[1]
  embedding[, 2:3] <- (embedding[, 2:3] - 0.5) * 200
  print("Optimizing cluster color mapping...")
  lab_to_hex <- apply(embedding, 1, lab_to_hex)
  color_mapping <- setNames(as.list(lab_to_hex), rownames(cluster_distance))
  color_mapping <- unlist(color_mapping)
  return(color_mapping)
}

#' @rdname mapping
#' @title map clusters between different clustering results
#' @name cluster_mapping_iou
#' @description
#' Function to map clusters between different clustering results based on cluster overlap
#' @param cluster_label_mapping List of cluster results for cells to be mapped
#' @param cluster_label_reference List of cluster results for cells to be mapped to
#' @return A list representing the mapping result of cluster_label_mapping
#' @export

cluster_mapping_iou <- function(
    cluster_label_mapping,
    cluster_label_reference
    ) {
  iou <- function(i, j) {
    I <- sum(cluster_label_mapping == i & cluster_label_reference == j)
    U <- sum(cluster_label_mapping == i | cluster_label_reference == j)
    return(I / U)
  }
  # Cells should be identical between different runs.
  stopifnot(length(cluster_label_mapping) == length(cluster_label_reference))
  ufunc_iou <- Vectorize(iou, vectorize.args = c("i", "j"))
  cluster_label_mapping <- as.character(cluster_label_mapping)
  cluster_label_reference <- as.character(cluster_label_reference)

  # Reference label types should be more than mapping label types
  mapping_label_list <- unique(cluster_label_mapping)
  reference_label_list <- unique(cluster_label_reference)
  stopifnot(length(mapping_label_list) <= length(reference_label_list))

  # Grid label lists for vectorized calculation
  mapping_vector_column <- matrix(rep(mapping_label_list, each = length(reference_label_list)),
                                  nrow = length(mapping_label_list), ncol = length(reference_label_list), byrow = TRUE)
  reference_vector_index <- matrix(rep(reference_label_list, times = length(mapping_label_list)),
                                   nrow = length(mapping_label_list), ncol = length(reference_label_list))
  iou_matrix <- outer(mapping_vector_column, reference_vector_index, Vectorize(iou))

  # Greedy mapping to the largest IOU of each label
  relationship <- list()
  index_not_mapped <- rep(TRUE, length(mapping_label_list))
  iou_matrix_backup <- iou_matrix
  while(sum(iou_matrix) != 0) {
    which_max <- which(iou_matrix == max(iou_matrix), arr.ind = TRUE)
    reference_index <- which_max[1, 1]
    mapping_index <- which_max[1, 2]
    relationship[[mapping_label_list[mapping_index]]] <- reference_label_list[reference_index]
    # Clear mapped labels to avoid duplicated mapping
    iou_matrix[reference_index, ] <- 0
    iou_matrix[, mapping_index] <- 0
    index_not_mapped[mapping_index] <- FALSE
  }
  # Check if every label is mapped to a reference
  duplicate_map_label <- rep(1, length(reference_label_list))
  for (mapping_index in seq_along(index_not_mapped)) {
    is_force_map <- index_not_mapped[mapping_index]
    if (is_force_map) {
      reference_index <- which.max(iou_matrix_backup[, mapping_index])
      relationship[[mapping_label_list[mapping_index]]] <- paste0(
        reference_label_list[reference_index], ".", duplicate_map_label[reference_index]
      )
      duplicate_map_label[reference_index] <- duplicate_map_label[reference_index] + 1
    }
  }
  mapped_cluster_label <- sapply(cluster_label_mapping, function(x) relationship[[x]])

  mapped_cluster_label <- as.list(mapped_cluster_label)
  return(mapped_cluster_label)
}




