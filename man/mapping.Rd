% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping.R
\name{map_graph}
\alias{map_graph}
\alias{embed_graph}
\alias{cluster_mapping_iou}
\title{map the vertices between two graph}
\usage{
map_graph(
  cluster_distance,
  color_distance,
  random_seed = 123,
  distance_metric = "mul_1",
  random_max_iter = 5000,
  verbose = FALSE
)

embed_graph(
  cluster_distance,
  transformation = "umap",
  l_range = c(30, 80),
  log_colors = FALSE,
  trim_fraction = 0.0125
)

cluster_mapping_iou(cluster_label_mapping, cluster_label_reference)
}
\arguments{
\item{cluster_distance}{a data.frame representing the dissimilarity between clusters}

\item{color_distance}{a data.frame representing the perceptual difference between colors}

\item{random_seed}{Integer for random seed in heuristic solver}

\item{distance_metric}{Metric used for matrix mapping, default is manhattan}

\item{random_max_iter}{optional parameter}

\item{verbose}{Boolean flag for outputting info, default is FALSE}

\item{cluster_label_mapping}{List of cluster results for cells to be mapped}

\item{cluster_label_reference}{List of cluster results for cells to be mapped to}
}
\value{
optimized color mapping for clusters including cluster names and hex colors

A list where keys are cluster names and values are hex colors, representing the optimized color mapping

A list representing the mapping result of cluster_label_mapping
}
\description{
Function to embed the cluster distance graph into chosen colorspace, while keeping distance
relationship. Currently only supports CIE Lab space. Proper colors are selected within whole
colorspace based on the embedding of each cluster

Function to map clusters between different clustering results based on cluster overlap
}
