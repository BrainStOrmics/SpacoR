% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colorization.R
\name{colorize}
\alias{colorize}
\title{Colorize cell clusters based on spatial distribution}
\usage{
colorize(
  cell_coordinates,
  cell_labels,
  colorblind_type = c("none", "protanopia", "deuteranopia", "tritanopia", "general"),
  palette = NULL,
  image_palette = NULL,
  manual_mapping = NULL,
  neighbor_weight = 0.5,
  radius = 90,
  n_neighbors = 16
)
}
\arguments{
\item{cell_coordinates}{a list like object containing spatial coordinates for each cell.}

\item{cell_labels}{a list like object containing cluster labels for each cell.}

\item{colorblind_type}{Optional parameter.}

\item{palette}{a list of colors (in hex). If given, \code{image_palette}will be ignored. See \verb{Mode 1} above. Defaults to None.}

\item{image_palette}{an image in numpy array format. Should be a typical RGB image of shape (x, y, 3).
Ignored if \code{palette} is given. See \verb{Mode 2}above. Defaults to None.}

\item{manual_mapping}{a data structure for manual color mapping including cluster names and manually assigned colors (in hex).}

\item{neighbor_weight}{Weight for calculating cell neighborhood.Defaults to 0.5.}

\item{radius}{radius used to calculate cell neighborhood.Defaults to 90.}

\item{n_neighbors}{k for KNN neighbor detection.Defaults to 16.}

\item{neighbor_args}{arguments passed to \code{spatial_distance} function.}

\item{mapping_args}{arguments passed to \code{map_graph} function.}

\item{embed_args}{arguments passed to \code{embed_graph} function.}
}
\value{
Optimized color mapping for clusters, including cluster names and corresponding hex
}
\description{
Colorize cell clusters based on spatial distribution, so that spatially interlaced and spatially neighboring clusters are assigned with more perceptually different colors.
SpacoR provides 3 basic color mapping mode:
\enumerate{
\item Optimize the mapping of a pre-defined color palette.
\item Extract colors from image.
\item Automatically generate colors within colorspace.
}
}
