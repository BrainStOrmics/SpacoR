#' @rdname utils
#' @name matrix_distance
#' @title Calculate the distance between two matrices
#' @param matrix_x matrix x
#' @param matrix_y matrix y
#' @param metric metric used to calculate distance. Defaults to manhattan
#' @return the distance between two matrices
#' @export
matrix_distance <- function(matrix_x, matrix_y, metric = "manhattan") {
  stopifnot(dim(matrix_x) == dim(matrix_y)) # 确保矩阵尺寸相同

  metric <- tolower(metric)

  if (metric == "euclidean") {
    return(norm(matrix_x - matrix_y)^2)  #return(norm(matrix_x - matrix_y, type = "F"))
  } else if (metric == "manhattan") {
    return(sum(abs(matrix_x - matrix_y)))
  } else if (metric == "log") {
    return(sum(matrix_x * log(matrix_y)))
  } else if (metric == "mul_1") {
    return(1000 / sum(matrix_x * matrix_y))
  }
}

#' @rdname utils
#' @name hex_to_rgb
#' @title Convert hex string to RGB value
#' @param hex_codehex string representing a RGB color
#' @return integer values for RGB channels
#' @export
hex_to_rgb <- function(hex_code) {
  # print(paste("Received hex_code in hex_to_rgb:", hex_code))  # 打印接收到的颜色代码
  hex_code <- gsub("#", "", hex_code)
  if (length(hex_code) != 1 || nchar(hex_code) != 6) {
    stop("Invalid color code in hex_to_rgb. Expected a single 6-character hex string.")
  }
  # print(paste("hex_code before parsing:", hex_code))
  r <- strtoi(substr(hex_code, 1, 2), 16L)
  g <- strtoi(substr(hex_code, 3, 4), 16L)
  b <- strtoi(substr(hex_code, 5, 6), 16L)
  # print(paste("Parsed RGB values:", r, g, b))
  return(c(r, g, b))
}

#' @rdname utils
#' @name rgb_to_hex
#' @title Convert RGB value (0~255) to hex string
#' @param rgb_code RGB channel value
#' @return color hex string
#' @export
rgb_to_hex <- function(rgb_code) {
  if (length(rgb_code) != 3 || !all(sapply(rgb_code, is.numeric))) {
    stop("rgb_code must be a numeric vector with three elements.")
  }
  rgb_code_scaled <- rgb_code / 255
  hex_color <- rgb(rgb_code_scaled[1], rgb_code_scaled[2], rgb_code_scaled[3], maxColorValue = 1)
  return(hex_color)
}

#' @rdname utils
#' @name lab_to_hex
#' @title Convert CIE Lab color value to hex string
#' @param lab_code CIE Lab color value
#' @return color hex string
#' @export
lab_to_hex <- function(lab_code) {
  lab_object <- LAB(lab_code[1], lab_code[2], lab_code[3])
  xyz_object <- convert_color(lab_object, XYZColor)
  rgb_object <- convert_color(xyz_object, sRGBColor)

  # Clip to legal RGB color
  rgb_object@coords[, "R"] <- pmin(pmax(rgb_object@coords[, "R"], 0), 1) #不超过1且不小于0
  rgb_object@coords[, "G"] <- pmin(pmax(rgb_object@coords[, "G"], 0), 1)
  rgb_object@coords[, "B"] <- pmin(pmax(rgb_object@coords[, "B"], 0), 1)
  hex_color <- rgb(rgb_object, maxColorValue = 255)
}

#' @rdname utils
#' @name rgb_to_lms_img
#' @title Convert RGB image matrix (0~255) to lms image matrix
#' @param img RGB image
#' @return lms image matrix
#' @export
rgb_to_lms_img <- function(img) {
  # Implementation from https://gitlab.com/FloatFlow
  lms_matrix <- matrix(c(0.3904725, 0.07092586, 0.02314268,
                         0.54990437, 0.96310739, 0.12801221,
                         0.00890159, 0.00135809, 0.93605194), byrow = TRUE, nrow = 3)
  lms_img <- apply(img, 3, function(x) x %*% lms_matrix[, ])
  return(lms_img)
}

#' @rdname utils
#' @name lms_to_rgb_img
#' @title Convert lms image matrix to RGB image matrix
#' @param img lms image
#' @return RGB image matrix
#' @export
lms_to_rgb_img <- function(img) {
  ## Implementation from https://gitlab.com/FloatFlow
  rgb_matrix <- matrix(c(2.85831110, -0.210434776, -0.0418895045,
                         -1.62870796, 1.15841493, -0.118154333,
                         -0.0248186967, 0.000320463334, 1.06888657), byrow = TRUE, nrow = 3)
  # Apply matrix multiplication across the 3rd dimension of the image array
  rgb_img =
  apply(img, c(1, 2), function(x) x %*% rgb_matrix[, ])
  rgb_img <- pmax(rgb_img, 0)
  rgb_img <- pmin(rgb_img, 255)
  return(rgb_img)
}

#' @rdname utils
#' @name color_difference_rgb
#' @title Calculate the perceptual difference between colors
#' @param color_x color_x
#' @param color_y color_x
#' @return the perceptual difference between two colors
#' @export
color_difference_rgb <- function(color_x, color_y) {
  color_x <- gsub("#", "", color_x)
  color_y <- gsub("#", "", color_y)
  # print(paste("color_x:", color_x, "color_y:", color_y))
  if (!is.character(color_x) || length(color_x) != 1 || nchar(color_x) != 6 ||
      !is.character(color_y) || length(color_y) != 1 || nchar(color_y) != 6) {
    stop("Invalid color code. Color codes must be single 6-character hex strings.")
  }
  color_x <- gsub("#", "", color_x)
  rgb_x <- hex_to_rgb(color_x)
  color_y <- gsub("#", "", color_y)
  rgb_y <- hex_to_rgb(color_y)
  rgb_r <- (rgb_x[1] + rgb_y[1]) / 2

  return(sqrt(
    (2 + rgb_r / 256) * (rgb_x[1] - rgb_y[1])^2 +
      4 * (rgb_x[2] - rgb_y[2])^2 +
      (2 + (255 - rgb_r) / 256) * (rgb_x[3] - rgb_y[3])^2
  ))
}

#' @rdname utils
#' @name get_bin_color
#' @title Revert bin number in `extract_palette` function to Lab values
#' @param bin_number numbered bin color
#' @return Lab values for the centroid color of this bin
#' @export
get_bin_color <- function(bin_number) {
  l <- bin_number %/% 400 * 5 + 2.5
  bin_number <- bin_number %% 400

  a <- bin_number %/% 20 * 12.75 + 6.375 - 127
  bin_number <- bin_number %% 20

  b <- bin_number * 12.75 + 6.375 - 127

  return(c(l, a, b))
}


#' @rdname utils
#' @name palette_min_distance
#' @title Calculate the minimal distance within a Lab palette
#' @param palette a Lab palette
#' @return the min distance
#' @export
palette_min_distance <- function(palette) {
  n_palette <- nrow(palette)
  distance_matrix <- matrix(1e10, nrow = n_palette, ncol = n_palette)
  for (i in 1:n_palette) {
    for (j in i:n_palette) {
      if (i == j) {
        next
      } else {
      lab_distance <- ciede2000(lab1 = palette[i,], lab2 = palette[j,])$delta_E_00
      distance_matrix[i, j] <- lab_distance
      distance_matrix[j, i] <- lab_distance
      }
    }
  }
  return(min(c(distance_matrix)))
}

#' @rdname utils
#' @name color_score
#' @title Score color replacement
#' @param lab_color lab_color
#' @param color_count color_count
#' @param palette palette
#' @param wn float
#' @return color_score
#' @export
color_score <- function(lab_color, color_count, palette, wn) {
  l <- lab_color[1]
  a <- lab_color[2]
  b <- lab_color[3]

  # na
  n <- color_count

  # da
  e_distance <- rep(1e10, nrow(palette))
  for (i in length(palette)) {
    point = palette[i]
    e_distance[i] <- ciede2000(lab1=c(l, a, b), lab2=point)$delta_E_00
  }
  dist_e2000 <- min(e_distance)

  # Sa
  # Calculate Luv color, drop L
  Lab_2_u_UF <- Vectorize(function(l, a, b) {
    Luv <- XYZ_to_Luv(Lab_to_XYZ(l, a, b))
    return(Luv[1, "u"])
  })

  Lab_2_v_UF <- Vectorize(function(l, a, b) {
    Luv <- XYZ_to_Luv(Lab_to_XYZ(l, a, b))
    return(Luv[1, "v"])
  })

  color_u <- XYZ_to_Luv(Lab_to_XYZ(c(l, a, b)))[1]
  color_v <- XYZ_to_Luv(Lab_to_XYZ(c(l, a, b)))[2]
  color_uv <- c(color_u, color_v)

  palette_u <- apply(palette, 1, function(p) XYZ_to_Luv(Lab_to_XYZ(p))[1])
  palette_v <- apply(palette, 1, function(p) XYZ_to_Luv(Lab_to_XYZ(p))[2])
  palette_uv <- cbind(palette_u, palette_v)

  convex_hull_palette <- convhulln(palette_uv)
  distance_to_convex_hull <- max(rowSums(convex_hull_palette$equations[,1:2] * color_uv) + convex_hull_palette$equations[,3])

  # S(ca)
  ws <- 0.15
  return(ws * distance_to_convex_hull + wn * color_count + dist_e2000)

}

#' @rdname utils
#' @name simulate_cvd
#' @param palette_hex palette_hex
#' @param colorblind_type optional parameter
#' @return rgb_to_hex
#' @export
simulate_cvd <- function(palette_hex, colorblind_type) {
  ## Modified from https://gitlab.com/FloatFlow
  palette_rgb_img <- t(sapply(palette_hex, hex_to_rgb))
  palette_lms_img <- rgb_to_lms_img(palette_rgb_img)

  if (tolower(colorblind_type) == 'protanopia') {
    sim_matrix <- matrix(c(0, 0, 0, 0.90822864, 1, 0, 0.008192, 0, 1), nrow = 3, byrow = TRUE)
  } else if (tolower(colorblind_type) == 'deuteranopia') {
    sim_matrix <- matrix(c(1, 0, 0, 0, 0, 1.10104433, 0, -0.00901975, 0, 0, 1), nrow = 3, byrow = TRUE)
  } else if (tolower(colorblind_type) == 'tritanopia') {
    sim_matrix <- matrix(c(1, 0, 0, 0, 1, 0, -0.15773032, 1.19465634, 0), nrow = 3, byrow = TRUE)
  } else {
    stop(paste(colorblind_type, "is an unrecognized colorblindness type."))
  }

  palette_lms_img <- t(palette_lms_img)
  palette_lms_img <- apply(palette_lms_img, c(1, 2), function(x) sum(x * sim_matrix))
  palette_rgb_img <- lms_to_rgb_img(simulated_lms_img)

  return(sapply(palette_rgb_img[[1]], rgb_to_hex))
}

#' @rdname utils
#' @name extract_palette
#' @title Extract palette from image
#' @param reference_image reference_image
#' @param n_colors colors
#' @param colorblind_type colorblind type
#' @param l_range intergrate
#' @param trim_percentile trim_percentile
#' @param max_iteration max_iteration
#' @param verbose verbose
#' @return lab_to_hex
#' @export
extract_palette <- function(reference_image, n_colors, colorblind_type,
                            l_range = c(20, 85), trim_percentile = 0.03,
                            max_iteration = 20, verbose = FALSE) {

  lab_image <- convertColor(reference_image, from = "sRGB", to = "LAB")
  # Make L, A, B into 20 bins each.
  print("Extracting color bins...")
  bin_index_l <- as.integer(lab_image[,,1] / 5)
  bin_index_a <- as.integer((lab_image[,,2] + 127) / 2.55 / 5)
  bin_index_b <- as.integer((lab_image[,,3] + 127) / 2.55 / 5)
  bin_index_l[bin_index_l == 20] <- 19
  bin_index_a[bin_index_a == 20] <- 19
  bin_index_b[bin_index_b == 20] <- 19

  # Get bin frequency
  numbered_bin_colors <- bin_index_l * 400 + bin_index_a * 20 + bin_index_b
  bin_color_table <- table(numbered_bin_colors)
  bin_color_set <- as.integer(names(bin_color_table))
  bin_color_count <- as.integer(bin_color_table)

  # Truncate color set by l_range, filter out infrequent colors
  filter <- (bin_color_set > (l_range[1] / 5 * 400)) &
    (bin_color_set < (l_range[2] / 5 * 400)) &
    (bin_color_count > quantile(bin_color_count, trim_percentile))
  bin_color_set <- bin_color_set[filter]
  bin_color_count <- bin_color_count[filter]
  lab_color_set <- vapply(bin_color_set, get_bin_color, FUN.VALUE = numeric(3))
  # Initiate palette, chosen by frequency and distinction
  print("Initiating palette...")
  sigma <- 80
  palette <- list()

  for (i in 1:n_colors) {
    selected_index <- which.max(bin_color_count)
    palette[[i]] <- lab_color_set[selected_index, ]

    # Update frequency to punish chosen colors and similar ones
    for (j in seq_along(bin_color_count)) {
      if (j == selected_index) {
        bin_color_count[j] <- -1
      } else {
        lab_distance <- ciede2000(lab1 = lab_color_set[selected_index, ],
                                  lab2 = lab_color_set[j, ])$delta_E_00
        if (colorblind_type != "none") {
          if (colorblind_type == "general") {
            lab_distance <- lab_distance / 4
            for (cb_t in c("protanopia", "deuteranopia", "tritanopia")) {
              color_cvd <- simulate_cvd(c(lab_to_hex(lab_color_set[selected_index, ]),
                                          lab_to_hex(lab_color_set[j, ])),
                                        colorblind_type = cb_t)
              lab_distance <- lab_distance +
                color_difference_rgb(color_cvd[1], color_cvd[2]) / 4
            }
          } else {
            color_cvd <- simulate_cvd(c(lab_to_hex(lab_color_set[selected_index, ]),
                                        lab_to_hex(lab_color_set[j, ])),
                                      colorblind_type = colorblind_type)
            lab_distance <- color_difference_rgb(color_cvd[1], color_cvd[2])
          }
        }
        bin_color_count[j] <- bin_color_count[j] * (1 - exp(-(lab_distance / sigma)^2))
      }
    }
  }

  # Convert the list to an array
  palette <- do.call(rbind, palette)

  # Optimize palette distance
  print("Optimizing extracted palette...")
  palette_distance <- numeric(max_iteration) #创建数值向量，默认值0

  for (i in 1:max_iteration) {
    # Update palette distance
    palette_distance[i] <- palette_min_distance(palette)
    if (verbose) {
      # print(paste("Palette distance:", palette_distance[i]))
    } else if(i >= 4 && sd(palette_distance[i - 4 : i]) < 1e-6) {
      break
    }

    # Replace color starting from the one with the least frequency
    for (j in n_colors:1) {
      sample_index <- farthest_point_sample(lab_color_set, palette[-j, ], n_sample=3)
      sample_score <- rep(1e10, 3)

      for (k in seq_along(sample_index)) {
        idx <- sample_index[k]
        sample_score[k] <- color_score(
          lab_color = lab_color_set[idx, ],
          color_count = bin_color_count[idx],
          palette = palette[-j, ],
          wn = prod(dim(reference_image)[1:2]) * 0.0003
        )
      }

      # Update the palette with the color with the maximum score
      palette[j, ] <- lab_color_set[sample_index[which.max(sample_score)], ]
    }
  }
  palette_hex <- apply(palette, 2, lab_to_hex)
  return(palette_hex)
}
