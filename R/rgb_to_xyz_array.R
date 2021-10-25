#' customize the rgb_to_xyz function (from schemr) to work on arrays ##
#' @export
rgb_to_xyz_array <- function (rgb, transformation = "sRGB", linear_func = NULL)
{
  if (class(transformation) == "Matrix") {
    if (all(dim(transformation) == c(3L, 3L))) {
      m <- transformation
    }
    else {
      stop("Transformation should be a 3x3 matrix.")
    }
  }
  else {
    transformation_match <- match.arg(transformation, names(schemr:::transformation_matrix),
                                      several.ok = FALSE)
    m <- schemr:::transformation_matrix[[transformation_match]]
  }
  temp_r <- unlist(rgb[ 1], FALSE, FALSE)/255
  temp_g <- unlist(rgb[ 2], FALSE, FALSE)/255
  temp_b <- unlist(rgb[ 3], FALSE, FALSE)/255
  if (transformation_match == "sRGB") {
    temp_r <- schemr:::srgb_transformation(temp_r)
    temp_g <- schemr:::srgb_transformation(temp_g)
    temp_b <- schemr:::srgb_transformation(temp_b)
  }
  else if (transformation_match == "Adobe") {
    temp_r <- adobe_transformation(temp_r)
    temp_g <- adobe_transformation(temp_g)
    temp_b <- adobe_transformation(temp_b)
  }
  else if (!is.null(linear_func)) {
    temp_r <- linear_func(temp_r)
    temp_g <- linear_func(temp_g)
    temp_b <- linear_func(temp_b)
  }
  else {
    temp_r <- temp_r * 100
    temp_g <- temp_g * 100
    temp_b <- temp_b * 100
  }
  xyz <- matrix(nrow = 1, ncol = 3)

  xyz[,1] <- temp_r * m[1, 1] + temp_g * m[1, 2] + temp_b * m[1, 3]
  xyz[,2] <- temp_r * m[2, 1] + temp_g * m[2, 2] + temp_b * m[2, 3]
  xyz[,3] <- temp_r * m[3, 1] + temp_g * m[3, 2] + temp_b * m[3, 3]
  return(t(xyz))
}
