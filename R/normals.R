## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Function to generate normals for a point cloud
#' @description Returns a function which estimates normals for a
#'   3D point cloud.
#' @param x \code{integer}. Number of neighbors used to estimate the normals.
#' @param method One of \code{"PCA"} to estimate the normal direction at
#'   each point by linear least squares fitting of a plane over its nearest
#'   neighbors, or \code{"Jet"} to estimate the normal direction at each
#'   point by fitting a jet surface over its nearest neighbors. See details.
#' @returns A function which takes one argument: a numeric matrix with
#'   three columns, each row represents a point, and the function returns a
#'   matrix of the same size as the input matrix, with each row giving one
#'   unit normal per point.
#' @details See \url{https://doc.cgal.org/latest/Point_set_processing_3/} for details.
#'   The \code{getNormalsFun} function is intended to be used in the
#'   \code{\link[SurfaceMesh]{reconstructPoisson}} function. If you want to use it for
#'   another purpose, be careful because the function it returns does not
#'   check the matrix it takes as argument.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' fun      <- getNormalsFun(6)
#' mesh_psr <- reconstructPoisson(mesh[["vertices"]],
#'                                normalsFun=fun,
#'                                smAngle=10,
#'                                smRadius=1.5,
#'                                smDistance=0.3)
#'
#' mesh_psr_rgl <- toRGL(mesh_psr)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_psr_rgl)
#'
#' @export
getNormalsFun <- function(x, method = c("PCA", "Jet")) {
  method_choices <- c("pca", "jet")
  method         <- match.arg(tolower(method), choices=method_choices)
  x <- as.integer(x)
  if(x < 2L) {
    stop("There must be at least two neighbors.", call. = TRUE)
  }
  fun <- function(points) {
    if(!is.matrix(points) || !is.numeric(points) || (ncol(points) != 3L) || (nrow(points) <= 3L)) {
      stop("The `points` argument must be a numeric matrix with 3 columns and at least 4 points.", call. = TRUE)
    }
    storage.mode(points) <- "double"
    if(anyNA(points)) {
      stop("`points` may not have missing values.", call. =TRUE)
    }
    normals_jet_pca_cpp(t(points), x, method)
  }
  class(fun) <- "CGALnormalsFun"
  fun
}
