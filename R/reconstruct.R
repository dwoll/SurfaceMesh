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

#' @title Compute average spacing input points
#' @description Computes the average spacing of all input points to their
#'   \code{nNeighbors} nearest neighbor points. This provides an order of a
#'   point set density, useful in a surface reconstruction pipeline.
#' @param x Numeric matrix with 3 columns which stores the points, one point per row.
#' @param nNeighbors \code{integer}. How many neighboring points to consider
#'   for distance calculation.
#' @returns The average spacing.
#' @seealso \code{\link[SurfaceMesh]{removeOutliers}},
#'    \code{\link[SurfaceMesh]{reconstructPoisson}},
#'    \code{\link[SurfaceMesh]{reconstructSSS}},
#'    \code{\link[SurfaceMesh]{alphaWrap}}
#' @details See \url{https://doc.cgal.org/latest/Point_set_processing_3/}
#'   for details.
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' pts <- matrix(rnorm(120L), ncol=3)
#' getAverageSpacing(pts, nNeighbors=6L)
#'
#' @export
getAverageSpacing <- function(x, nNeighbors = 6L) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L)) {
    stop("The `x` argument must be a numeric matrix with 3 columns.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Missing values in the `x` matrix are not allowed.", call. = TRUE)
  }
  if(nrow(x) <= ncol(x)) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  stopifnot(isStrictPositiveInteger(nNeighbors))
  getAverageSpacing_cpp(t(x), as.integer(nNeighbors))
}

#' @title Remove outliers from point cloud
#' @description Remove outliers from point cloud based on their distance to
#'   neighboring points.
#' @param x Numeric matrix with 3 columns which stores the points, one point per row.
#' @param nNeighbors \code{integer}. How many neighboring points to consider
#'   for distance calculation.
#' @param threshPerc \code{numeric}. The maximum percentage of points to remove.
#'   If specified, \code{threshDst} should not be specified, otherwise
#'   \code{threshDst} takes precedence.
#' @param threshDst \code{numeric}. The minimum distance for a point to be
#'   considered as outlier. If missing and \code{threshPerc} is missing as well,
#'   twice the average spacing is chosen as a reasonable default value.
#' @returns A numeric matrix of points with outliers removed.
#' @seealso \code{\link[SurfaceMesh]{reconstructAFS}},
#'    \code{\link[SurfaceMesh]{reconstructPoisson}},
#'    \code{\link[SurfaceMesh]{reconstructSSS}},
#'    \code{\link[SurfaceMesh]{alphaWrap}}
#' @details See \url{https://doc.cgal.org/latest/Point_set_processing_3/}
#'   for details.
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' ptsBefore   <- matrix(rnorm(120L, mean=0, sd=1), ncol=3L)
#' (avgSpacing <- getAverageSpacing(ptsBefore, nNeighbors=6L))
#' ## contaminate
#' idx              <- sample(seq_len(nrow(ptsBefore)), 20L, replace=FALSE)
#' ptsBefore[idx, ] <- matrix(rnorm(5L, mean=10, sd=1), ncol=3L)
#' ptsAfter         <- removeOutliers(ptsBefore, nNeighbors=6L)
#' dim(ptsBefore)
#' dim(ptsAfter)
#'
#' @export
removeOutliers <- function(x,
                           nNeighbors = 6L,
                           threshPerc,
                           threshDst) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L)) {
    stop("The `x` argument must be a numeric matrix with 3 columns.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Missing values in the `x` matrix are not allowed.", call. = TRUE)
  }
  if(nrow(x) <= ncol(x)) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  stopifnot(isStrictPositiveInteger(nNeighbors))
  if(missing(threshPerc)) {
    threshPerc <- -1.0
  } else {
    stopifnot(isPositiveNumber(threshPerc), threshPerc < 100)
  }
  if(missing(threshDst)) {
    threshDst <- -1.0
  } else {
    stopifnot(isPositiveNumber(threshDst))
  }
  pointsOut <- removeOutliers_cpp(
    t(x), as.integer(nNeighbors), threshPerc, threshDst)
  t(pointsOut)
}

#' @title Advancing front surface reconstruction
#' @description Reconstruction of a surface mesh from a cloud of 3D points.
#'
#' @param x Numeric matrix with 3 columns which stores the points, one point per row.
#' @param jetSmoothing Optional integer >= 2. If specified,
#'   the point cloud is smoothed before the reconstruction, using
#'   this integer as the number of neighbors for the smoothing. Note that this
#'   smoothing preprocessing relocates the points and then should not be used
#'   if the points have been sampled without noise on the surface.
#' @param repairSoup Boolean. Clean the mesh (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices)?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/}
#'   for details.
#' @seealso \code{\link[SurfaceMesh]{removeOutliers}},
#'    \code{\link[SurfaceMesh]{reconstructSSS}},
#'    \code{\link[SurfaceMesh]{reconstructPoisson}}, \code{\link[SurfaceMesh]{alphaWrap}}
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh         <- dataHeart1
#' mesh_rgl     <- toRGL(mesh)
#' mesh_afs     <- reconstructAFS(mesh[["vertices"]], jetSmoothing=30)
#' mesh_afs_rgl <- toRGL(mesh_afs)
#'
#' mfrow3d(1, 2)
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_rgl, col="gray")
#' next3d()
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_afs_rgl, col="gray")
#'
#' @export
reconstructAFS <- function(x, jetSmoothing, repairSoup=TRUE, normals=FALSE) {
  if(!is.matrix(x) || !is.numeric(x)) {
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of points in `x`.", call. = TRUE)
  }
  if(anyNA(x)) {
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(!missing(jetSmoothing)) {
    stopifnot(isPositiveInteger(jetSmoothing), jetSmoothing >= 2L)
  } else {
    jetSmoothing <- 0L
  }
  stopifnot(isBoolean(repairSoup))
  mesh_cpp <- reconstructAFS_cpp(t(x), as.integer(jetSmoothing), repairSoup, normals)
  fromCPP(mesh_cpp)
}

#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'   Performs well if the aim is to approximate a noisy point cloud with a smooth surface.
#'   Not appropriate if the surface is expected to interpolate the input points.
#' @param x Numeric matrix with 3 columns which stores the points, one point per row.
#' @param normalsFun A function to generate normals, e.g., as returned by
#'   \code{\link[SurfaceMesh]{getNormalsFun}}.
#' @param spacing Optional size parameter. Smaller values increase the precision
#'   of the output mesh at the cost of higher computation time. If missing,
#'   an average spacing is chosen as a reasonable default value.
#' @param smAngle Bound for the minimum facet angle in degrees.
#' @param smRadius Relative bound for the radius of the surface Delaunay balls.
#' @param smDistance Relative bound for the center-center distances.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/}
#'   for details.
#' @seealso \code{\link[SurfaceMesh]{removeOutliers}},
#'    \code{\link[SurfaceMesh]{reconstructAFS}},
#'    \code{\link[SurfaceMesh]{reconstructSSS}},
#'    \code{\link[SurfaceMesh]{alphaWrap}}
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh     <- makeMesh(dataSeptuaginta)
#' mesh_rgl <- toRGL(mesh)
#' mesh_psr <- reconstructPoisson(mesh[["vertices"]],
#'                                normalsFun=getNormalsFun(6L),
#'                                smAngle=20,
#'                                smRadius=1,
#'                                smDistance=20)
#' mesh_psr_rgl <- toRGL(mesh_psr)
#'
#' mfrow3d(1, 2)
#' shade3d(mesh_rgl, col="gray")
#' next3d()
#' shade3d(mesh_psr_rgl, col="gray")
#'
#' @export
reconstructPoisson <- function(
  x,
  normalsFun= getNormalsFun(6L),
  spacing,
  smAngle   = 20,
  smRadius  = 30,
  smDistance= 0.375,
  normals   = FALSE) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L)) {
    stop("The `x` argument must be a numeric matrix with 3 columns.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Missing values in the `x` matrix are not allowed.", call. = TRUE)
  }
  if(nrow(x) <= ncol(x)) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(!is.function(normalsFun) || !inherits(normalsFun, "CGALnormalsFun")) {
    stop("Invalid argument `normalsFun`must be a function as ",
         "returned by `getNormalsFun`")
  }
  if(missing(spacing)) {
    spacing <- -1.0
  } else {
    stopifnot(isPositiveNumber(spacing))
  }
  stopifnot(isPositiveNumber(smAngle))
  stopifnot(isPositiveNumber(smRadius))
  stopifnot(isPositiveNumber(smDistance))
  ## normalsIn matrix does not need to be transposed because it comes vom C++
  normalsIn <- normalsFun(x)
  mesh_cpp <- reconstructPoisson_cpp(
    t(x), normalsIn, spacing, smAngle, smRadius, smDistance, normals)
  fromCPP(mesh_cpp)
}

#' @title Scale-space surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'   A good choice if the input point cloud is noisy but the user still wants the
#'   surface to pass exactly through the points.
#' @param x Numeric matrix with 3 columns which stores the points, one point per row.
#' @param scaleIterations Positive integer. Number of iterations used to increase the scale.
#' @param neighbors Positive integer. Number of neighbors used to smooth the point cloud.
#' @param samples Positive integer. Number of samples used to smooth the point cloud.
#' @param separateShells Boolean. Separate the shells?
#' @param forceManifold Boolean. Force a manifold output mesh?
#' @param borderAngle Bound on the angle in degrees used to detect border edges.
#' @param repairSoup Boolean. Clean the mesh (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices)?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#' @details See \url{https://doc.cgal.org/latest/Scale_space_reconstruction_3/}
#'   for details.
#' @seealso \code{\link[SurfaceMesh]{removeOutliers}},
#'    \code{\link[SurfaceMesh]{reconstructAFS}},
#'    \code{\link[SurfaceMesh]{reconstructPoisson}}, \code{\link[SurfaceMesh]{alphaWrap}}
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' mesh_sss <- reconstructSSS(mesh[["vertices"]],
#'                            scaleIterations=2L,
#'                            forceManifold  =TRUE,
#'                            neighbors      =12L)
#' mesh_sss_rgl <- toRGL(mesh_sss)
#'
#' mfrow3d(1, 2)
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_rgl, col="gray")
#' next3d()
#' view3d(0, -90, zoom=0.7)
#' shade3d(mesh_sss_rgl, col="gray")
#'
#' @export
reconstructSSS <- function(
  x,
  scaleIterations=1L,
  neighbors      =12L,
  samples        =300L,
  separateShells =FALSE,
  forceManifold  =TRUE,
  borderAngle    =45,
  repairSoup     =TRUE,
  normals        =FALSE) {
  if(!is.matrix(x) || !is.numeric(x)) {
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(anyNA(x)) {
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  stopifnot(isStrictPositiveInteger(scaleIterations))
  stopifnot(isStrictPositiveInteger(neighbors), neighbors >= 2)
  stopifnot(isStrictPositiveInteger(samples))
  stopifnot(isBoolean(separateShells))
  stopifnot(isBoolean(forceManifold))
  stopifnot(isNonNegativeNumber(borderAngle))
  stopifnot(isBoolean(repairSoup))
  stopifnot(isBoolean(normals))
  mesh_cpp <- reconstructSSS_cpp(
    t(x),
    as.integer(scaleIterations),
    as.integer(neighbors),
    as.integer(samples),
    separateShells,
    forceManifold,
    as.double(borderAngle),
    repairSoup,
    normals)
  fromCPP(mesh_cpp)
}
