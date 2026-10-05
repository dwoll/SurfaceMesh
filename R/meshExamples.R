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

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Ellipsoid mesh
#' @description Create a triangle mesh of a 3D ellipsoid.
#' @param ctr \code{numeric}. Vector with the x, y, z-coordinates of the centroid.
#' @param r \code{numeric}. Vector with the 3 semi-axis lengths.
#' @param dirs \code{numeric}. 3-by-3 matrix. Columns are the unit
#'   direction vectors of the ellipsoid's semi-axes.
#'   CAVE: Y-axis direction is assumed to be front-back, corresponding to
#'   \strong{rgl} (OpenGL) convention.
#' @param euler \code{numeric}. Euler rotation angles phi, theta, psi. If
#'  given, takes precedence over \code{dirs}.
#' @param nIter \code{integer}. Number of iterations (the mesh is obtained
#'   by iteratively subdividing the faces of an icosahedron).
#' @returns A \code{\link[rgl]{mesh3d}} object.
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#' mesh_rgl <- meshEllipsoid(c(0, 0, 0), r=c(1, 0.5, 0.25))
#' wire3d(mesh_rgl)
#'
#' @importFrom rgl subdivision3d cube3d scale3d rotate3d translate3d
#' @export
meshEllipsoid <- function(ctr = c(0, 0, 0),
                          r = c(1, 1, 1),
                          dirs = diag(nrow=3L),
                          euler,
                          nIter = 3L) {
  stopifnot(isNumberVector(ctr))
  stopifnot(isPositiveNumberVector(r))
  stopifnot(isStrictPositiveInteger(nIter))
  sphere <- rgl::subdivision3d(rgl::icosahedron3d(), depth = nIter)
  class(sphere) <- c("mesh3d", "shape3d")

  verts       <- sphere[["vb"]][1L:3L, ]
  norm        <- sqrt(colSums(verts^2))
  verts[1L, ] <- verts[1L, ] / norm
  verts[2L, ] <- verts[2L, ] / norm
  verts[3L, ] <- verts[3L, ] / norm

  sphere[["normals"]] <- verts
  sphere[["vb"]]      <- rbind(verts, 1)
  ellipsoid <- rgl::scale3d(sphere, r[1L], r[2L], r[3L])
  eulerUse  <- if(!missing(euler)) {
    stopifnot(isNumberVector(euler))
    euler
  } else {
    ## dirs matrix is checked in getEulerAngles()
    getEulerAngles(dirs)
  }

  ## see ?rotate3d() for convention
  rot_mat <- t(rgl::rotationMatrix(eulerUse[1], 0, 0, 1) %*%
               rgl::rotationMatrix(eulerUse[2], 0, 1, 0) %*%
               rgl::rotationMatrix(eulerUse[3], 1, 0, 0))
  ellipsoid <- rgl::rotate3d(ellipsoid, matrix=rot_mat)
  rgl::translate3d(ellipsoid, ctr[1L], ctr[2L], ctr[3L])
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Sphere mesh
#' @description Triangle mesh of a 3D sphere.
#' @param ctr \code{numeric}. Vector with the x, y, z-coordinates of the centroid.
#' @param r \code{numeric}. Radius.
#' @param nIter \code{integer}. Number of iterations (the mesh is obtained
#'   by iteratively subdividing the faces of an icosahedron).
#' @returns A \code{\link[rgl]{mesh3d}} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#' mesh_rgl <- meshSphere(c(0, 0, 0), r=2)
#' wire3d(mesh_rgl)
#'
#' @export
#' @importFrom rgl subdivision3d icosahedron3d scale3d translate3d
meshSphere <- function(ctr = c(0, 0, 0), r = 1, nIter = 3L) {
    stopifnot(isNumberVector(ctr))
    stopifnot(isPositiveNumber(r))
    stopifnot(isStrictPositiveInteger(nIter))
    sphere <- rgl::subdivision3d(rgl::icosahedron3d(), depth = nIter)
    class(sphere) <- c("mesh3d", "shape3d")

    verts       <- sphere[["vb"]][1L:3L, ]
    norm        <- sqrt(colSums(verts^2))
    verts[1L, ] <- verts[1L, ] / norm
    verts[2L, ] <- verts[2L, ] / norm
    verts[3L, ] <- verts[3L, ] / norm

    sphere[["normals"]] <- verts
    sphere[["vb"]]      <- rbind(verts, 1)
    rgl::translate3d(rgl::scale3d(sphere, r, r, r),
                     ctr[1], ctr[2], ctr[3])
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Iso-oriented cuboid
#' @description Mesh of an iso-oriented cuboid, i.e. with edges parallel
#'   to the axes.
#'
#' @param lo \code{numeric} 3-vector. Lower left front corner, a point
#'   whose coordinates must all be lower than those of the upper corner.
#' @param up \code{numeric} 3-vector. Upper right back corner, a point
#'   whose coordinates must all be higher than those of the lower corner.
#' @details CAVE: Y-axis direction is assumed to be front-back,
#'   corresponding to \strong{rgl} (OpenGL) convention.
#'
#' @return A \strong{rgl} mesh, i.e. a \code{mesh3d} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @exmaples
#' library(SurfaceMesh)
#' library(rgl)
#' lo <- c(-4, -2, -1)
#' up <- c( 4,  2,  1)
#' cuboid_rgl <- meshIsoCuboid(lo, up)
#' wire3d(cuboid_rgl)
#'
#' @importFrom rgl cube3d translate3d scale3d
#' @export
meshIsoCuboid <- function(lo, up) {
  stopifnot(all(lo <= up))
  vd  <- up - lo
  ctr <- lo + 0.5*vd
  ax  <- vd[1L]
  ay  <- vd[2L]
  az  <- vd[3L]
  rgl::translate3d(rgl::scale3d(rgl::cube3d(), ax/2, ay/2, az/2),
                   ctr[1L], ctr[2L], ctr[3L])
}
