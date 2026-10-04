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
#'   direction vectors of the ellipsoid's semi-axes. Y-axis is assumed to be
#'   front-back, according to rgl (OpenGL) convention.
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
  eulerUse <- if(!missing(euler)) {
    stopifnot(isNumberVector(euler))
    euler
  } else {
    stopifnot(isOrthonormal(dirs, tol=1e-6))
    ## if reflection: flip 1 axis
    if(isTRUE(all.equal(det(dirs1), -1, tol=0.0001))) {
      dirs1[ , 1] <- -1*dirs1[ , 1]
    }

    getEulerAngles(dirs1)
  }

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

    ## sphere[["vb"]] <- rbind(r*verts + c(ctr[1], ctr[2], ctr[3]), 1)
    sphere[["normals"]] <- verts
    sphere[["vb"]] <- rbind(verts, 1)
    sphere <- rgl::scale3d(sphere, r, r, r)
    rgl::translate3d(sphere, ctr[1], ctr[2], ctr[3])
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Iso-oriented cuboid
#' @description Mesh of an iso-oriented cuboid, i.e. a cuboid with
#'   edges parallel to the axes.
#'
#' @param lcorner lower corner, a point whose coordinates must be
#'   lower than those of the upper corner.
#' @param ucorner upper corner, a point whose coordinates must be
#'   greater than those of the lower corner.
#'
#' @return A \strong{rgl} mesh, i.e. a \code{mesh3d} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @importFrom rgl cube3d translate3d scale3d
#' @export
meshIsoCuboid <- function(lcorner, ucorner) {
  stopifnot(all(lcorner <= ucorner))
  center <- (lcorner + ucorner) / 2
  ax <- ucorner[1L] - lcorner[1L]
  ay <- ucorner[2L] - lcorner[2L]
  az <- ucorner[3L] - lcorner[3L]
  rgl::translate3d(rgl::scale3d(rgl::cube3d(), ax/2, ay/2, az/2),
                   center[1L], center[2L], center[3L])
}
