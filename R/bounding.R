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
#' @title Get bounding box (axis-parallel or oriented)
#' @description Get the axis-parallel or optimal (oriented) bounding box of a 3D surface mesh.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param oriented Boolean. Get the optimal (oriented) bounding box?
#' @param triangulate Boolean. Triangulate the faces of the bounding box?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @seealso See \code{\link[SurfaceMesh]{getConvexHull}} for the convex hull,
#'   \code{\link[SurfaceMesh]{getBoundingSphere}} for the bounding sphere,
#'   \code{\link[SurfaceMesh]{getBoundingEll}} for the bounding ellipse.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh1     <- dataOloid
#' mesh1_rgl <- toRGL(mesh1)
#' bb1       <- getBoundingBox(mesh1)
#' bb1_rgl   <- toRGL(bb1)
#' view3d(0, 30, zoom=0.9)
#' wire3d(mesh1_rgl)
#' wire3d(bb1_rgl)
#'
#' mesh2     <- dataHeart1
#' mesh2_rgl <- toRGL(mesh2)
#' bb2       <- getBoundingBox(mesh2, oriented=TRUE)
#' bb2_rgl   <- toRGL(bb2)
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh2_rgl)
#' wire3d(bb2_rgl)

#' @export
#' @importFrom rgl translate3d scale3d cube3d
getBoundingBox <- function(x,
                           oriented    = FALSE,
                           triangulate = FALSE,
                           normals     = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(oriented))
  stopifnot(isBoolean(triangulate))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- if(oriented) {
    outL <- getBoundingBoxOptimal_cpp(meshCPP, triangulate, normals)
    fromCPP(outL[["mesh"]])
  } else {
    outL    <- getBoundingBox_cpp(meshCPP)
    lcorner <- outL[["lcorner"]]
    ucorner <- outL[["ucorner"]]
    center  <- (lcorner + ucorner) / 2
    ax      <- ucorner[1L] - lcorner[1L]
    ay      <- ucorner[2L] - lcorner[2L]
    az      <- ucorner[3L] - lcorner[3L]
    m_rgl   <- rgl::translate3d(rgl::scale3d(rgl::cube3d(), ax/2, ay/2, az/2),
                                center[1L], center[2L], center[3L])

    makeMesh(m_rgl, repairSoup=FALSE, triangulate=triangulate, normals=normals)
  }
  meshOut
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Bounding ellipsoid of a set of 3D points
#' @description Get the approximate bounding ellipsoid of a given set of 3D points.
#' @param x \code{numeric} matrix with 3 columns with one point per row.
#' @param out \code{character}. Type of return value. One of \code{"CtrRadDir"}
#'   for a list giving the sphere's center, semi-axis lengths, and semi-axis
#'   directions, \code{"CGALmesh"} for a \code{CGALmesh} object, and
#'   \code{"rgl"} for a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#' @param eps \code{numeric}. For the approximation ratio \code{1 + eps}.
#'   See details.
#' @param nIter \code{integer}. Number of iterations (the mesh is obtained
#'   by iteratively subdividing the faces of an icosahedron).
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{list} with components \code{"center"}, \code{"lengths"},
#'   and \code{directions};
#'   a \code{CGALmesh} object; or a \code{\link[rgl]{mesh3d}} object -
#'   depending on argument \code{out}.
#' @seealso See \code{\link[SurfaceMesh]{getBoundingBox}} for the bounding box,
#'   \code{\link[SurfaceMesh]{getBoundingSphere}} for the bounding sphere,
#'   \code{\link[SurfaceMesh]{getConvexHull}} for the convex hull.
#' @details See \url{https://doc.cgal.org/latest/Bounding_volumes/} for details.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh     <- makeMesh(dataHeart1)
#' mesh_rgl <- toRGL(mesh)
#' bell_rgl <- getBoundingEll(mesh[["vertices"]], out="rgl", eps=0.01)
#'
#' shade3d(mesh_rgl, color="lightgreen")
#' wire3d(bell_rgl)
#'
#' @export
getBoundingEll <- function(x,
                           out = c("CtrRadDir", "CGALmesh", "rgl"),
                           nIter = 3L,
                           eps = 0.01,
                           normals = FALSE) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L) || (nrow(x) <= 3L)) {
    stop("`x` must be a numeric matrix with 3 columns and at least 3 points.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Points in `x` with missing values are not allowed.", call. = TRUE)
  }
  out_choices <- tolower(c("CtrRadDir", "CGALmesh", "rgl"))
  out         <- match.arg(tolower(out), choices=out_choices)
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isPositiveNumber(eps))
  storage.mode(eps) <- "double"
  stopifnot(isBoolean(normals))

  bellL <- getBoundingEllipsoid_cpp(t(x), eps)
  if(out == "CtrRadDir") {
    bellL
  } else {
    ctr  <- bellL[["center"]]
    r    <- bellL[["lengths"]]
    dirs <- bellL[["directions"]]

    ## CAVE: input from CGAL y is up-down,
    ## but rgl y is front-back
    dirsUse        <- dirs
    dirsUse[ , 2L] <- dirs[ , 3L]
    dirsUse[ , 3L] <- dirs[ , 2L]

    bell_rgl <- meshEllipsoid(ctr, r=r, dirs=dirsUse, nIter=nIter)
    if(out == "CGALmesh") {
      makeMeshValid(bell_rgl, normals=normals)
    } else if(out == "rgl") {
      if(!normals) {
        bell_rgl[["normals"]] <- NULL
      }
      bell_rgl
    } else {
      stop("Wrong output format.")
    }
  }
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Bounding sphere of a set of 3D points
#' @description Get the bounding sphere of a given set of 3D points.
#' @param x \code{numeric} matrix with 3 columns with one point per row.
#' @param out \code{character}. Type of return value. One of \code{"CtrRad"}
#'   for a list giving the sphere's center and radius, \code{"CGALmesh"} for
#'   a \code{CGALmesh} object, and \code{"rgl"} for a \code{\link[rgl]{mesh3d}}
#'   object from package \strong{rgl}.
#' @param nIter \code{integer}. Number of iterations (the mesh is obtained
#'   by iteratively subdividing the faces of an icosahedron)
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{list} with components \code{"center"} and \code{"radius"},
#'   a \code{CGALmesh} object; or a \code{\link[rgl]{mesh3d}} object -
#'   depending on argument \code{out}.
#' @seealso See \code{\link[SurfaceMesh]{getBoundingBox}} for the bounding box,
#'   \code{\link[SurfaceMesh]{getBoundingEll}} for the bounding ellipse,
#'   \code{\link[SurfaceMesh]{getConvexHull}} for the convex hull.
#' @details See \url{https://doc.cgal.org/latest/Bounding_volumes/} for details.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh     <- makeMesh(dataToroHelix)
#' mesh_rgl <- toRGL(mesh)
#' bs_rgl   <- getBoundingSphere(mesh[["vertices"]], out="rgl")
#'
#' shade3d(mesh_rgl, color="seashell")
#' wire3d(bs_rgl)
#'
#' @export
getBoundingSphere <- function(x,
                              out = c("CtrRad", "CGALmesh", "rgl"),
                              nIter = 3L,
                              normals = FALSE) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L) || (nrow(x) <= 3L)) {
    stop("`x` must be a numeric matrix with 3 columns and at least 3 points.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Points in `x` with missing values are not allowed.", call. = TRUE)
  }
  out_choices <- tolower(c("CtrRad", "CGALmesh", "rgl"))
  out         <- match.arg(tolower(out), choices=out_choices)
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(normals))

  bsL <- getBoundingSphere_cpp(t(x))
  if(out == "CtrRad") {
    bsL
  } else {
    bs_rgl <- meshSphere(bsL[["center"]], r=bsL[["radius"]], nIter=nIter)
    if(out == "CGALmesh") {
      makeMeshValid(bs_rgl, normals=normals)
    } else if(out == "rgl") {
      if(!normals) {
        bs_rgl[["normals"]] <- NULL
      }
      bs_rgl
    } else {
      stop("Wrong output format.")
    }
  }
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Convex hull of a set of 3D points
#' @description Get the convex hull of a given set of 3D points
#'   using the quickhull algorithm.
#' @param x \code{numeric} matrix with 3 columns with one point per row.
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object.
#' @seealso See \code{\link[SurfaceMesh]{getBoundingBox}} for the bounding box,
#'   \code{\link[SurfaceMesh]{getBoundingSphere}} for the bounding sphere,
#'   \code{\link[SurfaceMesh]{getBoundingEll}} for the bounding ellipse.
#' @details See \url{https://doc.cgal.org/latest/Convex_hull_3/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh     <- makeMesh(dataHopfTorus)
#' mesh_rgl <- toRGL(mesh)
#' hull     <- getConvexHull(mesh[["vertices"]])
#' hull_rgl <- toRGL(hull)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(hull_rgl)
#'
#' @export
getConvexHull <- function(x, normals = FALSE) {
  if(!is.matrix(x) || !is.numeric(x) || (ncol(x) != 3L) || (nrow(x) <= 3L)) {
    stop("`x` must be a numeric matrix with 3 columns and at least 3 points.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)) {
    stop("Points in `x` with missing values are not allowed.", call. = TRUE)
  }
  stopifnot(isBoolean(normals))
  meshCPP <- getConvexHull_cpp(t(x), normals)
  fromCPP(meshCPP)
}
