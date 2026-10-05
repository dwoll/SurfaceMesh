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
#' @param x A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}},
#'   a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   or a numeric matrix with 3 columns which stores the point coordinates,
#'   one point per row, and at least 4 points.
#' @param oriented Boolean. Get the optimal (oriented) bounding box?
#' @param out \code{character}. Type of return value. One of \code{"Points"}
#'   for a matrix giving the corner points, \code{"CGALmesh"} for a \code{CGALmesh}
#'   object, and \code{"rgl"} for a \code{\link[rgl]{mesh3d}} object from
#'   package \strong{rgl}.
#' @param triangulate Boolean. Triangulate the faces of the bounding box?
#' @param normals Boolean. Return vertex normals?
#' @returns A \code{CGALmesh} object, a \code{\link[rgl]{mesh3d}} object,
#'   or a matrix with the corner points as rows.
#' @seealso See \code{\link[SurfaceMesh]{getConvexHull}} for the convex hull,
#'   \code{\link[SurfaceMesh]{getBoundingSphere}} for the bounding sphere,
#'   \code{\link[SurfaceMesh]{getBoundingEll}} for the bounding ellipse.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#'
#' mesh1     <- makeMesh(dataTubularKnot)
#' mesh1_rgl <- toRGL(mesh1)
#' bb1_rgl   <- getBoundingBox(mesh1[["vertices"]], out="rgl")
#' open3d(windowRect=50 + c(0, 0, 400, 400))
#' view3d(0, 30, zoom=0.9)
#' shade3d(mesh1_rgl, color="seashell")
#' wire3d(bb1_rgl)
#'
#' mesh2     <- dataHeart1
#' mesh2_rgl <- toRGL(mesh2)
#' bb2_rgl   <- getBoundingBox(mesh2[["vertices"]], oriented=TRUE, out="rgl")
#' open3d(windowRect=50 + c(0, 0, 400, 400))
#' shade3d(mesh2_rgl, color="seashell")
#' wire3d(bb2_rgl)
#'
#' @export
#' @importFrom rgl translate3d scale3d cube3d
getBoundingBox <- function(x,
                           oriented    = FALSE,
                           out         = c("CGALmesh", "rgl", "Points"),
                           triangulate = FALSE,
                           normals     = FALSE) {
  xIn <- getVertsMat(x, nPtsMin=4L)
  stopifnot(isBoolean(oriented))
  out_choices <- tolower(c("CGALmesh", "rgl", "Points"))
  out         <- match.arg(tolower(out), choices=out_choices)
  stopifnot(isBoolean(triangulate))
  stopifnot(isBoolean(normals))
  if(oriented) {
    outL <- getBoundingBoxOptimal_cpp(t(xIn), triangulate, normals)
    if(out == "points") {
      t(outL[["vertices"]])
    } else {
      obb <- fromCPP(outL[["mesh"]])
      if(out == "cgalmesh") {
        obb
      } else if(out == "rgl") {
        toRGL(obb)
      } else {
        stop("Wrong output format.")
      }
    }
  } else {
    outL <- getBoundingBox_cpp(t(xIn))
    ptLo <- outL[["lo"]]
    ptUp <- outL[["up"]]
    if(out == "points") {
      ax <- ptLo[1L]
      ay <- ptLo[2L]
      az <- ptLo[3L]
      bx <- ptUp[1L]
      by <- ptUp[2L]
      bz <- ptUp[3L]
      rbind(c(ax, ay, az),
            c(bx, ay, az),
            c(bx, ay, bz),
            c(ax, ax, bz),
            c(ax, by, az),
            c(bx, by, az),
            c(bx, by, bz),
            c(ax, by, bz))
    } else {
      m_rgl <- meshIsoCuboid(ptLo, ptUp)
      if(out == "cgalmesh") {
        makeMesh(m_rgl, repairSoup=FALSE, triangulate=triangulate, normals=normals)
      } else if(out == "rgl") {
        m_rgl
      } else {
        stop("Wrong output format.")
      }
    }
  }
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Bounding ellipsoid of a set of 3D points
#' @description Get the approximate bounding ellipsoid of a given set of 3D points.
#' @param x A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}},
#'   a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   or a numeric matrix with 3 columns which stores the point coordinates,
#'   one point per row, and at least 4 points.
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
#'   and \code{"directions"};
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
                           out = c("CGALmesh", "rgl", "CtrRadDir"),
                           nIter = 3L,
                           eps = 0.01,
                           normals = FALSE) {
  xIn         <- getVertsMat(x, nPtsMin=4L)
  out_choices <- tolower(c("CGALmesh", "rgl", "CtrRadDir"))
  out         <- match.arg(tolower(out), choices=out_choices)
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isPositiveNumber(eps))
  storage.mode(eps) <- "double"
  stopifnot(isBoolean(normals))

  bellL <- getBoundingEllipsoid_cpp(t(xIn), eps)
  if(out == "ctrraddir") {
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
    if(out == "cgalmesh") {
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
#' @param x A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}},
#'   a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   or a numeric matrix with 3 columns which stores the point coordinates,
#'   one point per row, and at least 4 points.
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
                              out = c("CGALmesh", "rgl", "CtrRad"),
                              nIter = 3L,
                              normals = FALSE) {
  xIn         <- getVertsMat(x, nPtsMin=4L)
  out_choices <- tolower(c("CGALmesh", "rgl", "CtrRad"))
  out         <- match.arg(tolower(out), choices=out_choices)
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(normals))

  bsL <- getBoundingSphere_cpp(t(xIn))
  if(out == "ctrrad") {
    bsL
  } else {
    bs_rgl <- meshSphere(bsL[["center"]], r=bsL[["radius"]], nIter=nIter)
    if(out == "cgalmesh") {
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
  stopifnot(isBoolean(normals))
  xIn     <- getVertsMat(x, nPtsMin=4L)
  meshCPP <- getConvexHull_cpp(t(xIn), normals)
  fromCPP(meshCPP)
}
