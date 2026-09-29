## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Get the vertex normals of a mesh
#' @description Get the vertex normals of a 3D surface mesh
#'   Note: This function is currently not doing anything since vertex
#'   normals that may be present in the R input object are not imported
#'   to the C++ side in \code{\link[SurfaceMesh]{makeMesh}}.
#'   In the future, this function may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @return The matrix of normals attached to the vertices of the mesh.
#'   \code{NULL} if there are no vertex normals.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' mesh <- makeMeshValid(dataSphere, normals=TRUE)
#' getVertexNormals(mesh)
#'
#' @export
getVertexNormals <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
          " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    getVertexNormals_cpp(meshCPP)
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Assign given normal vectors to mesh vertices
#' @description Assign given per-vertex normal vectors to a 3D surface mesh.
#'   Note: This function is currently doing anything except adding a list
#'   component \code{normals} to the \code{CGALmesh} object after setting
#'   the property map on the C++ side, and then exporting it back to R.
#'   In the future, this may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param normals A numeric matrix with three columns and as many rows as
#'   the number of vertices.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE, normals=FALSE)
#'
#' ## random unit normal vectors
#' nVerts  <- nrow(mesh[["vertices"]])
#' nrmls0  <- matrix(runif(nVerts*3), ncol=3)
#' lens    <- sqrt(diag(tcrossprod(nrmls0)))
#' nrmls   <- diag(1/lens) %*% nrmls0
#' mesh_vn <- setVertexNormals(mesh, nrmls)
#' head(mesh_vn[["normals"]])
#' head(nrmls)
#' @export
setVertexNormals <- function(x, normals) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(is.matrix(normals), is.numeric(normals), ncol(normals) == 3L)
  storage.mode(normals) <- "double"
  if(anyNA(normals)) {
    stop("Vectors in `normals` with missing values are not allowed.", call. = TRUE)
  }
  meshCPP <- fromR(x)
  meshOut <- setVertexNormals_cpp(meshCPP, t(normals))
  fromCPP(meshOut)
}
