#' @title Get the face colors of a mesh
#' @description Get the face colors of a 3D surface mesh
#'   Note: This function is currently doing anything since face colors
#'   that may be present in the R input object are not imported to
#'   the C++ side in \code{\link[SurfaceMesh]{makeMesh}}. In the
#'   future, this function may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @return The matrix of colors (integer RGB values) attached to
#'   the faces of the mesh. \code{NA_integer_} if there are no
#'   face colors.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' mesh <- makeMeshValid(dataSphere)
#' getFaceColors(mesh)
#'
#' @export
getFaceColors <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
          " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    getFaceColors_cpp(meshCPP)
}

#' @title Get the vertex colors of a mesh
#' @description Get the vertex colors of a 3D surface mesh
#'   Note: This function is currently doing anything since vertex colors
#'   that may be present in the R input object are not imported to
#'   the C++ side in \code{\link[SurfaceMesh]{makeMesh}}. In the
#'   future, this function may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @return The matrix of colors (integer RGB values) attached to
#'   the vertices of the mesh. \code{NA_integer_} if there are no
#'   vertex colors.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' mesh <- makeMeshValid(dataSphere)
#' getVertexColors(mesh)
#'
#' @export
getVertexColors <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
          " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    getVertexColors_cpp(meshCPP)
}

#' @title Get the vertex normals of a mesh
#' @description Get the vertex normals of a 3D surface mesh
#'   Note: This function is currently doing anything since vertex normals
#'   that may be present in the R input object are not imported to
#'   the C++ side in \code{\link[SurfaceMesh]{makeMesh}}. In the
#'   future, this function may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @return The matrix of normals attached to
#'   the vertices of the mesh. \code{NA_integer_} if there are no
#'   vertex normals
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
#' @title Assign colors to mesh faces
#' @description Assign given colors to the faces of a 3D surface mesh.
#'   Note: This function is currently doing anything except adding a list
#'   component \code{normals} to the \code{CGALmesh} object after setting
#'   the property map on the C++ side, and then exporting back to R. In the
#'   future, this may be useful.
#' @param colors Either a \code{character} vector with color names or an
#'   integer matrix with 3 columns for RGB values in [0, 255]. When the
#'   vector has a single color name, or when the matrix has a single row,
#'   this color is assigned to all faces. Otherwise, the vector must have
#'   as many elements as there are faces, or the matrix must have as many
#'   rows as there are faces.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(SurfaceMesh)
#' mesh     <- makeMeshValid(dataTorus)
#' mesh_wfc <- setFaceColors(mesh, "navyblue")
#' str(mesh_wfc)
#' @export
setFaceColors <- function(x, colors) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
             " (i.e., the output of the `makeMesh()` function).")
    }
    colorsMat <- if(isStringVector(colors)) {
        col2rgb(colors, alpha=FALSE) # 3 rows -> already transposed
    } else if(is.matrix(colors)) {
        # matrix with 3 columns R, G, B
        stopifnot(ncol(colors) == 3L, is.numeric(colors))
        mode(colors) <- "integer"
        stopifnot(all(colors >= 0L), all(colors <= 255L))
        t(colors)
    } else {
        # vector with 3 values R, G, B
        stopifnot(length(colors) == 3L, is.numeric(colors))
        mode(colors) <- "integer"
        stopifnot(all(colors >= 0L), all(colors <= 255L))
        as.matrix(colors)  # column vector
    }
    meshCPP <- fromR(x)
    meshOut <- setFaceColors_cpp(meshCPP, colorsMat)
    fromCPP(meshCPP)
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Assign colors to mesh vertices
#' @description Assign given colors to the vertices of a 3D surface mesh.
#'   Note: This function is currently doing anything except adding a list
#'   component \code{normals} to the \code{CGALmesh} object after setting
#'   the property map on the C++ side, and then exporting back to R. In the
#'   future, this may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param colors Either a \code{character} vector with color names or an
#'   integer matrix with 3 columns for RGB values in from 0 to 255. The vector
#'   must have as many elements as there are faces, or the matrix must have as
#'   many rows as there are faces.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' ## TODO
#' library(SurfaceMesh)
#' @export
setVertexColors <- function(x, colors) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  colorsMat <- if(isStringVector(colors)) {
      col2rgb(colors, alpha=FALSE) # 3 rows -> already transposed
  } else if(is.matrix(colors)) {
      # matrix with 3 columns R, G, B
      stopifnot(ncol(colors) == 3L, is.numeric(colors))
      mode(colors) <- "integer"
      stopifnot(all(colors >= 0L), all(colors <= 255L))
      t(colors)
  } else {
      stop("Wrong format for `colors`.")
  }
  meshCPP <- fromR(x)
  meshOut <- setVertexColors_cpp(meshCPP, colorsMat)
  fromCPP(meshCPP)
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @title Assign given normal vectors to mesh vertices
#' @description Assign given per-vertex normal vectors to a 3D surface mesh.
#'   Note: This function is currently doing anything except adding a list
#'   component \code{normals} to the \code{CGALmesh} object after setting
#'   the property map on the C++ side, and then exporting back to R. In the
#'   future, this may be useful.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[SurfaceMesh]{makeMesh}}.
#' @param normals A numeric matrix with three columns and as many rows as
#'   the number of vertices.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' ## TODO
#' library(SurfaceMesh)
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
