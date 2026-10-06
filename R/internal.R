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

#' @noRd
isFalsy <- function(x){
  isFALSE(x) || is.null(x) || is.na(x)
}

#' @noRd
makeTriangle <- function(vertices, indices) {
  vertices[indices, ]
}

#' @noRd
isAtomicVector <- function(x) {
  is.atomic(x) && is.vector(x)
}

#' @noRd
isNumber <- function(x) {
	is.numeric(x) && (length(x) == 1L) && !is.na(x)
}

#' @noRd
isNumberVector <- function(x) {
  is.vector(x) && !anyNA(x) && is.numeric(x)
}

#' @noRd
isPositiveNumber <- function(x) {
	is.numeric(x) && (length(x) == 1L) && (x > 0) && !is.na(x)
}

#' @noRd
isPositiveNumberVector <- function(x) {
  is.vector(x) && !anyNA(x) && is.numeric(x) && all(x > 0)
}

#' @noRd
isNonNegativeNumber <- function(x){
	is.numeric(x) && (length(x) == 1L) && (x >= 0) && !is.na(x)
}

#' @noRd
isPositiveInteger <- function(x){
	is.numeric(x) && (length(x) == 1L) && !is.na(x) && (floor(x) == x) && (x >= 0)
}

#' @noRd
isStrictPositiveInteger <- function(x){
	isPositiveInteger(x) && (x > 0)
}

#' @noRd
isBoolean <- function(x){
	is.logical(x) && (length(x) == 1L) && !is.na(x)
}

#' @noRd
isString <- function(x){
  is.character(x) && (length(x) == 1L) && !is.na(x)
}

#' @noRd
isStringVector <- function(x) {
  is.vector(x) && is.character(x) && !anyNA(x)
}

#' @noRd
isOrthonormal <- function(x, tol=1e-6) {
  isMat <- is.matrix(x)  &&
           is.numeric(x) &&
           identical(dim(x), c(3L, 3L)) &&
           all(is.finite(x))

  isOn <- isTRUE(all.equal(crossprod(x),
                           diag(3L),
                           check.attributes=FALSE,
                           tolerance=tol))

  isMat && isOn
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @noRd
getVFT <- function(x, beforeCheck = FALSE) {
  transposed <- !beforeCheck
  i0 <- as.integer(transposed)
  if(inherits(x, "mesh3d")) {
    triangles <- x[["it"]]
    if(!is.null(triangles)) {
      triangles <- lapply(seq_len(ncol(triangles)), function(i) { triangles[, i] - i0 })
    }
    quads       <- x[["ib"]]
    mIsTriangle <- is.null(quads)
    mIsQuad     <- is.null(triangles)
    if(!mIsTriangle) {
      quads <- lapply(seq_len(ncol(quads)), function(i) { quads[, i] - i0 })
    }
    ## triangles, quads may be NULL, automatically removed by c()
    faces    <- c(triangles, quads)
    vertices <- x[["vb"]][-4L, ]
    if(!transposed) {
      vertices <- t(vertices)
    }
    rmesh <- list(vertices=vertices, faces=faces)
  } else if(inherits(x, "CGALmesh")) {
    mIsTriangle <- attr(x, "toRGL") == 3L
    mIsQuad     <- attr(x, "toRGL") == 4L
    vertices    <- x[["vertices"]]
    if(transposed) {
      vertices <- t(vertices)
    }
    faces <- x[["faces"]]
    if(is.matrix(faces)) {
      faces <- lapply(seq_len(nrow(faces)), function(i) { faces[i, ] - i0 })
    } else if(!beforeCheck) {
      faces <- lapply(faces, function(face) { face - 1L })
    }
    rmesh <- list(vertices=vertices, faces=faces)
  } else if(is.list(x)) {
    rmesh <- checkMesh(x[["vertices"]], x[["faces"]], aslist=TRUE)
    mIsTriangle <- rmesh[["isTriangle"]]
    mIsQuad     <- rmesh[["isQuad"]]
    if(beforeCheck) {
      rmesh[["vertices"]] <- t(rmesh[["vertices"]])
      rmesh[["faces"]]    <- lapply(rmesh[["faces"]], function(face) { face + 1L })
    }
  } else {
    stop("Invalid `x` argument.", call. = FALSE)
  }
  list(rmesh=rmesh, isTriangle=mIsTriangle, isQuad=mIsQuad)
}

## ----------------------------------------------------------------------- //
## convert R mesh to format required for C++
## ----------------------------------------------------------------------- //
#' @noRd
fromR <- function(x) {
  vertices <- t(x[["vertices"]])
  stopifnot(is.numeric(vertices))
  storage.mode(vertices) <- "double"

  ## create list of faces, reduce index by 1 for C++ counting
  faces <- if(is.matrix(x[["faces"]])) {
    stopifnot(is.numeric(x[["faces"]]))
    n_faces <- nrow(x[["faces"]])
    lapply(seq_len(n_faces), function(i) { as.integer(x[["faces"]][i, ] - 1L) })
  } else {
    stopifnot(all(vapply(x[["faces"]], is.numeric, logical(1))))
    lapply(x[["faces"]], function(face) { as.integer(face - 1L) })
  }

  if(hasName(x, "normals")) {
    normals <- t(x[["normals"]])
    storage.mode(normals) <- "double"
    list(vertices=vertices, faces=faces, normals=normals)
  } else {
    list(vertices=vertices, faces=faces)
  }
}

## ----------------------------------------------------------------------- //
## get matrix of vertices from input: matrix, CGALmesh object, mesh3d object
## ----------------------------------------------------------------------- //
#' @noRd
getVertsMat <- function(x, nPtsMin=4L) {
  if(inherits(x, "CGALmesh")) {
    x[["vertices"]]
  } else if(is.matrix(x)) {
    if(!is.numeric(x) || (ncol(x) != 3L) || (nrow(x) < nPtsMin)) {
      msg <- paste0("`x` must be a numeric matrix with 3 columns and at least ", nPtsMin, " points.")
      stop(msg, call. = TRUE)
    }
    storage.mode(x) <- "double"
    if(anyNA(x)) {
      stop("Points in `x` with missing values are not allowed.", call. = TRUE)
    }
    x
  } else if(inherits(x, "mesh3d")) {
    vft <- getVFT(x, beforeCheck = TRUE)
    vft[["rmesh"]][["vertices"]]
  } else {
    msg <- paste0("The `x` argument must be a 'CGALmesh' object, a 'mesh3d' object,",
                  " or a numeric matrix with 3 columns and at least ", nPtsMin , " points")
    stop(msg)
  }
}

## ----------------------------------------------------------------------- //
## convert CPP mesh to format required in R
## ----------------------------------------------------------------------- //
#' @importFrom utils hasName
#' @noRd
fromCPP <- function(x) {
  x[["vertices"]] <- t(x[["vertices"]])
  if(hasName(x, "edges")) {
    edgesDF                 <- x[["edges"]]
    x[["edgesDF"]]          <- edgesDF
    x[["edges"]]            <- as.matrix(edgesDF[, c("i1", "i2")])
    edgesExterior           <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
    x[["edgesExterior"]]    <- edgesExterior
    x[["verticesExterior"]] <- which(table(edgesExterior) != 2L)
  }

  if(hasName(x, "normals")) {
    x[["normals"]] <- t(x[["normals"]])
  }

  attr(x, "toRGL") <- FALSE
  if(is.matrix(x[["faces"]])) {
    x[["faces"]] <- t(x[["faces"]])
    if(ncol(x[["faces"]]) %in% c(3L, 4L)) {
      attr(x, "toRGL") <- ncol(x[["faces"]])
    }
  } else {
    n_verts   <- lengths(x[["faces"]])
    n_verts_u <- unique(n_verts)
    if(length(n_verts_u) == 1L) {
      x[["faces"]] <- do.call(rbind, x[["faces"]])
      if(n_verts_u %in% c(3L, 4L)) {
        attr(x, "toRGL") <- n_verts_u
      }
    } else {
      if(all(n_verts_u %in% c(3L, 4L))) {
        attr(x, "toRGL") <- 34L
      }
    }
  }

  class(x) <- "CGALmesh"
  x
}

## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
#' @noRd
checkMesh <- function(vertices, faces, aslist) {
  if(!is.matrix(vertices) || (ncol(vertices) != 3L) || !is.numeric(vertices)) {
    stop("The `vertices` argument must be a numeric matrix with 3 columns.")
  }
  storage.mode(vertices) <- "double"
  if(anyNA(vertices)) {
    stop("Found missing values in `vertices`.")
  }

  homFaces    <- FALSE
  mIsTriangle <- FALSE
  mIsQuad     <- FALSE
  mToRGL      <- FALSE
  if(is.matrix(faces)) {
    if(ncol(faces) < 3L) {
      stop("Faces must be given by at least 3 indices.")
    }
    storage.mode(faces) <- "integer"
    if(anyNA(faces)) {
      stop("Found missing values in `faces`.")
    }
    if(any(faces < 1L)) {
      stop("Faces cannot contain indices lower than 1.")
    }
    if(any(faces > nrow(vertices))) {
      stop("Faces cannot contain indices higher than the number of vertices.")
    }

    homFaces <- ncol(faces)
    if(homFaces %in% c(3L, 4L)) {
      mIsTriangle <- homFaces == 3L
      mIsQuad     <- homFaces == 4L
      mToRGL      <- homFaces
    }

    if(aslist) {
      faces <- lapply(seq_len(nrow(faces)), function(i) { faces[i, ] - 1L })
    } else {
      faces <- t(faces - 1L)
    }
  } else if(is.list(faces)) {
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check) {
      stop("The `faces` argument must be a list of integer vectors.")
    }

    check <- any(vapply(faces, anyNA, logical(1L)))
    if(check) {
      stop("Found missing values in `faces`.")
    }

    faces <- lapply(faces, function(x) { as.integer(x) - 1L })
    sizes <- lengths(faces)
    if(any(sizes < 3L)) {
      stop("Faces must be given by at least 3 indices.")
    }

    check <- any(vapply(faces, function(f) {
        any(f < 0L) || any(f >= nrow(vertices)) }, logical(1L)))
    if(check) {
      stop("Faces cannot contain indices lower than 1 or higher ",
          "than the number of vertices.")
    }
    usizes <- length(unique(sizes))
    if(usizes == 1L) {
      homFaces    <- sizes[1L]
      mIsTriangle <- homFaces == 3L
      mIsQuad     <- homFaces == 4L
      if(homFaces %in% c(3L, 4L)) {
        mToRGL <- homFaces
      }
    } else if((usizes == 2L) && all(sizes %in% c(3L, 4L))) {
      mToRGL <- 34L
    }
  } else {
    stop("The `faces` argument must be a list or a matrix.")
  }
  list("vertices"  =t(vertices),
       "faces"     =faces,
       "homFaces"  =homFaces,
       "isTriangle"=mIsTriangle,
       "isQuad"    =mIsQuad,
       "toRGL"     =mToRGL)
}

## ----------------------------------------------------------------------- //
# vertices = numeric matrix with 3 columns
# faces = integer matrix or list
# no missings
## ----------------------------------------------------------------------- //
#' @noRd
checkMeshValid <- function(vertices, faces, aslist) {
  if(!is.matrix(vertices) || (ncol(vertices) != 3L) || !is.numeric(vertices)) {
    stop("The `vertices` argument must be a numeric matrix with three columns.")
  }
  storage.mode(vertices) <- "double"

  homFaces    <- FALSE
  mIsTriangle <- FALSE
  mIsQuad     <- FALSE
  mToRGL      <- FALSE
  if(is.matrix(faces)) {
    storage.mode(faces) <- "integer"
    homFaces <- ncol(faces)
    if(homFaces %in% c(3L, 4L)) {
      mIsTriangle <- homFaces == 3L
      mIsQuad     <- homFaces == 4L
      mToRGL      <- homFaces
    }

    if(aslist) {
      faces <- lapply(seq_len(nrow(faces)), function(i) { faces[i, ] - 1L })
    } else {
      faces <- t(faces - 1L)
    }
  } else if(is.list(faces)) {
    faces  <- lapply(faces, function(x) { as.integer(x) - 1L })
    sizes  <- lengths(faces)
    usizes <- length(unique(sizes))
    if(usizes == 1L) {
      homFaces    <- sizes[1L]
      mIsTriangle <- homFaces == 3L
      mIsQuad     <- homFaces == 4L
      if(homFaces %in% c(3L, 4L)) {
        mToRGL <- homFaces
      }
    } else if((usizes == 2L) && all(sizes %in% c(3L, 4L))) {
      mToRGL <- 34L
    }
  } else {
    stop("The `faces` argument must be a list or a matrix.")
  }
  list("vertices"  =t(vertices),
       "faces"     =faces,
       "homFaces"  =homFaces,
       "isTriangle"=mIsTriangle,
       "isQuad"    =mIsQuad,
       "toRGL"     =mToRGL)
}

## ----------------------------------------------------------------------- //
## options for point sampling on mesh
## see internal_utils.cpp -> CGAL sample_triangle_mesh()
## ----------------------------------------------------------------------- //
#' @noRd
checkSampleOpts <- function(x) {
  if(!hasName(x, "method") || is.null(x[["method"]])) {
    x[["method"]]  <- 1L
  } else {
    method_choices <- c("random", "grid", "mc")
    method         <- match.arg(x[["method"]], choices=method_choices)
    x[["method"]]  <- match(method, method_choices)
  }
  if(!hasName(x, "sampleVerts") || is.null(x[["sampleVerts"]])) {
    x[["sampleVerts"]] <- TRUE
  } else {
    stopifnot(isBoolean(x[["sampleVerts"]]))
  }
  if(!hasName(x, "sampleEdges") || is.null(x[["sampleEdges"]])) {
    x[["sampleEdges"]] <- TRUE
  } else {
    stopifnot(isBoolean(x[["sampleEdges"]]))
  }
  if(!hasName(x, "sampleFaces") || is.null(x[["sampleFaces"]])) {
    x[["sampleFaces"]] <- TRUE
  } else {
    stopifnot(isBoolean(x[["sampleFaces"]]))
  }
  if(!hasName(x, "ptsOnEdges") || is.null(x[["ptsOnEdges"]])) {
    x[["ptsOnEdges"]] <- 0L
  } else {
    stopifnot(isStrictPositiveInteger(x[["ptsOnEdges"]]))
  }
  if(!hasName(x, "ptsOnFaces") || is.null(x[["ptsOnFaces"]])) {
    x[["ptsOnFaces"]] <- 0L
  } else {
    stopifnot(isStrictPositiveInteger(x[["ptsOnFaces"]]))
  }
  if(!hasName(x, "gridSpacing") || is.null(x[["gridSpacing"]])) {
    x[["gridSpacing"]] <- 0.0
  } else {
    stopifnot(isPositiveNumber(x[["gridSpacing"]]))
  }
  if(!hasName(x, "ptsPerDist") || is.null(x[["ptsPerDist"]])) {
    x[["ptsPerDist"]] <- 0.0
  } else {
    stopifnot(isPositiveNumber(x[["ptsPerDist"]]))
  }
  if(!hasName(x, "ptsPerEdge") || is.null(x[["ptsPerEdge"]])) {
    x[["ptsPerEdge"]] <- 0L
  } else {
    stopifnot(isStrictPositiveInteger(x[["ptsPerEdge"]]))
  }
  if(!hasName(x, "ptsPerArea") || is.null(x[["ptsPerArea"]])) {
    x[["ptsPerArea"]] <- 0.0
  } else {
    stopifnot(isPositiveNumber(x[["ptsPerArea"]]))
  }
  if(!hasName(x, "ptsPerFace") || is.null(x[["ptsPerFace"]])) {
    x[["ptsPerFace"]] <- 0L
  } else {
    stopifnot(isStrictPositiveInteger(x[["ptsPerFace"]]))
  }

  x[["method"]]     <- as.integer(x[["method"]])
  x[["ptsOnEdges"]] <- as.integer(x[["ptsOnEdges"]])
  x[["ptsOnFaces"]] <- as.integer(x[["ptsOnFaces"]])
  x[["ptsPerEdge"]] <- as.integer(x[["ptsPerEdge"]])
  x[["ptsPerFace"]] <- as.integer(x[["ptsPerFace"]])

  storage.mode(x[["gridSpacing"]]) <- "double"
  storage.mode(x[["ptsPerDist"]])  <- "double"
  storage.mode(x[["ptsPerArea"]])  <- "double"

  x
}

## ----------------------------------------------------------------------- //
## Euler angles from 3D rotation matrix
## ----------------------------------------------------------------------- //
#' @noRd
getEulerAngles <- function(x, tol=1e-6) {
    stopifnot(isOrthonormal(x, tol=tol))
    ## if reflection: flip 1 axis
    if(det(x) < 0) {
        x[, 3L] <- -x[, 3L]
    }
    theta <- acos(max(-1, min(1, x[3L, 3L])))
    if(abs(sin(theta)) > 1e-8) {
        phi <- atan2(x[2L, 3L],  x[1L, 3L])
        psi <- atan2(x[3L, 2L], -x[3L, 1L])
    } else if(x[3L, 3L] > 0) {      # gimbal lock: only phi -/+ psi is determined
        phi <- atan2(x[2L, 1L], x[1L, 1L])
        psi <- 0
    } else {
        phi <- atan2(-x[2L, 1L], -x[1L, 1L])
        psi <- 0
    }
    c(phi, theta, psi)
}
