## ----------------------------------------------------------------------- //
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of a human heart from a CT scan
#'
#' @description An object of class \code{CGALmesh} representing a human heart
#'     based on a CT scan.
#'
#' @format A \code{CGALmesh} object.
#' @examples
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(-5, -90, zoom=0.8)
#' shade3d(toRGL(dataHeart2), color="gray")
"dataHeart2"
