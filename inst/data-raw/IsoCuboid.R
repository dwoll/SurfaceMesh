## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## Iso cuboid
## ----------------------------------------------------------------------- //

#' @title Iso-oriented cuboid
#' @description Mesh of an iso-oriented cuboid, i.e. a cuboid with 
#'   edges parallel to the axes.
#'
#' @param lcorner lower corner, a point whose coordinates must be 
#'   lower than those of the upper corner
#' @param ucorner upper corner, a point whose coordinates must be 
#'   greater than those of the lower corner
#'
#' @return A \strong{rgl} mesh, i.e. a \code{mesh3d} object.
#' @export
#' @importFrom rgl cube3d translate3d scale3d
meshIsoCuboid <- function(lcorner, ucorner) {
    stopifnot(all(lcorner <= ucorner))
    center <- (lcorner + ucorner) / 2
    ax <- ucorner[1L] - lcorner[1L]
    ay <- ucorner[2L] - lcorner[2L]
    az <- ucorner[3L] - lcorner[3L]
    translate3d(
        scale3d(cube3d(), ax/2, ay/2, az/2),
        center[1L], center[2L], center[3L]
    )
}
