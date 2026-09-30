## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title A mesh of a Hopf torus
#' @description A \code{CGALmesh} object, i.e., the output of
#'   \code{\link[SurfaceMesh]{makeMesh}} representing a Hopf torus.
#' @format A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent.
#' @examples
#' ## The following files were used to generate this mesh.
#' system.file("data-raw", "HopfTorus.R",         package="SurfaceMesh")
#' system.file("data-raw", "00_GenerateMeshes.R", package="SurfaceMesh")
#'
#' library(SurfaceMesh)
#' library(rgl)
#' view3d(35, -20, zoom=0.8)
#' shade3d(toRGL(dataHopfTorus), color="darkred")
"dataHopfTorus"
