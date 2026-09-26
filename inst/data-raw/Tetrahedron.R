## ----------------------------------------------------------------------- //
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## Quartic surface
## this does not look right: https://math.stackexchange.com/a/2999723
## x^4 + 2x^2y^2 + 2x^2z^2 + y^4 + 2y^2z^2 + z^4 + 8xyz + 10x^2 - 10y^2 - 10z^2 + 20 = 0
## ----------------------------------------------------------------------- //

#' @title A mesh of a quartic surface
#' @description Triangle mesh of the implicit quartic surface defined by
#'   \eqn{(x^2+y^2+z^2)^2 + 8xyz + 10x^2 - 10y^2 - 10z^2 + 20 = 0}.
#' @param nx,ny,nz numbers of grid subdivisions along the x, y, z axes used
#'   to sample the implicit function, integers (at least 3)
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @author Daniel Wollschlaeger.
#'
#' @examples
#' library(rgl)
#' mesh <- meshTetrahedron(nx=100L, ny=100L, nz=100L)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' shade3d(mesh, color="steelblue")
#' wire3d(mesh)
#'
#' @export
#' @importFrom rgl tmesh3d
meshTetrahedron <- function(nx=100L, ny=100L, nz=100L) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    f <- function(x, y, z) {
        x2 <- x*x
        y2 <- y*y
        z2 <- z*z
        (x2 + y2 + z2)^2 + 8*x*y*z + 10*x2 - 10*y2 - 10*z2 + 20
    }
    # the surface is contained in roughly x in [-1.2, 1.2], y,z in [-2.7, 2.7]
    x <- seq(-1.5, 1.5, length.out=nx)
    y <- seq(-3,   3,   length.out=ny)
    z <- seq(-3,   3,   length.out=nz)
    G <- expand.grid(x=x, y=y, z=z)
    voxel <- array(with(G, f(x, y, z)), dim=c(nx, ny, nz))
    surf  <- misc3d::computeContour3d(voxel, level=0, x=x, y=y, z=z)
    tmesh3d(vertices=t(surf),
            indices =matrix(seq_len(nrow(surf)), nrow=3L))
}

# dataTetrahedron <- meshTetrahedron(nx=100L, ny=100L, nz=100L)
# wire3d(dataTetrahedron)
