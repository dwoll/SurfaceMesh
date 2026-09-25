## ----------------------------------------------------------------------- //
## Code adapted from package
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## Spider Cage
## ----------------------------------------------------------------------- //

meshSpiderCage <- function(a, nx, ny, nz) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    f <- function(x, y, z, a) {
        (sqrt((x^2 - y^2)^2/(x^2 + y^2) + 3*(z*sin(a))^2) - 3)^2 +
            6*(sqrt((x*y)^2/(x^2 + y^2) + (z*cos(a))^2) - 1.5)^2
    }
    x <- seq(-4.7, 4.7, length.out=nx)
    y <- seq(-4.7, 4.7, length.out=ny)
    z <- seq(-2.8, 2.8, length.out=nz)
    G <- expand.grid(x=x, y=y, z=z)
    voxel <- array(with(G, f(x, y, z, a)), dim=c(nx, ny, nz))
    surf  <- misc3d::computeContour3d(voxel, max(voxel), 0.5, x=x, y=y, z=z)
    tmesh3d(vertices=t(surf),
            indices =matrix(seq_len(nrow(surf)), nrow=3L))
}

# dataSpiderCage <- meshSpiderCage(a=0.9, nx=250L, ny=250L, nz=250L)
# wire3d(dataSpiderCage)
