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
## Ortho Circle
## ----------------------------------------------------------------------- //

# a=0.075, b=3, nx=100L, ny=100L, nz=100L
meshOrthoCircle <- function(a, b, nx, ny, nz) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    f <- function(x, y, z, a, b) {
        x2 <- x*x
        y2 <- y*y
        z2 <- z*z
        xy2 <- x2 + y2 - 1
        yz2 <- y2 + z2 - 1
        zx2 <- z2 + x2 - 1
        (xy2*xy2 + z2) * (yz2*yz2 + x2) * (zx2*zx2 + y2) - a*a*(1 + b*(x2 + y2 + z2))
    }
    x <- seq(-1.3, 1.3, length.out=nx)
    y <- seq(-1.3, 1.3, length.out=ny)
    z <- seq(-1.3, 1.3, length.out=nz)
    G <- expand.grid(x=x, y=y, z=z)
    voxel <- array(with(G, f(x, y, z, a, b)), c(nx, ny, nz))
    surf  <- misc3d::computeContour3d(voxel, level=0, x=x, y=y, z=z)
    tmesh3d(vertices=t(surf),
            indices =matrix(seq_len(nrow(surf)), nrow=3L))
}

# dataOrthoCircle <- meshOrthoCircle(a=0.075, b=3, nx=100L, ny=100L, nz=100L)
# wire3d(dataOrthoCircle)
