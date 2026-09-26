## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/SurfaceReconstruction/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## ICN5D eight
## from user ICN5D on forum http://hi.gher.space/forum/
## ----------------------------------------------------------------------- //

meshICN5Deight <- function(a, nx, ny, nz) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    # equation f(x, y, z, a) = 0.5
    f <- function(x, y, z, a) {
        (sqrt((x^2 - a^2)^2/(x^2 + a^2) + y^2) - 3)^2 +
            (sqrt((x*a)^2/(x^2 + a^2) + z^2) - 1.5)^2
    }
    x <- seq(-5, 5, length.out=nx)
    y <- seq(-4, 4, length.out=ny)
    z <- seq(-2.5, 2.5, length.out=nz)
    G <- expand.grid(x=x, y=y, z=z)
    voxel   <- array(with(G, f(x, y, z, a)), dim=c(nx, ny, nz))
    surface <- misc3d::computeContour3d(
      voxel, maxvol=max(voxel), level=0.5, x=x, y=y, z=z)
    
    mesh <- misc3d:::t2ve(misc3d::makeTriangles(surface))
    tmesh3d(vertices=mesh[["vb"]], indices=mesh[["ib"]])
}

# dataICN5deight <- meshICN5Deight(a=2.4, nx=200L, ny=200L, nz=200L)
# wire3d(dataICN5deight)
