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
## Solid Mobius Strip
## ----------------------------------------------------------------------- //

# a = 0.4, b = 0.1
meshSolidMobiusStrip <- function(x, y, z, a, b, nx, ny, nz) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    # : f(x, y, z)=0
    f <- function(x, y, z, a, b) {
        ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2 -
            4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2
    }
    x <- seq(-1.4, 1.4, length.out=nx)
    y <- seq(-1.7, 1.7, length.out=ny)
    z <- seq(-0.7, 0.7, length.out=nz)
    G <- expand.grid(x=x, y=y, z=z)
    voxel <- array(with(G, f(x, y, z, a, b)), dim=c(nx, ny, nz))
    surf  <- misc3d::computeContour3d(
        voxel, maxvol=max(voxel), level=0, x=x, y=y, z=z)
    
    mesh <- misc3d:::t2ve(misc3d::makeTriangles(surf))
    tmesh3d(vertices=mesh[["vb"]], indices=mesh[["ib"]])
}


# dataSolidMobiusStrip <- meshSolidMobiusStrip(a=0.4, b=0.1, nx=100L, ny=100L, nz=100L)
# wire3d(toRGL(m))
