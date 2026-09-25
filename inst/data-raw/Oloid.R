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
library(SurfaceMesh)

## ----------------------------------------------------------------------- //
## Oloid
## ----------------------------------------------------------------------- //

meshOloid <- function(len) {
    tt     <- seq(0, 2, length.out=len)[-1L]
    crcl1  <- cbind(cospi(tt)+0.5, sinpi(tt), 0)
    crcl2  <- cbind(cospi(tt)-0.5, 0,         sinpi(tt))
    points <- rbind(crcl1, crcl2)
    getConvexHull(points)
}

# mOloid <- meshOloid(len=500L)
# wire3d(toRGL(mOloid))
