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

library(rgl)
library(SurfaceMesh)

tt     <- seq(0, 2, length.out=1000L)[-1L]
crcl1  <- cbind(cospi(tt)+0.5, sinpi(tt), 0)
crcl2  <- cbind(cospi(tt)-0.5, 0,         sinpi(tt))
points <- rbind(crcl1, crcl2)

# dataOloid <- getConvexHull(points, normals=TRUE)
dataOloid <- alphaWrap(points,alphaRel = 0.1, offsetRel = 1000)

mesh_rgl  <- toRGL(dataOloid)
wire3d(mesh_rgl)
