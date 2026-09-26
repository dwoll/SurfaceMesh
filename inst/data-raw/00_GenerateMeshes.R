## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
## Generate sample meshes
## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //

library(rgl)
library(SurfaceMesh)

path_in  <- system.file("data-raw", package="SurfaceMesh")
path_out <- "c:/Users/dwollsch/Documents/SurfaceMesh/data/"

source(paste0(path_in, "/Cyclide.R"))
source(paste0(path_in, "/HopfTorus.R"))
source(paste0(path_in, "/ICN5Deight.R"))
source(paste0(path_in, "/IsoCuboid.R"))
source(paste0(path_in, "/Oloid.R"))
source(paste0(path_in, "/OrthoCircle.R"))
source(paste0(path_in, "/PentaPrism.R"))
source(paste0(path_in, "/Septuaginta.R"))
source(paste0(path_in, "/SolidMobiusStrip.R"))
source(paste0(path_in, "/Sphere.R"))
source(paste0(path_in, "/SpiderCage.R"))
source(paste0(path_in, "/ToroidalHelix.R"))
source(paste0(path_in, "/Torus.R"))
source(paste0(path_in, "/TruncIcosahedron.R"))
source(paste0(path_in, "/TubularKnot.R"))

## ----------------------------------------------------------------------- //
## Cyclide
## ----------------------------------------------------------------------- //

dataCyclide <- meshCyclide(a = 97, c = 32, mu = 57)
view3d(0, -35, zoom=0.8)
shade3d(dataCyclide, color="orangered")

save(dataCyclide, file=paste0(path_out, "dataCyclide.rda"))

## ----------------------------------------------------------------------- //
## Hopf Torus
## ----------------------------------------------------------------------- //

dataHopfTorus <- meshHopfTorus(nu=200L, nv=150L, nlobes=3L, A=0.44)
view3d(35, -20, zoom=0.8)
shade3d(dataHopfTorus, color="darkred")

# dataHopfTorus <- remeshIsotropic(makeMeshValid(dataHopfTorus),
#                                  method="adaptive",
#                                  edgeMin=0.7,
#                                  edgeMax=5)
wire3d(toRGL(dataHopfTorus))
save(dataHopfTorus, file=paste0(path_out, "dataHopfTorus.rda"))

## ----------------------------------------------------------------------- //
## ICN5D eight
## ----------------------------------------------------------------------- //

dataICN5Deight <- try(meshICN5Deight(a=2.4, nx=200L, ny=200L, nz=200L))
if(!inherits(dataICN5Deight, "try-error")) {
    view3d(30, -30, zoom=0.8)
    shade3d(dataICN5Deight, color="violetred")

    # dataICN5Deight <- try(meshICN5Deight(a=2.4, nx=200L, ny=200L, nz=200L))
    # dataICN5Deight <- Rvcg::vcgUniformRemesh(dataICN5Deight,
    #                                          multiSample=TRUE,
    #                                          mergeClost=TRUE)
    # wire3d(toRGL(dataICN5Deight))
    #
    save(dataICN5Deight, file=paste0(path_out, "dataICN5Deight.rda"))
}

## ----------------------------------------------------------------------- //
## Oloid
## ----------------------------------------------------------------------- //

dataOloid <- meshOloid(len=500L)
view3d(20, 15, zoom=0.8)
shade3d(toRGL(dataOloid), color="darkviolet")

# dataOloid <- remeshIsotropic(dataOloid,
#                              method="adaptive",
#                              edgeMin=0.2,
#                              edgeMax=2,
#                              dihedralAngle=60)

save(dataOloid, file=paste0(path_out, "dataOloid.rda"))

## ----------------------------------------------------------------------- //
## OrthoCircle
## ----------------------------------------------------------------------- //

dataOrthoCircle <- try(meshOrthoCircle(a=0.075, b=3, nx=100L, ny=100L, nz=100L))
if(!inherits(dataOrthoCircle, "try-error")) {
    view3d(-15, -10, zoom=0.8)
    shade3d(dataOrthoCircle, color="darkolivegreen4")

    # dataOrthoCircle <- try(meshOrthoCircle(a=0.075, b=3, nx=100L, ny=100L, nz=100L))
    # dataOrthoCircle <- Rvcg::vcgUniformRemesh(dataOrthoCircle,
    #                                          multiSample=TRUE,
    #                                          mergeClost=TRUE)
    # wire3d(toRGL(dataOrthoCircle))

    save(dataOrthoCircle, file=paste0(path_out, "dataOrthoCircle.rda"))
}

## ----------------------------------------------------------------------- //
## Pentagrammic Prism
## ----------------------------------------------------------------------- //

dataPentaPrism <- meshPentaPrism(raw=TRUE,  triangulate=TRUE)
meshPentaPrism <- meshPentaPrism(raw=FALSE, triangulate=TRUE)

view3d(-25, -25, zoom=0.8)
shade3d(toRGL(meshPentaPrism), color="forestgreen")

save(dataPentaPrism, file=paste0(path_out, "dataPentaPrism.rda"))

## ----------------------------------------------------------------------- //
## Septuaginta
## ----------------------------------------------------------------------- //

dataSeptuaginta <- meshSeptuaginta()
view3d(-10, -20, zoom=0.8)
shade3d(dataSeptuaginta, color="orangered")

save(dataSeptuaginta, file=paste0(path_out, "dataSeptuaginta.rda"))

## ----------------------------------------------------------------------- //
## Solid Möbius strip
## ----------------------------------------------------------------------- //

dataSolidMobiusStrip <- try(meshSolidMobiusStrip(a=0.4, b=0.1, nx=100L, ny=100L, nz=100L))
if(!inherits(dataSolidMobiusStrip, "try-error")) {
    view3d(-10, -20, zoom=0.8)
    shade3d(dataSolidMobiusStrip, color="darkred")

    # dataMobiusStrip      <- try(meshSolidMobiusStrip(a=0.4, b=0.1, nx=100L, ny=100L, nz=100L))
    # dataSolidMobiusStrip <- Rvcg::vcgUniformRemesh(dataSolidMobiusStrip,
    #                          multiSample=FALSE,
    #                          mergeClost =TRUE)
    # wire3d(dataSolidMobiusStrip)

    save(dataSolidMobiusStrip, file=paste0(path_out, "dataSolidMobiusStrip.rda"))
}

## ----------------------------------------------------------------------- //
## Sphere
## ----------------------------------------------------------------------- //

dataSphere <- meshSphere(r=1, nIter=3)
view3d(-20, -20, zoom = 0.75)
shade3d(dataSphere, color="darkgoldenrod")
wire3d(dataSphere)

save(dataSphere, file=paste0(path_out, "dataSphere.rda"))

## ----------------------------------------------------------------------- //
## Spider Cage
## ----------------------------------------------------------------------- //

dataSpiderCage <- try(meshSpiderCage(a=0.9, nx=250L, ny=250L, nz=250L))
if(!inherits(dataSpiderCage, "try-error")) {
    view3d(-25, -45, zoom=0.8)
    shade3d(dataSpiderCage, color="darkviolet")

    # dataSpiderCage <- try(meshSpiderCage(a=0.9, nx=250L, ny=250L, nz=250L))
    # dataSpiderCage <- Rvcg::vcgUniformRemesh(dataSpiderCage,
    #                                           multiSample=TRUE,
    #                                           mergeClost=TRUE)
    # wire3d(toRGL(dataSpiderCage))

    save(dataSpiderCage, file=paste0(path_out, "dataSpiderCage.rda"))
}

## ----------------------------------------------------------------------- //
## Toroidal Helix
## ----------------------------------------------------------------------- //

dataToroHelix <- meshToroHelix(
    R=5, r=1.35, w=10, a=0.7, nu=10L*20L, nv=30L, alpha=1, twists=1)

# dataToroHelix <- Rvcg::vcgIsotropicRemeshing(mToroHelix,
#                                              TargetLen = 0.4,
#                                              FeatureAngleDeg = 40)

view3d(15, -20, zoom=0.8)
shade3d(dataToroHelix, color="violetred")

save(dataToroHelix, file=paste0(path_out, "dataToroHelix.rda"))

## ----------------------------------------------------------------------- //
## Torus
## ----------------------------------------------------------------------- //

dataTorus <- meshTorus(R=3, r=1)
view3d(-20, -20, zoom = 0.75)
shade3d(dataTorus, color = "green")
wire3d(dataTorus)

save(dataTorus, file=paste0(path_out, "dataTorus.rda"))

## ----------------------------------------------------------------------- //
## Truncated Icosahedron
## ----------------------------------------------------------------------- //

dataTruncIcosahedron <- meshTruncIcosahedron(raw=TRUE, triangulate=TRUE)
meshTruncIcosahedron <- meshTruncIcosahedron(raw=FALSE, triangulate=TRUE)
view3d(-15, -10, zoom=0.8)
shade3d(toRGL(meshTruncIcosahedron), color="darkolivegreen4")

save(dataTruncIcosahedron, file=paste0(path_out, "dataTruncIcosahedron.rda"))

## ----------------------------------------------------------------------- //
## Tubular Knot
## ----------------------------------------------------------------------- //

dataTubularKnot <- meshTubularKnot(p=3, q=11, a=0.55, 10*60, 6)
view3d(-15, -15, zoom=0.8)
shade3d(dataTubularKnot, color="forestgreen")

save(dataTubularKnot, file=paste0(path_out, "dataTubularKnot.rda"))
