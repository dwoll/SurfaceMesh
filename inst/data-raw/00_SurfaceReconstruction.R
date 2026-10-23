## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
## Surface Reconstruction Demo
## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //

library(SurfaceMesh)
library(rgl)

## unfinished

## ----------------------------------------------------------------------- //
## Dupin Cyclide
## ----------------------------------------------------------------------- //

mCyclide   <- dataCyclide
ptsCyclide <- t(mCyclide$vb[-4L, ])

cyclide_afs  <- reconstructAFS(ptsCyclide) # jetSmoothing=4L

cyclide_sss  <- reconstructSSS(ptsCyclide,
                               scaleIterations=2L,
                               neighbors=12L,
                               samples=250L,
                               separateShells=FALSE,
                               forceManifold=TRUE,
                               borderAngle=60)

cyclide_pois <- reconstructPoisson(ptsCyclide,
                                   # normalsFun=getNormalsFun(12L),
                                   spacing=0.01,
                                   smAngle=10,     # 20
                                   smRadius=5,    # 30
                                   smDistance=0.5)  # 0.375

wire3d(mCyclide)
shade3d(toRGL(cyclide_afs),  color="darkred")
shade3d(toRGL(cyclide_sss),  color="orangered")
shade3d(toRGL(cyclide_pois), color="orangered")

## ----------------------------------------------------------------------- //
## Hopf Torus
## ----------------------------------------------------------------------- //

mHopfTorus   <- dataHopfTorus
ptsHopfTorus <- mHopfTorus$vertices

ht_afs  <- reconstructAFS(ptsHopfTorus)
ht_pois <- reconstructPoisson(ptsHopfTorus,
                              spacing=0.01,
                              smAngle=10,
                              smDistance=0.4)

wire3d(toRGL(mHopfTorus))
shade3d(toRGL(ht_afs),  color="darkred")
shade3d(toRGL(ht_pois), color="orangered")

## ----------------------------------------------------------------------- //
## Solid Möbius strip
## ----------------------------------------------------------------------- //

mSolidMobiusStrip   <- meshSolidMobiusStrip(a=0.4, b=0.1, nx=100L, ny=100L, nz=100L)
ptsSolidMobiusStrip <- t(mSolidMobiusStrip$vb[-4L, ])

sms_afs  <- reconstructAFS(ptsSolidMobiusStrip)
sms_pois <- reconstructPoisson(ptsSolidMobiusStrip,
                               spacing=0.01)

wire3d(mSolidMobiusStrip)
shade3d(toRGL(sms_afs),  color="darkred")
shade3d(toRGL(sms_pois), color="orangered")

## ----------------------------------------------------------------------- //
## Toroidal Helix
## ----------------------------------------------------------------------- //

## ----------------------------------------------------------------------- //
## Tubular Knot
## ----------------------------------------------------------------------- //

mTubularKnot   <- dataTubularKnot
ptsTubularKnot <- t(mTubularKnot$vb[-4L, ])

tk_afs  <- reconstructAFS(mTubularKnot)
tk_pois <- reconstructPoisson(mTubularKnot,
                              spacing=0.01,
                              smDistance=0.9)

wire3d(mTubularKnot)
shade3d(toRGL(tk_afs),  color="violetred")
shade3d(toRGL(tk_pois), color="darkviolet")

## ----------------------------------------------------------------------- //
## Bunny
## ----------------------------------------------------------------------- //

data(bunny, package="onion")
bunny_afs  <- reconstructAFS(bunny)
bunny_sss  <- reconstructSSS(bunny, scaleIterations=2L, neighbors=6L)
bunny_pois <- reconstructPoisson(bunny,
                                 spacing=0.0001,
                                 smDistance=0.9)

shade3d(toRGL(bunny_afs),  color="violetred")
shade3d(toRGL(bunny_pois), color="darkviolet")

## ----------------------------------------------------------------------- //
## Stanford Dragon
## ----------------------------------------------------------------------- //

StanfordDragon_path <- "c:/users/dwollsch/Documents/SurfaceReconstruction/data-raw/StanfordDragon.txt"
ptsStanfordDragon   <- data.matrix(read.table(StanfordDragon_path))
sd_afs  <- reconstructAFS(ptsStanfordDragon)
sd_pois <- reconstructPoisson(ptsStanfordDragon,
                              spacing=0.0003)

shade3d(toRGL(sd_afs),  color="darkolivegreen4")
shade3d(toRGL(sd_pois), color="forestgreen")
