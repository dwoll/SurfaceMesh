## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //
## Surface Reconstruction Demo
## ----------------------------------------------------------------------- //
## ----------------------------------------------------------------------- //

library(SurfaceMesh)

## unfinished

## ----------------------------------------------------------------------- //
## Cyclide
## ----------------------------------------------------------------------- //

mCyclide   <- meshCyclide(aa=0.94, cc=0.34, dd=0.56)
ptsCyclide <- t(mCyclide$vb[-4L, ])

cyclide_afs  <- reconstructAFS(ptsCyclide)
cyclide_pois <- reconstructPoisson(ptsCyclide,
                                   spacing=0.01)

wire3d(mCyclide)
shade3d(toRGL(cyclide_afs),  color="darkred")
shade3d(toRGL(cyclide_pois), color="orangered")

## ----------------------------------------------------------------------- //
## Hopf Torus
## ----------------------------------------------------------------------- //

mHopfTorus   <- meshHopfTorus(nu=200L, nv=150L, nlobes=3L, A=0.44)
ptsHopfTorus <- t(mHopfTorus$vb[-4L, ])

ht_afs  <- reconstructAFS(ptsHopfTorus)
ht_pois <- reconstructPoisson(ptsHopfTorus,
                              spacing=0.01,
                              smAngle=10,
                              smDistance=0.4)

wire3d(mHopfTorus)
shade3d(toRGL(ht_afs),  color="darkred")
shade3d(toRGL(ht_pois), color="orangered")

## ----------------------------------------------------------------------- //
## ICN5D eight
## ----------------------------------------------------------------------- //

mICN5deight <- meshICN5Deight(a=2.4, nx=200L, ny=200L, nz=200L)
wire3d(mICN5deight)

## ----------------------------------------------------------------------- //
## Oloid
## ----------------------------------------------------------------------- //

## ----------------------------------------------------------------------- //
## OrthoCircle
## ----------------------------------------------------------------------- //

mOrthoCircle   <- meshOrthoCircle(a=0.075, b=3, nx=100L, ny=100L, nz=100L)
ptsOrthoCircle <- t(mOrthoCircle$vb[-4L, ])
wire3d(mOrthoCircle)

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
## Spider Cage
## ----------------------------------------------------------------------- //

mSpiderCage   <- meshSpiderCage(a=0.9, nx=250L, ny=250L, nz=250L)
ptsSpiderCage <- t(mSpiderCage$vb[-4L, ])

dsc_afs  <- reconstructAFS(ptsSpiderCage)
dsc_pois <- reconstructPoisson(ptsSpiderCage,
                               spacing=0.01,
                               smDistance=0.9)

wire3d(mSpiderCage)
shade3d(toRGL(dsc_afs),  color="violetred")
shade3d(toRGL(dsc_pois), color="darkviolet")

## ----------------------------------------------------------------------- //
## Tubular Knot
## ----------------------------------------------------------------------- //

mTubularKnot   <- meshTubularKnot(p=3, q=11, a=0.55, 10*60, 6)
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

## ----------------------------------------------------------------------- //
## Sphere
## ----------------------------------------------------------------------- //

Sphere_path <- "c:/apps/CGAL-6.2_/data/points_3/sphere_20k.off"
ptsSphere0  <- readMeshFile(Sphere_path)$vertices
ptsSphere   <- removeOutliers(ptsSphere0)
sphere_pois <- reconstructPoisson(ptsSphere,
                                  normalsFun=getNormalsFun(12L, method="PCA"),
                                  # spacing=0.5,
                                  smAngle=10,
                                  smRadius=10,
                                  smDistance=0.5)

wire3d(toRGL(sphere_pois))
