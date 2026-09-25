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
## Cyclide
## ----------------------------------------------------------------------- //

# aa = 0.94, cc = 0.34, dd = 0.56
meshCyclide <- function(aa, cc, dd) {
    if(!requireNamespace("misc3d", quietly=TRUE)) {
        stop("Package `misc3d` required but not found.")
    }
    bb = sqrt(aa^2-cc^2)
    fx <- function(u, v) {
        (dd*(cc-aa*cos(u)*cos(v)) + bb^2*cos(u)) / (aa-cc*cos(u)*cos(v))
    }
    
    fy <- function(u, v) {
        (bb*sin(u)*(aa-dd*cos(v))) / (aa-cc*cos(u)*cos(v))
    }
    
    fz <- function(u, v) {
        (bb*sin(v)*(cc*cos(u)-dd)) / (aa-cc*cos(u)*cos(v))
    }

    tris <- misc3d::parametric3d(fx, fy, fz,
                                 umin=0, umax=2*pi, vmin=0, vmax=2*pi,
                                 n=50L, engine="none")
    
    mesh <- misc3d:::t2ve(tris)
    tmesh3d(vertices=mesh[["vb"]], indices=mesh[["ib"]])
}

# dataCyclide <- meshCyclide(aa=0.94, cc=0.34, dd=0.56)
# wire3d(dataCyclide)

       