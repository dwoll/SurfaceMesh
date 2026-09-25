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
## Toroidal Helix
## ----------------------------------------------------------------------- //

meshToroHelix <- function(R, r, w, a, nu, nv, alpha=1, twists=2) {
    # cross product
    xprod <- function(v, w) {
        c(v[2] * w[3] - v[3] * w[2], 
          v[3] * w[1] - v[1] * w[3], 
          v[1] * w[2] - v[2] * w[1])
    }
    
    # helix curve
    helix <- function(t, R, r, w) {
        c((R + r*cos(t)) * cos(t/w),
          (R + r*cos(t)) * sin(t/w),
          r*sin(t))
    }
    
    # derivative (tangent)
    dhelix <- function(t, R, r, w) {
        v <- c(
            -r*sin(t)*cos(t/w) - (R+r*cos(t))/w*sin(t/w),
            -r*sin(t)*sin(t/w) + (R+r*cos(t))/w*cos(t/w),
            r*cos(t))
        v / sqrt(c(crossprod(v)))
    }
    
    # second derivative (normal)
    ddhelix <- function(t, R, r, w) {
        v <- c(
            -r*cos(t)*cos(t/w) + r*sin(t)/w*sin(t/w) +
                r*sin(t)/w*sin(t/w) - (R+r*cos(t))/w^2*cos(t/w),
            -r*cos(t)*sin(t/w) - r*sin(t)/w*cos(t/w) -
                r*sin(t)/w*cos(t/w) - (R+r*cos(t))/w^2*sin(t/w),
            -r*sin(t))
        v / sqrt(c(crossprod(v)))
    }
    
    # binormal
    bnrml <- function(t, R, r, w) {
        v <- xprod(dhelix(t, R, r, w), ddhelix(t, R, r, w))
        v / sqrt(c(crossprod(v)))
    }
    
    # mesh maker
    scos     <- function(x,alpha) sign(cos(x)) * abs(cos(x))^alpha
    ssin     <- function(x,alpha) sign(sin(x)) * abs(sin(x))^alpha
    vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
    u_ <- seq(0, w*2*pi, length.out = nu+1)[-1L]
    v_ <- seq(0, 2*pi, length.out = nv+1)[-1L]
    for(i in seq_len(nu)) {
      u <- u_[i]
      for(j in seq_len(nv)) {
        v <- v_[j]
        h <- helix(u, R, r, w)
        vertices[,(i-1)*nv+j] <-
          h +
          a*(scos(v,alpha) *
               (cos(twists*u)*ddhelix(u, R, r, w) +
                  sin(twists*u)*bnrml(u, R, r, w)) +
               ssin(v,alpha) *
               (-sin(twists*u)*ddhelix(u, R, r, w) +
                  cos(twists*u)*bnrml(u, R, r, w)))
      }
    }
    tris1 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
    tris2 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
    nv <- as.integer(nv)
    for(i in seq_len(nu)) {
      ip1 <- ifelse(i == nu, 1L, i+1L)
      for(j in seq_len(nv)) {
        jp1 <- ifelse(j == nv, 1L, j+1L)
        tris1[,(i-1L)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)
        tris2[,(i-1L)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)
      }
    }
    tmesh3d(vertices=vertices, indices=cbind(tris1, tris2), homogeneous=FALSE)
}


# dataToroidalHelix <- meshToroidalHelix(
#   R=5, r=1.35, w=10, a=0.7, nu=10*20, nv=30, alpha=1, twists=1)
# wire3d(dataToroidalHelix)
