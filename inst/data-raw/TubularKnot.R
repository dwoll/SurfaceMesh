## ----------------------------------------------------------------------- //
## Code adapted from package
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## Tubular Knot
## ----------------------------------------------------------------------- //

meshTubularKnot <- function(p, q, a, nu, nv) {
    # cross product
    xprod <- function(v, w) {
        c(v[2] * w[3] - v[3] * w[2], 
          v[3] * w[1] - v[1] * w[3], 
          v[1] * w[2] - v[2] * w[1])
    }
    
    # knot curve
    knot <- function(t, p, q) {
        r <- cos(q*t)+2
        c(r * cos(p*t),
          r * sin(p*t),
          -sin(q*t))
    }
    
    # derivative (tangent)
    dknot <- function(t, p, q) {
        v <- c(-q*sin(q*t)*cos(p*t) - p*sin(p*t)*(cos(q*t)+2),
               -q*sin(q*t)*sin(p*t) + p*cos(p*t)*(cos(q*t)+2),
               -q*cos(q*t))
        v / sqrt(c(crossprod(v)))
    }
    
    # second derivative (normal)
    ddknot <- function(t, p, q) {
        v <- c(-q*(q*cos(q*t)*cos(p*t)-p*sin(p*t)*sin(q*t)) - 
                   p*(p*cos(p*t)*cos(q*t+2)-q*sin(q*t)*sin(p*t)),
               -q*(q*cos(q*t)*sin(p*t)+p*cos(p*t)*sin(q*t)) + 
                   p*(-p*sin(p*t)*cos(q*t+2)-q*sin(q*t)*cos(p*t)),
               q*q*sin(q*t))
        v / sqrt(c(crossprod(v)))
    }
    
    # binormal
    bnrml <- function(t, p, q) {
        v <- xprod(dknot(t, p, q), ddknot(t, p, q))
        v / sqrt(c(crossprod(v)))
    }
    nu       <- as.integer(nu)
    nv       <- as.integer(nv)
    vertices <- matrix(NA_real_, nrow=3L, ncol=nu*nv)
    u_ <- seq(0, 2*pi, length.out = nu+1L)[-1L]
    v_ <- seq(0, 2*pi, length.out = nv+1L)[-1L]
    for(i in seq_len(nu)) {
      u <- u_[i]
      for(j in seq_len(nv)) {
        v <- v_[j]
        h <- knot(u, p, q)
        vertices[, (i-1)*nv+j] <- 
          h + 
          a*(cos(v) * 
               (cos(u)*ddknot(u, p, q) + 
                  sin(u)*bnrml(u, p, q)) + 
               sin(v) * 
               (-sin(u)*ddknot(u, p, q) + 
                  cos(u)*bnrml(u, p, q)))
      }
    }
    tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
    tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
    for(i in seq_len(nu)) {
      ip1 <- ifelse(i == nu, 1L, i+1L)
      for(j in seq_len(nv)) {
        jp1 <- ifelse(j==nv, 1L, j+1L)
        tris1[,(i-1)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)   
        tris2[,(i-1)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)   
      }
    }
    tmesh3d(vertices=vertices, indices=cbind(tris1, tris2), homogeneous=FALSE)
}

# dataTubularKnot <- meshTubularKnot(p=3, q=11, a=0.55, 10*60, 6)
# wire3d(dataTubularKnot)
