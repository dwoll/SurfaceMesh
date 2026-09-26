## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

library(rgl)

## ----------------------------------------------------------------------- //
## Hopf Torus
## ----------------------------------------------------------------------- //

meshHopfTorusHelper <- function(u, cos_v, sin_v, nlobes, A, alpha) {
  B <- pi/2 - (pi/2 - A)*cos(u*nlobes)
  C <- u + A*sin(2*u*nlobes)
  y1 <- 1 + cos(B)
  y23 <- sin(B) * c(cos(C), sin(C))
  y2 <- y23[1L]
  y3 <- y23[2L]
  x1 <- cos_v*y3 + sin_v*y2
  x2 <- cos_v*y2 - sin_v*y3
  x3 <- sin_v*y1
  x4 <- cos_v*y1
  yden <- sqrt(2*y1)
  if(is.null(alpha)) {
    t(cbind(x1, x2, x3) / (yden-x4))
  } else {
    t(acos(x4/yden)/(yden^alpha-abs(x4)^alpha)^(1/alpha) * cbind(x1, x2, x3))
  }
}

#' @title Hopf torus mesh
#' @description Triangle mesh of a Hopf torus.
#'
#' @param nlobes number of lobes of the Hopf torus, a positive integr
#' @param A parameter of the Hopf torus, number strictly between
#'   \code{0} and \code{pi/2}
#' @param alpha if not \code{NULL}, this is the exponent of a modified
#'   stereographic projection, a positive number; otherwise the ordinary
#'   stereographic projection is used
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- meshHopfTorus(nu = 90, nv = 90)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "forestgreen")
#' wire3d(mesh)
#' mesh <- mesnHopfTorus(nu = 90, nv = 90, alpha = 1.5)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "yellowgreen")
#' wire3d(mesh)
#'
#' @export
#' @importFrom rgl tmesh3d addNormals
meshHopfTorus <- function(nlobes = 3, A = 0.44, alpha = NULL, nu, nv) {
  # stopifnot(isStrictPositiveInteger(nlobes))
  # stopifnot(isPositiveNumber(A))
  # stopifnot(A < pi/2)
  # stopifnot(is.null(alpha) || isPositiveNumber(alpha))
  # stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vs    <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  u_ <- seq(0, 2*pi, length.out = nu + 1L)[-1L]
  v_ <- seq(0, 2*pi, length.out = nv + 1L)[-1L]
  cos_v <- cos(v_)
  sin_v <- sin(v_)
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  for(i in 1L:(nu-1L)) {
    i_nv <- i*nv
    vs[, (i_nv - nv + 1L):i_nv] <-
      meshHopfTorusHelper(u_[i], cos_v, sin_v, nlobes, A, alpha)
    k1 <- i_nv - nv
    k_ <- k1 + j_
    l_ <- k1 + jp1_
    m_ <- i_nv + j_
    tris1[, k_] <- rbind(k_, l_, m_)
    tris2[, k_] <- rbind(l_, i_nv + jp1_, m_)
  }
  i_nv <- nu*nv
  vs[, (i_nv - nv + 1L):i_nv] <-
    meshHopfTorusHelper(2*pi, cos_v, sin_v, nlobes, A, alpha)
  k1 <- i_nv - nv
  k_ <- k1 + j_
  l_ <- k1 + jp1_
  tris1[, k_] <- rbind(k_, l_, j_)
  tris2[, k_] <- rbind(l_, jp1_, j_)
  tmesh3d(vertices    = vs,
          indices     = cbind(tris1, tris2),
          normals     = NULL,
          homogeneous = FALSE)
}

# dataHopfTorus <- meshHopfTorus(nu=200L, nv=150L, nlobes=3L, A=0.44)
# wire3d(dataHopfTorus)
