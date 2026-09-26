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

path_in <- system.file("data-raw", package="SurfaceMesh")
source(paste0(path_in, "/Torus.R"))

## ----------------------------------------------------------------------- //
## Dupin Cyclide
## ----------------------------------------------------------------------- //

#' @title Cyclide mesh
#' @description Triangle mesh of a Dupin cyclide.
#' @param a,c,mu cyclide parameters, positive numbers such that
#'   \code{c < mu < a}
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- meshCyclide(a = 97, c = 32, mu = 57)
#' sphere <- meshSphere(x = 32, y = 0, z = 0, r = 40)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "chartreuse")
#' wire3d(mesh)
#' shade3d(sphere, color = "red")
#' wire3d(sphere)
#'
#' @export
#' @importFrom rgl tmesh3d
meshCyclide <- function(a, c, mu, nu = 90L, nv = 40L){
  stopifnot(c > 0, a > mu, mu > c)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  normals  <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  b2 <- a * a - c * c
  bb <- sqrt(b2 * (mu * mu - c * c))
  omega <- (a * mu + bb) / c
  Omega0 <- c(omega, 0, 0)
  inversion <- function(M) {
    OmegaM <- M - Omega0
    k <- c(crossprod(OmegaM))
    OmegaM / k + Omega0
  }
  h <- (c * c) / ((a - c) * (mu - c) + bb)
  r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
  R <- (h * (a - c)) / ((a - c) * (mu + c) + bb)
  # a' = a/c; mu' = mu/c
  # R/r = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1) + sqrt((a'^2-1)*(mu'^2-1))) / ((a'-1)*(mu'+1)+sqrt((a'^2-1)*(mu'^2-1)))
  # = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1)/(a'-1) + sqrt((mu'^2-1)*(a'+1)/(a'-1))) / ((mu'+1) + sqrt((mu'^2-1)*(a'+1)/(a'-1)))
  # = ((aa+1) * (1 + sqrt((aa-1)/(aa+1))*sqrt((muu+1)/(muu-1))) / ((muu+1) * (1  + sqrt((aa+1)/(aa-1))*sqrt((muu-1)/(muu+1)))))
  #(Wolfram) muu = sqrt(1 + (aa^2-1)/ratio^2)
  bb2 <- b2 * (mu * mu - c * c)
  denb1 <- c * (a*c - mu*c + c*c - a*mu - bb)
  b1 <- (a*mu*(c-mu)*(a+c) - bb2 + c*c + bb*(c*(a-mu+c) - 2*a*mu))/denb1
  denb2 <- c * (a*c - mu*c - c*c + a*mu + bb)
  b2 <- (a*mu*(c+mu)*(a-c) + bb2 - c*c + bb*(c*(a-mu-c) + 2*a*mu))/denb2
  omegaT <- (b1 + b2)/2
  OmegaT <- c(omegaT, 0, 0)
  tormesh <- meshTorus(R, r, nu = nu, nv = nv)
  rtnormals <- r * tormesh[["normals"]][1L:3L, ]
  xvertices <- tormesh[["vb"]][1L:3L, ] + OmegaT
  for(i in 1L:nu){
    k0 <- i * nv - nv
    for(j in 1L:nv){
      k <- k0 + j
      rtnormal <- rtnormals[, k]
      xvertex <- xvertices[, k]
      vertex <- inversion(xvertex)
      vertices[, k] <- vertex
      foo <- vertex - inversion(rtnormal + xvertex)
      normals[k, ] <- foo / sqrt(c(crossprod(foo)))
    }
  }
  tmesh3d(
    vertices    = vertices,
    indices     = tormesh[["it"]],
    normals     = normals,
    homogeneous = FALSE
  )
}

#' @title Conformal cyclide mesh
#' @description Cyclide mesh of a Dupin cyclide constructed from a
#'   conformal parameterization.
#'
#' @param a,c cyclide parameters, \code{a > c > 0}; see the figure in the
#'   documentation of \code{\link{cyclideMesh}}
#' @param aspectRatio the aspect ratio, positive number
#' @param normals Boolean, whether to add normals to the mesh
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @export
#' @importFrom rgl tmesh3d addNormals
meshConformalCyclide <- function(
    a, c, aspectRatio, normals = TRUE, nu = 90L, nv = 40L) {
  stopifnot(c > 0, a > c, aspectRatio > 0)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  b2 <- a*a - c*c
  ar2plus1 <- aspectRatio*aspectRatio + 1
  ratio <- sqrt(ar2plus1)
  mu <- sqrt(c*c + b2/ar2plus1)
  msg <- sprintf("The value of `mu` is %s.", formatC(mu))
  bb <- b2/ratio
  bb2 <- b2 * b2 / ar2plus1
  omega <- (a * mu + bb) / c
  Omega0 <- c(omega, 0, 0)
  inversion <- function(M) {
    OmegaM <- M - Omega0
    k <- c(crossprod(OmegaM))
    OmegaM / k + Omega0
  }
  h <- (c * c) / ((a - c) * (mu - c) + bb)
  r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
  R <- ratio * r
  # a' = a/c; mu' = mu/c
  # R/r = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1) + sqrt((a'^2-1)*(mu'^2-1))) / ((a'-1)*(mu'+1)+sqrt((a'^2-1)*(mu'^2-1)))
  # = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1)/(a'-1) + sqrt((mu'^2-1)*(a'+1)/(a'-1))) / ((mu'+1) + sqrt((mu'^2-1)*(a'+1)/(a'-1)))
  # = ((aa+1) * (1 + sqrt((aa-1)/(aa+1))*sqrt((muu+1)/(muu-1))) / ((muu+1) * (1  + sqrt((aa+1)/(aa-1))*sqrt((muu-1)/(muu+1)))))
  #(Wolfram) muu = sqrt(1 + (aa^2-1)/ratio^2)
  denb1 <- c * (a*c - mu*c + c*c - a*mu - bb)
  b1 <- (a*mu*(c-mu)*(a+c) - bb2 + c*c + bb*(c*(a-mu+c) - 2*a*mu))/denb1
  denb2 <- c * (a*c - mu*c - c*c + a*mu + bb)
  b2 <- (a*mu*(c+mu)*(a-c) + bb2 - c*c + bb*(c*(a-mu-c) + 2*a*mu))/denb2
  omegaT <- (b1 + b2)/2
  OmegaT <- c(omegaT, 0, 0)
  tormesh <- meshTorus2(R, r, nu = nu, nv = nv)
  xvertices <- tormesh[["vb"]][1L:3L, ] + OmegaT
  for(i in 1L:nu){
    k0 <- i * nv - nv
    for(j in 1L:nv){
      k <- k0 + j
      vertices[, k] <- inversion(xvertices[, k])
    }
  }
  mesh <- tmesh3d(
    vertices    = vertices,
    indices     = tormesh[["it"]],
    homogeneous = FALSE
  )
  if(normals) {
    mesh <- addNormals(mesh)
  }
  mesh
}

## ----------------------------------------------------------------------- //
## Dupin Cyclide
## ----------------------------------------------------------------------- //

# aa = 0.94, cc = 0.34, dd = 0.56
# meshCyclide <- function(aa, cc, dd) {
#     if(!requireNamespace("misc3d", quietly=TRUE)) {
#         stop("Package `misc3d` required but not found.")
#     }
#     bb = sqrt(aa^2-cc^2)
#     fx <- function(u, v) {
#         (dd*(cc-aa*cos(u)*cos(v)) + bb^2*cos(u)) / (aa-cc*cos(u)*cos(v))
#     }
#
#     fy <- function(u, v) {
#         (bb*sin(u)*(aa-dd*cos(v))) / (aa-cc*cos(u)*cos(v))
#     }
#
#     fz <- function(u, v) {
#         (bb*sin(v)*(cc*cos(u)-dd)) / (aa-cc*cos(u)*cos(v))
#     }
#
#     tris <- misc3d::parametric3d(fx, fy, fz,
#                                  umin=0, umax=2*pi, vmin=0, vmax=2*pi,
#                                  n=50L, engine="none")
#
#     mesh <- misc3d:::t2ve(tris)
#     tmesh3d(vertices=mesh[["vb"]], indices=mesh[["ib"]])
# }

# dataCyclide <- meshCyclide(aa=0.94, cc=0.34, dd=0.56)
# wire3d(dataCyclide)
