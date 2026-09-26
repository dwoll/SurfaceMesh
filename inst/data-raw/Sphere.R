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
## Torus
## ----------------------------------------------------------------------- //

#' @title Sphere mesh
#' @description Triangle mesh of a sphere.
#' @param x,y,z coordinates of the center
#' @param r radius
#' @param iterations number of iterations (the mesh is obtained by iteratively
#'   subdividing the faces of an icosahedron)
#' @return A \strong{rgl} mesh (class \code{mesh3d}).
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @export
#' @importFrom rgl subdivision3d icosahedron3d
meshSphere <- function(x = 0, y = 0, z = 0, r = 1, nIter = 3L) {
    # stopifnot(isNumber(x), isNumber(y), isNumber(z))
    # stopifnot(isPositiveNumber(r))
    # stopifnot(isStrictPositiveInteger(nIter))
    sphere <- subdivision3d(icosahedron3d(), depth = nIter)
    vs <- sphere[["vb"]][1L:3L, ]
    h <- sqrt(apply(vs, 2L, function(x) sum(x * x)))
    vs[1L, ] <- vs[1L, ] / h
    vs[2L, ] <- vs[2L, ] / h
    vs[3L, ] <- vs[3L, ] / h
    sphere[["vb"]] <- rbind(r*vs + c(x, y, z), 1)
    sphere[["normals"]] <- vs
    sphere
}