library(rgl)

## ----------------------------------------------------------------------- //
## Entzensberger Star
## ----------------------------------------------------------------------- //

## constant = 400L
meshEntzenStar <- function(n = 60L, constant=400L) {
    if (length(n) != 1L || !is.numeric(n) || is.na(n) ||
        !is.finite(n) || n < 3 || n != floor(n)) {
        stop("`n` must be a single integer of at least 3.", call. = FALSE)
    }
    n <- as.integer(n)
    
    m <- 2L * n
    theta <- 2 * pi * (0:(m - 1L)) / m
    phi <- rep(pi * seq_len(n - 1L) / n, each = m)
    theta <- rep(theta, times = n - 1L)
    
    ux <- sin(phi) * cos(theta)
    uy <- sin(phi) * sin(theta)
    uz <- cos(phi)
    q <- ux^2 * uy^2 + uy^2 * uz^2 + ux^2 * uz^2
    
    # Solve for t = r^2 in constant * t^2 * q = (1 - t)^3.
    lo <- numeric(length(q))
    hi <- rep(1, length(q))
    for (i in seq_len(60L)) {
        mid <- (lo + hi) / 2
        above <- constant * mid^2 * q > (1 - mid)^3
        hi[above] <- mid[above]
        lo[!above] <- mid[!above]
    }
    r <- sqrt((lo + hi) / 2)
    
    vertices <- cbind(
        c(0, r * ux, 0),
        c(0, r * uy, 0),
        c(1, r * uz, -1),
        1
    )
    
    ring_index <- function(j, k) 2L + (j - 1L) * m + (k %% m)
    k <- 0:(m - 1L)
    triangles <- list(
        rbind(rep(1L, m), ring_index(1L, k), ring_index(1L, k + 1L))
    )
    
    for (j in seq_len(n - 2L)) {
        a <- ring_index(j, k)
        a_next <- ring_index(j, k + 1L)
        b <- ring_index(j + 1L, k)
        b_next <- ring_index(j + 1L, k + 1L)
        triangles[[length(triangles) + 1L]] <-
            rbind(c(a, a_next), c(b, b), c(a_next, b_next))
    }
    
    bottom <- 2L + (n - 1L) * m
    triangles[[length(triangles) + 1L]] <-
        rbind(ring_index(n - 1L, k), rep(bottom, m),
              ring_index(n - 1L, k + 1L))
    
    rgl::tmesh3d(
        vertices = t(vertices),
        indices = do.call(cbind, triangles),
        homogeneous = TRUE
    )
}

# m <- meshEntzenStar(n=60, constant=100)
# shade3d(m, col="hotpink3")
# wire3d(m)