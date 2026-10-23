# SurfaceMesh 0.1.1 (2026-10-19)

 * Add `isOutwardOriented()` to check whether mesh has positive orientation.
 * Add `reverseFaceOrientations()` to reverse the orientation of (selected) faces.
 * Add `remeshACVD()` for approximated centroidal Voronoi diagram remeshing.
 * `reconstructAFS()`: change argument name `jetSmoothing` to `nNeighbors` for consistency.
 * `reconstructAFS()`: change argument names `neighbors`, `scaleIterations`, `samples` to `nNeighbors`, `nScaleIter`, `nSamples` for consistency.
 * `remeshSmoothShape()`: change argument name `indices` to `faces` for clarity.

# SurfaceMesh 0.1.0 (2026-10-18)

Initial release to CRAN.
