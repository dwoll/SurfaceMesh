# TODO

  * vignette
  * `remeshIsotropic()`
      * does not work well for pentagrammic prism - why? (`Rvcg::vcgIsotropicRemeshing()` works)

# Wishlist

  * smoothing
      * `PMP::angle_and_area_smoothing()`
      * `PMP::tangential_relaxation()`
  * remeshing
      * https://www.cgal.org/2025/05/22/Surface_remeshing/
      * https://doc.cgal.org/latest/PMP_Remeshing/index.html#Chapter_PMPRemeshing
      * `approximated_centroidal_Voronoi_diagram_remeshing()`
      * `surface_Delaunay_remeshing()`
      * `PMP::remesh_planar_patches()`
      * `PMP::remesh_almost_planar_patches()`
  * bounding meshes
      * approximate bounding ellipsoid
      * bounding spheres
  * `fill_boundary_hole()`
      * pass more parameters (small holes)
  * clip
  * isosurfacing https://doc.cgal.org/latest/Isosurfacing_3/
      * Marching cubes
      * Topologically correct marching cubes
      * Dual contouring
  * color wash for mesh distances as in `Rvcg::vcgMetro()`
