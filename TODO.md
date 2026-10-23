# TODO

# Wishlist

  * `meshEllipsoid()` check using dirs directly
  * read, write more supported file formats
  * read, write normals
  * `removeSelfIntersections()` does not work for example
  * `toRGL()`, `makeMesh()`, `makeMeshValid()`
      *  handle switching y <-> z coordinates between CGAL, rgl conventions?
  * face normals
  * https://doc.cgal.org/latest/Polygon_mesh_processing/Polygon_mesh_processing_2interpolated_corrected_curvatures_PH_8cpp-example.html#a5
  * https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/Advancing_front_surface_reconstruction_2reconstruction_structured_8cpp-example.html
  * remeshing: https://doc.cgal.org/latest/PMP_Remeshing/
      * `PMP::remesh_planar_patches()`
      * `PMP::remesh_almost_planar_patches()`
  * `fill_boundary_hole()`
      * pass more parameters (small holes)
  * clip (plane, ...)
  * isosurfacing https://doc.cgal.org/latest/Isosurfacing_3/
      * Marching cubes
      * Topologically correct marching cubes
      * Dual contouring
  * color wash for mesh distances as in `Rvcg::vcgMetro()`
  * alpha shapes
