# TODO

  * `toRGL()`, `makeMesh()`, `makeMeshValid()`
      *  handle switching y <-> z coordinates between CGAL, rgl conventions?
  * `meshEllipsoid()` check using dirs directly
  * read, write more supported file formats
  * read, write normals
  * `setVertexNormals()` check normals are unit vectors
  * `removeSelfIntersections()` does not work for example

# Wishlist

  * functions that accept a matrix should also accept a `CGALmesh` object, `mesh3d` object and extract vertices - via S3 methods
  * face normals
  * https://doc.cgal.org/latest/Polygon_mesh_processing/Polygon_mesh_processing_2interpolated_corrected_curvatures_PH_8cpp-example.html#a5
  * https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/Advancing_front_surface_reconstruction_2reconstruction_structured_8cpp-example.html
  * remeshing: https://doc.cgal.org/latest/PMP_Remeshing/
      * `approximated_centroidal_Voronoi_diagram_remeshing()`
      * `surface_Delaunay_remeshing()`
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
