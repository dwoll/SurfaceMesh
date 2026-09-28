# TODO

  * `removeSelfIntersections()` does not work for example

# Wishlist

  * https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/Advancing_front_surface_reconstruction_2reconstruction_structured_8cpp-example.html
  * remeshing
      * https://doc.cgal.org/latest/PMP_Remeshing/
      * `approximated_centroidal_Voronoi_diagram_remeshing()`
      * `surface_Delaunay_remeshing()`
      * `PMP::remesh_planar_patches()`
      * `PMP::remesh_almost_planar_patches()`
  * bounding meshes
      * approximate bounding ellipsoid
      * bounding spheres
  * `fill_boundary_hole()`
      * pass more parameters (small holes)
  * clip (plane, ...)
  * isosurfacing https://doc.cgal.org/latest/Isosurfacing_3/
      * Marching cubes
      * Topologically correct marching cubes
      * Dual contouring
  * color wash for mesh distances as in `Rvcg::vcgMetro()`
  * alpha shapes
  * `assignFaceColors()`, `assignVertexColors()`, `assignVertexNormals()`
      * currently, there is no value in bringing these property maps to C++ side as nothing is done with them, except later exporting back to R side
      * `make_surf_mesh()` etc. need to copy face colors, vertex colors, vertex normals
      * support alpha for vertex colors, face colors in `set*()`
      * mesh operations invalidate face colors, vertex colors
      * add to vignette vertex colors, face colors, assigned normals
