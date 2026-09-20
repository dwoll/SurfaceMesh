// ----------------------------------------------------------------------- //
// Code adapted from packages
// https://github.com/stla/Boov/
// https://github.com/stla/PolygonSoup/
// https://github.com/stla/cgalMeshes/
// developed and copyright by
// Stéphane Laurent <laurent_step@outlook.fr>
// adapted by
// Daniel Wollschlaeger
// License: GPL-3
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <cmath>

// ----------------------------------------------------------------------- //
// initial mesh generation - EPIC kernel - TODO make parameter
// [[Rcpp::export]]
Rcpp::List makeMesh_cpp(const Rcpp::List rmesh,
                        const bool triangulate,
                        const bool repairSoup,
                        const bool removeIntersections,
                        const int removeMethod,
                        const bool fillHoles,
                        const bool fairHole,
                        const unsigned int maxNumHoles,
                        const bool normals,
                        const bool verbose) {
  if(verbose) { rmessage("Processing mesh..."); }
  Mesh3 mesh = make_surf_mesh<K, Mesh3, Point3>(
      rmesh,
      triangulate,         // triangulate
      repairSoup,          // repair_soup
      removeIntersections, // remove_intersections
      removeMethod,        // remove_method
      fillHoles,           // fill_holes
      fairHole,            // fair hole
      maxNumHoles,         // max_num_holes
      verbose);            // verbose
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// initial mesh generation assuming valid input - EPIC kernel - TODO make parameter
// [[Rcpp::export]]
Rcpp::List makeMeshFF_cpp(const Rcpp::String filename,
                          const bool triangulate,
                          const bool repairSoup,
                          const bool removeIntersections,
                          const int removeMethod,
                          const bool fillHoles,
                          const bool fairHole,
                          const unsigned int maxNumHoles,
                          const bool normals,
                          const bool verbose) {
  if(verbose) { rmessage("Processing mesh..."); }
  Mesh3 mesh = make_surf_mesh_ff<K, Mesh3, Point3>(
      filename,
      triangulate,         // triangulate
      repairSoup,          // repair_soup
      removeIntersections, // remove_intersections
      removeMethod,        // remove_method
      fillHoles,           // fill_holes
      fairHole,            // fair hole
      maxNumHoles,         // max_num_holes
      verbose);            // verbose
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// initial mesh generation assuming valid input - EPIC kernel - TODO make parameter
// [[Rcpp::export]]
Rcpp::List makeMeshValid_cpp(const Rcpp::List rmesh,
                             const bool soup,
                             const bool triangulate,
                             const bool normals,
                             const bool verbose) {
  if(verbose) { rmessage("Processing mesh..."); }
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      soup,
      triangulate,
      false,       // repairSoup
      verbose);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// initial mesh generation assuming valid input - EPIC kernel - TODO make parameter
// [[Rcpp::export]]
Rcpp::List makeMeshValidFF_cpp(const Rcpp::String filename,
                               const bool soup,
                               const bool triangulate,
                               const bool normals,
                               const bool verbose) {
  if(verbose) { rmessage("Processing mesh..."); }
  Mesh3 mesh = make_surf_mesh_valid_ff<Mesh3, Point3>(
      filename,
      soup,
      triangulate,
      false,     // repair_soup
      verbose);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List addVNormals_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
    rmesh,
    false,       // soup
    false,       // triangulate
    false,       // repair_soup
    false);      // verbose
 return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, true);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesBoundVolume_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   if(!CGAL::is_closed(mesh)) {
      Rcpp::warning("Mesh is not closed.");
      return false;
  }
  if(PMP::does_self_intersect(mesh)) {
      Rcpp::warning("Mesh has self-intersections.");
      return false;
  }
  return PMP::does_bound_a_volume(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool doesSelfIntersect_cpp(
  const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   return PMP::does_self_intersect(mesh);
}

// ----------------------------------------------------------------------- //
// use EPEC kernel for fill_boundary_holes()
// [[Rcpp::export]]
Rcpp::List fillBoundaryHoles_cpp(
  const Rcpp::List rmesh,
  const bool fairHole,
  const unsigned int maxNumHoles,
  const bool normals,
  const bool verbose) {
  EMesh3 mesh = make_surf_mesh<EK, EMesh3, EPoint3>(
      rmesh,
      true,         // triangulate - must be triangle
      true,         // repair_soup
      false,        // remove_intersections
      1,            // remove_method
      true,         // fill_holes
      fairHole,     // fair hole
      maxNumHoles,  // max_num_holes
      verbose);     // verbose
   return get_rmesh<EK, EMesh3, EPoint3, EVector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getArea_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   if(PMP::does_self_intersect(mesh)) {
    Rcpp::warning("The mesh self-intersects.");
    return Rcpp::NumericVector::get_na();
  }
  const K::FT a = PMP::area(mesh);
  return CGAL::to_double<K::FT>(a);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List getBoundingBox_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      false,       // triangulate
      false,       // repair_soup
      false);      // verbose
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  Rcpp::NumericVector lcorner = { bbox.xmin(), bbox.ymin(), bbox.zmin() };
  Rcpp::NumericVector ucorner = { bbox.xmax(), bbox.ymax(), bbox.zmax() };
  return Rcpp::List::create(
    Rcpp::Named("lcorner") = lcorner,
    Rcpp::Named("ucorner") = ucorner);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List getBoundingBoxOptimal_cpp(
  const Rcpp::List rmeshIn, const bool triangulate, const bool normals) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmeshIn,
      false,       // soup
      false,       // triangulate
      false,       // repair_soup
      false);      // verbose
   std::array<Point3, 8> obb_pts;
  CGAL::oriented_bounding_box(mesh, obb_pts,
                              CGAL::parameters::use_convex_hull(true));
  // make mesh out of oriented bounding box
  Mesh3 obb_mesh;
  CGAL::make_hexahedron(
    obb_pts[0], obb_pts[1], obb_pts[2], obb_pts[3],
    obb_pts[4], obb_pts[5], obb_pts[6], obb_pts[7],
    obb_mesh);
  Rcpp::List rmesh_obb = get_rmesh<K, Mesh3, Point3, Vector3>(obb_mesh, triangulate, normals);
  Rcpp::NumericMatrix hex_verts(3, 8);
  for(int i = 0; i < 8; i++) {
    Point3 pt = obb_pts[i];
    Rcpp::NumericVector v =
      Rcpp::NumericVector::create(pt.x(), pt.y(), pt.z());
    hex_verts(Rcpp::_, i) = v;
  }
  return Rcpp::List::create(
    Rcpp::Named("mesh")       = rmesh_obb,
    Rcpp::Named("hxVertices") = hex_verts);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericVector getCentroid_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   Rcpp::NumericVector ctr(3);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::warning("The mesh is not triangle.");
      ctr(0) = Rcpp::NumericVector::get_na();
      ctr(1) = Rcpp::NumericVector::get_na();
      ctr(2) = Rcpp::NumericVector::get_na();
  } else {
      const Point3 centroid = PMP::centroid(mesh);
      ctr(0) = CGAL::to_double<K::FT>(centroid.x());
      ctr(1) = CGAL::to_double<K::FT>(centroid.y());
      ctr(2) = CGAL::to_double<K::FT>(centroid.z());
  }
  return ctr;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List getConvexHull_cpp(const Rcpp::NumericMatrix rpoints, const bool normals) {
  const std::size_t nPts = rpoints.ncol();
  std::vector<Point3> points;
  points.reserve(nPts);
  for(std::size_t i = 0; i < nPts; i++) {
    Rcpp::NumericVector pt = rpoints(Rcpp::_, i);
    points.emplace_back(Point3(pt(0), pt(1), pt(2)));
  }
  Mesh3 mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericVector getDistance_cpp(
    const Rcpp::List rmesh, const Rcpp::NumericMatrix rpoints) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh3> Primitive;
  typedef CGAL::AABB_traits_3<K, Primitive> Tree_Traits;
  typedef CGAL::AABB_tree<Tree_Traits> Tree;

  const std::size_t nPts = rpoints.ncol();
  Rcpp::NumericVector distances(nPts);
  if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::warning("The mesh is not triangle.");
      for(std::size_t i = 0; i < nPts; i++) {
          distances(i) = Rcpp::NumericVector::get_na();
      }
  } else {
      Tree tree(faces(mesh).first, faces(mesh).second, mesh);
      for(std::size_t i = 0; i < nPts; i++) {
          Rcpp::NumericVector point_i = rpoints(Rcpp::_, i);
          const Point3 pt = Point3(point_i(0), point_i(1), point_i(2));
          double dsq = CGAL::to_double<typename K::FT>(tree.squared_distance(pt));
          distances(i) = std::sqrt(dsq);
      }
  }
  return distances;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
double getVolume_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   if(!CGAL::is_closed(mesh)) {
    Rcpp::warning("The mesh is not closed.");
    return Rcpp::NumericVector::get_na();
  }
  if(PMP::does_self_intersect(mesh)) {
    Rcpp::warning("The mesh self-intersects.");
    return Rcpp::NumericVector::get_na();
  }
  const K::FT vol = PMP::volume(mesh);
  return CGAL::to_double<K::FT>(vol);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool hasGarbage_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      false,       // triangulate
      false,       // repair_soup
      false);      // verbose
  return mesh.has_garbage();
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isClosed_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      false,       // triangulate
      false,       // repair_soup
      false);      // verbose
   return CGAL::is_closed(mesh);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
bool isValid_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      false,       // triangulate
      false,       // repair_soup
      false);      // verbose
  return mesh.is_valid(false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List orientToBoundVolume_cpp(
  const Rcpp::List rmesh, const bool normals) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
      rmesh,
      false,       // soup
      true,        // triangulate - must be triangle
      false,       // repair_soup
      false);      // verbose
   if(!CGAL::is_triangle_mesh(mesh)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  PMP::orient_to_bound_a_volume(mesh);
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// use EPEC kernel for autorefine_triangle_soup()
// [[Rcpp::export]]
Rcpp::List removeSelfIntersections_cpp(
  const Rcpp::List rmesh,
  const int method,
  const bool normals,
  const bool verbose) {
  EMesh3 mesh = make_surf_mesh<EK, EMesh3, EPoint3>(
      rmesh,
      true,        // triangulate - must be triangle
      true,        // repair_soup
      true,        // remove_intersections
      method,      // remove_method
      false,       // fill_holes
      false,       // fair hole
      0,           // max_num_holes
      verbose);    // verbose
   return get_rmesh<EK, EMesh3, EPoint3, EVector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::NumericMatrix samplePoints_cpp(const Rcpp::List rmesh, const Rcpp::List ropts) {
  Rcpp::List ropts_l = Rcpp::as<Rcpp::List>(ropts);
  sample_opts opts = ropts_to_sample_opts(ropts_l);
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
    rmesh,
    false,       // soup
    true,        // triangulate - must be triangle
    false,       // repair_soup
    false);      // verbose
  std::vector<Point3> points;
  sample_points<Mesh3, Point3>(mesh, points, opts);
  const Rcpp::NumericMatrix rpoints = points3_to_matrix<K, Point3>(points);
  return Rcpp::transpose(rpoints);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List triangulateMesh_cpp(const Rcpp::List rmesh, const bool normals) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
    rmesh,
    false,       // soup
    true,        // triangulate - must be triangle
    false,       // repair_soup
    false);      // verbose
 return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
