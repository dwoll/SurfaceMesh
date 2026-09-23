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

#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshSmoothShape_cpp(
  const Rcpp::List rmesh,
  const Rcpp::IntegerVector indices,
  const unsigned int nIter,
  const double time,
  const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    std::set<Mesh3::Vertex_index> constrained_vertices;
    for(Mesh3::Vertex_index v : vertices(mesh)) {
      if(is_border(v, mesh)) {
          constrained_vertices.insert(v);
      }
    }
    CGAL::Boolean_property_map<std::set<Mesh3::Vertex_index>> vcmap(constrained_vertices);
    const size_t nIdx = indices.size();
    if(nIdx == 0) {
        PMP::smooth_shape(mesh, time,
                          CGAL::parameters::number_of_iterations(nIter)
                          .vertex_is_constrained_map(vcmap));
    } else {
        std::list<fc_dscrptr> selectedFaces;
        const size_t nFaces = mesh.number_of_faces();
        for(std::size_t i = 0; i < nIdx; i++) {
          const size_t idx = indices(i);
          if(idx >= nFaces) {
            Rcpp::stop("Face index too large.");
          }
          // TODO use `CGAL::SM_Face_index(idx)`
          selectedFaces.push_back(*(mesh.faces().begin() + idx));
        }
        PMP::smooth_shape<Mesh3>(selectedFaces, mesh, time,
                                 PMP::parameters::number_of_iterations(nIter)
                                 .vertex_is_constrained_map(vcmap));
    }
    mesh.collect_garbage();
    // subdivision requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// https://doc.cgal.org/latest/PMP_Remeshing/PMP_Remeshing_2mesh_smoothing_example_8cpp-example.html
// [[Rcpp::export]]
Rcpp::List remeshSmoothAA_cpp(
  const Rcpp::List rmesh,
  const unsigned int nIter,
  const double dihedralAngle,
  const bool useAngleSmooth,
  const bool useAreaSmooth,
  const bool useSafeConstr,
  const bool useDelaunay,
  const bool doProject,
  const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    // constrain edges with a dihedral angle over given value
    typedef boost::property_map<Mesh3, CGAL::edge_is_feature_t>::type EIFMap;
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    PMP::detect_sharp_edges(mesh, dihedralAngle, eif);  // dihedralAngle = 60
    unsigned int sharp_counter = 0;
    for(dg_dscrptr e : edges(mesh)) {
        if(get(eif, e)) {
            ++sharp_counter;
        }
    }
    // smooth with both angle and area criteria + Delaunay flips
    PMP::angle_and_area_smoothing(mesh,
                                  CGAL::parameters::number_of_iterations(nIter)
                                                   .use_angle_smoothing(useAngleSmooth)
                                                   .use_area_smoothing(useAreaSmooth)
                                                   .use_safety_constraints(useSafeConstr) // false: authorize all moves
                                                   .use_Delaunay_flips(useDelaunay)
                                                   .do_project(doProject)
                                                   .edge_is_constrained_map(eif));

    mesh.collect_garbage();
    // subdivision requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// https://doc.cgal.org/latest/PMP_Remeshing/PMP_Remeshing_2tangential_relaxation_example_8cpp-example.html
// [[Rcpp::export]]
Rcpp::List remeshSmoothTR_cpp(
  const Rcpp::List rmesh,
  const unsigned int nIter,
  const bool relaxConstr,
  const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    PMP::tangential_relaxation(mesh,
                               CGAL::parameters::number_of_iterations(nIter)
                                                .relax_constraints(relaxConstr));
    // remeshing requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
