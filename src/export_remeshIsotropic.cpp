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

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>

// ----------------------------------------------------------------------- //
// code adapted from
// https://doc.cgal.org/latest/PMP_Remeshing/PMP_Remeshing_2isotropic_remeshing_example_8cpp-example.html
// ----------------------------------------------------------------------- //
struct halfedge2edge {
  halfedge2edge(const Mesh3& m, std::vector<dg_dscrptr>& edges)
    : m_mesh(m), m_edges(edges)
  {}

  void operator()(const hlfdg_dscrptr& h) const {
    m_edges.push_back(edge(h, m_mesh));
  }

  const Mesh3& m_mesh;
  std::vector<dg_dscrptr>& m_edges;
};

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsoUniform_cpp(
    const Rcpp::List rmesh,
    const double targetEdgeLen,
    const unsigned int nIter,
    const unsigned int nRelaxSteps,
    const bool protectConstraints,
    const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose

    std::vector<dg_dscrptr> border;
    PMP::Uniform_sizing_field<Mesh3> sizing_field(targetEdgeLen, mesh);
    CGAL::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, targetEdgeLen, mesh);
    PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                             CGAL::parameters::number_of_iterations(nIter)
                                 .number_of_relaxation_steps(nRelaxSteps)
                                 .protect_constraints(protectConstraints));

    mesh.collect_garbage();
    // remeshing requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List remeshIsoAdapt_cpp(
    const Rcpp::List rmesh,
    const double tol,
    const double edgeMin,
    const double edgeMax,
    const unsigned int nIter,
    const unsigned int nRelaxSteps,
    const bool protectConstraints,
    const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        true,        // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    const std::pair edge_min_max{ edgeMin, edgeMax };
    PMP::Adaptive_sizing_field<Mesh3> sizing_field(
        tol,
        edge_min_max,
        faces(mesh),
        mesh);
    PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                             CGAL::parameters::number_of_iterations(nIter)
                                 .number_of_relaxation_steps(nRelaxSteps)
                                 .protect_constraints(protectConstraints));

    mesh.collect_garbage();
    // remeshing requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
