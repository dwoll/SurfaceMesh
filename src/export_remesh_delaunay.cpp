// ----------------------------------------------------------------------- //
// code adapted from
// https://doc.cgal.org/latest/PMP_Remeshing/PMP_Remeshing_2delaunay_remeshing_example_8cpp-example.html
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Mesh_constant_domain_field_3.h>

// ----------------------------------------------------------------------- //
// const double facetSize,
// [[Rcpp::export]]
Rcpp::List remeshDelaunay_cpp(
    const Rcpp::List rmesh,
    const double edgeSize,
    const double facetAngle,
    const double facetDist,
    const bool   protectConstraints,
    const double dihedralAngle,
    const bool   normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3, Vector3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose

    using Sizing_field = CGAL::Mesh_constant_domain_field_3<K, int>;
    using EIFMap = boost::property_map<Mesh3, CGAL::edge_is_feature_t>::type;
    if(PMP::does_self_intersect(mesh)) {
        Rcpp::stop("`x` has self-intersections.");
    }

    // Sizing_field size_facet(facetSize);
    Mesh3 mesh_out;
    Sizing_field size_edge(edgeSize);
    if(protectConstraints) {
        EIFMap eif = get(CGAL::edge_is_feature, mesh);
        PMP::detect_sharp_edges(mesh, dihedralAngle, eif);
        Mesh3 mesh_tmp = PMP::surface_Delaunay_remeshing(mesh,
                                                         CGAL::parameters::protect_constraints(true)
                                                                           .mesh_edge_size(size_edge)
                                                                           // .mesh_facet_size(size_facet)
                                                                           .mesh_facet_angle(facetAngle)
                                                                           .mesh_facet_distance(facetDist)
                                                                           .edge_is_constrained_map(eif));
         mesh_out = std::move(mesh_tmp);
    } else {
        Mesh3 mesh_tmp = PMP::surface_Delaunay_remeshing(mesh,
                                                         CGAL::parameters::protect_constraints(false)
                                                                           .mesh_edge_size(size_edge)
                                                                           // .mesh_facet_size(size_facet)
                                                                           .mesh_facet_angle(facetAngle)
                                                                           .mesh_facet_distance(facetDist));
         mesh_out = std::move(mesh_tmp);
    }

    // remeshing requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh_out, false, normals);
}
