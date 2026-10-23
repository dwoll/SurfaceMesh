// ----------------------------------------------------------------------- //
// code adapted from
// https://doc.cgal.org/latest/PMP_Remeshing/PMP_Remeshing_2acvd_remeshing_example_8cpp-example.html
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/approximated_centroidal_Voronoi_diagram_remeshing.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

// ----------------------------------------------------------------------- //
// const double facetSize,
// [[Rcpp::export]]
Rcpp::List remeshACVD_cpp(
    const Rcpp::List rmesh,
    const Rcpp::String method,
    const unsigned int nb_clusters,
    const double vert_count_ratio,
    const double gradation_factor,
    const bool normals) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3, Vector3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose

    if(PMP::does_self_intersect(mesh)) {
        Rcpp::stop("`x` has self-intersections.");
    }

    // ensure input is single connected component
    auto fcm = mesh.add_property_map<Mesh3::Face_index, std::size_t>("fcm", 0).first;
    std::size_t nbcc = PMP::connected_components(mesh, fcm);
    if(nbcc != 1) {
      Rcpp::stop("The input is not a single connected component mesh.");
    }

    if(method == "uniform_iso") {
        // Uniform Isotropic ACVD
        PMP::approximated_centroidal_Voronoi_diagram_remeshing(
            mesh, nb_clusters,
            CGAL::parameters::vertex_count_ratio(vert_count_ratio));
    } else if(method == "uniform_iso_qem") {
        // With Post-Processing QEM Optimization
        PMP::approximated_centroidal_Voronoi_diagram_remeshing(
            mesh, nb_clusters,
            CGAL::parameters::use_postprocessing_qem(true)
                             .vertex_count_ratio(vert_count_ratio));
    } else if(method == "uniform_qem") {
        // With QEM Energy Minimization
        PMP::approximated_centroidal_Voronoi_diagram_remeshing(
            mesh, nb_clusters,
            CGAL::parameters::use_qem_based_energy(true)
                             .vertex_count_ratio(vert_count_ratio));
    } else if(method == "adapt_iso") {
        // Adaptive Isotropic ACVD
        PMP::approximated_centroidal_Voronoi_diagram_remeshing(
            mesh, nb_clusters,
            CGAL::parameters::gradation_factor(gradation_factor)
                           .vertex_count_ratio(vert_count_ratio));
    } else {
        Rcpp::stop("Wrong method.");
    }

    // remeshing requires triangle mesh -> output is triangle
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
