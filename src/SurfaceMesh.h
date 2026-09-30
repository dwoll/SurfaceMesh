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
#define _CGALMESHHEADER_

// -------------------------------------------------------------------------- //
#include <Rcpp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Vector_3.h>

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
#define CGAL_EIGEN3_ENABLED 1
#define PIA_TAG CGAL::Parallel_if_available_tag
#define SEQ_TAG CGAL::Sequential_tag

// -------------------------------------------------------------------------- //
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel   EK;

typedef K::Point_3  Point3;
typedef EK::Point_3 EPoint3;

typedef K::Vector_3  Vector3;
typedef EK::Vector_3 EVector3;

typedef CGAL::Surface_mesh<Point3>  Mesh3;
typedef CGAL::Surface_mesh<EPoint3> EMesh3;

// -------------------------------------------------------------------------- //
// triangle sample options for sample_points()
// PMP::parameters::use_random_uniform_sampling(true)     // true
// PMP::parameters::use_grid_sampling(false)              // false
// PMP::parameters::use_monte_carlo_sampling(false)       // false
// PMP::parameters::do_sample_vertices(true)              // true
// PMP::parameters::do_sample_edges(true)                 // true
// PMP::parameters::do_sample_faces(true)                 // true
// PMP::parameters::grid_spacing(n)                       // double, for grid sampling
// PMP::parameters::number_of_points_on_edges(n)          //  unsigned int, for random sampling
// PMP::parameters::number_of_points_on_faces(n)          // *unsigned int, for random sampling
// PMP::parameters::number_of_points_per_distance_unit(n) // double, for random sampling and Monte Carlo sampling
// PMP::parameters::number_of_points_per_area_unit(n)     // double, for random sampling and Monte Carlo sampling
// PMP::parameters::number_of_points_per_edge(n)          // unsigned int, for Monte-Carlo sampling
// PMP::parameters::number_of_points_per_face(n)          // unsigned int, for Monte-Carlo sampling
struct sample_opts {
    unsigned int method; // 1: random uniform, 2: grid, 3: Monte Carlo
    bool sampleVerts;
    bool sampleEdges;
    bool sampleFaces;
    double gridSpacing;
    unsigned int ptsOnEdges;
    unsigned int ptsOnFaces;
    double ptsPerDist;
    unsigned int ptsPerEdge;
    double ptsPerArea;
    unsigned int ptsPerFace;
};

// -------------------------------------------------------------------------- //
template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix&);

template <typename KernelT, typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT>&);

template <typename KernelT, typename MeshT, typename VectorT>
std::vector<VectorT> compute_vnormals(MeshT);

template <typename KernelT, typename MeshT, typename VectorT>
std::optional<std::vector<VectorT>> get_vnormals(const MeshT&);

template <typename MeshT, typename VectorT>
void set_vnormals(MeshT&, const std::vector<VectorT>&);

template <typename KernelT, typename MeshT, typename PointT>
MeshT soup_to_mesh(
    std::vector<PointT>&,                   // points
    std::vector<std::vector<std::size_t>>&, // faces
    const bool,                             // triangulate
    const bool,                             // repair_soup
    const bool,                             // remove_intersections
    const int,                              // remove_method
    const bool,                             // fill_holes
    const bool,                             // fair_hole
    const unsigned int,                     // max_num_holes
    const bool);                            // verbose

template <typename KernelT, typename MeshT, typename PointT>
MeshT make_surf_mesh(
    const Rcpp::List&,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

template <typename KernelT, typename MeshT, typename PointT>
MeshT make_surf_mesh_ff(
    const Rcpp::String,
    const bool,
    const bool,
    const bool,
    const int,
    const bool,
    const bool,
    const unsigned int,
    const bool);

template <typename MeshT, typename PointT, typename VectorT>
MeshT make_surf_mesh_valid(
    const Rcpp::List&, const bool, const bool, const bool, const bool);

template <typename MeshT, typename PointT>
MeshT make_surf_mesh_valid_ff(
    const Rcpp::String, const bool, const bool, const bool, const bool);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh1(const MeshT&, const bool);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh2(const MeshT&, const bool, const std::size_t);

template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List get_rmesh(MeshT&, const bool, const bool);

template <typename KernelT, typename MeshT, typename PointT>
MeshT remove_selfint_mesh(const MeshT&, const int, const bool);

template <typename MeshT, typename PointT>
MeshT fill_boundary_holes(
    MeshT&, const bool, const double, const int, const unsigned int, const bool);

template <typename MeshT, typename VectorT>
void remove_properties(MeshT&, const std::vector<std::string>&);

template <typename MeshT>
void run_mesh_checks(const MeshT&);

template <typename MeshT, typename PointT>
void sample_points(const MeshT&, std::vector<PointT>&, const sample_opts&);

template <typename KernelT, typename MeshT, typename PointT>
std::tuple<double, double, double> get_metro(
    const MeshT&,
    const MeshT&,
    std::vector<double>&,
    std::vector<double>&,
    const bool,
    const double,
    const sample_opts&);

// -------------------------------------------------------------------------- //
// no template
sample_opts ropts_to_sample_opts(const Rcpp::List&);

void rmessage(std::string);

          std::vector<std::vector<std::size_t>>        list_to_faces1(const Rcpp::List&);
std::pair<std::vector<std::vector<std::size_t>>, bool> list_to_faces2(const Rcpp::List&);

// -------------------------------------------------------------------------- //
#endif
