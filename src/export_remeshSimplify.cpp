// ----------------------------------------------------------------------- //
// Daniel Wollschlaeger
// https://doc.cgal.org/latest/Surface_mesh_simplification/examples.html
// License: GPL-3
// ----------------------------------------------------------------------- //

#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::GarlandHeckbert_plane_policies<Mesh3, K>                  Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_plane_policies<Mesh3, K>    Prob_plane;
typedef SMS::GarlandHeckbert_triangle_policies<Mesh3, K>               Classic_tri;
typedef SMS::GarlandHeckbert_probabilistic_triangle_policies<Mesh3, K> Prob_tri;
typedef SMS::GarlandHeckbert_plane_and_line_policies<Mesh3, K>         Plane_and_line;

// ----------------------------------------------------------------------- //
// cost strategy: Lindstrom-Turk
// [[Rcpp::export]]
Rcpp::List simplifyLT_cpp(const Rcpp::List rmesh,
                          const bool repairSoup,
                          const bool repairMesh,
                          const Rcpp::String method,  // LT-R (ratio), LT-C (count), LT-BNCF bounded normal change filter
                          const double ueRatio,       // undirected edge count ratio relative to start, 0.1
                          const unsigned int ueCount, // undirected edge count, 1000 or num_edges(mesh)/2 - 1
                          const bool normals,
                          const bool verbose) {
  typedef SMS::LindstromTurk_placement<Mesh3> Placement;
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
    rmesh,
    true,        // soup
    true,        // triangulate - must be triangle
    repairSoup,  // repair_soup
    verbose);    // verbose

  if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("Mesh must be triangle.");
  }
  if((ueRatio <= 0.0) || (ueRatio >= 1.0)) {
      Rcpp::stop("ueRatio not in (0, 1).");
  }
  if(ueCount >= num_edges(mesh)) {
      Rcpp::stop("ueCount >= number of edges.");
  }

  int r;
  if(method == "LT-R") {
      // simplification stops when the ratio of undirected edges
      // left in the surface mesh relative to the input drops below the specified number
      SMS::Edge_count_ratio_stop_predicate<Mesh3> stop(ueRatio);
      r = SMS::edge_collapse(mesh, stop);
  } else if(method == "LT-C") {
      // simplification stops when the number of undirected edges
      // left in the surface mesh drops below the specified number
      SMS::Edge_count_stop_predicate<Mesh3> stop(ueCount);
      r = SMS::edge_collapse(mesh, stop);
  } else if(method == "LT-BNCF") {
      // Bounded_normal_change_filter checks if a placement would invert the normal of a face around
      // the stars of the two vertices of an edge that is candidate for an edge collapse.
      // It then rejects this placement
      SMS::Edge_count_stop_predicate<Mesh3> stop(ueCount);
      SMS::Bounded_normal_change_filter<> Filter;
      r = SMS::edge_collapse(mesh, stop,
                             CGAL::parameters::get_cost(SMS::LindstromTurk_cost<Mesh3>())
                                               .filter(Filter)
                                               .get_placement(Placement()));
  } else {
      Rcpp::stop("Wrong Lindstrom-Turk method.");
  }

  // report
  if(verbose) {
      std::string msg;
      std::string method_str = method;
      msg = std::to_string(r) + " edges removed with method " + method_str + ".";
      rmessage(msg);
  }
  // update faces
  mesh.collect_garbage();

  // surface mesh simplification does not guarantee that result has no self intersections
  if(repairMesh) {
      if(PMP::does_self_intersect(mesh)) {
          Mesh3 mesh_tmp = remove_selfint_mesh<K, Mesh3, Point3>(mesh, 1, verbose);
          mesh = std::move(mesh_tmp);
      }
  }
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// cost strategy: Garland-Heckbert
template <typename GHPolicies, typename MeshT>
int collapse_gh(MeshT& mesh, const double ueRatio) {
  SMS::Edge_count_ratio_stop_predicate<MeshT> stop(ueRatio);

  // Garland & Heckbert simplification policies
  typedef typename GHPolicies::Get_cost       GH_cost;
  typedef typename GHPolicies::Get_placement  GH_placement;
  typedef SMS::Bounded_normal_change_filter<> Filter;

  GHPolicies gh_policies(mesh);
  const GH_cost& gh_cost = gh_policies.get_cost();
  const GH_placement& gh_placement = gh_policies.get_placement();
  Filter filter;
  int r = SMS::edge_collapse(mesh, stop,
                             CGAL::parameters::get_cost(gh_cost)
                                 .filter(filter)
                                 .get_placement(gh_placement));

  return r;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List simplifyGH_cpp(const Rcpp::List rmesh,
                          const bool repairSoup,
                          const bool repairMesh,
                          const double ueRatio,      // undirected edge count ratio relative to start, 0.1, 0.2
                          const Rcpp::String policy, // "CP", "CT", "PP", "PT", "PL"
                          const bool normals,
                          const bool verbose) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
    rmesh,
    true,        // soup
    true,        // triangulate - must be triangle
    repairSoup,  // repair_soup
    verbose);    // verbose

  if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("Mesh must be triangle.");
  }
  if((ueRatio <= 0.0) || (ueRatio >= 1.0)) {
      Rcpp::stop("ueRatio not in (0, 1).");
  }
  int r;
  if(policy == "CP") {
    r = collapse_gh<Classic_plane, Mesh3>(mesh, ueRatio);
  } else if(policy == "CT") {
    r = collapse_gh<Classic_tri, Mesh3>(mesh, ueRatio);
  } else if(policy == "PP") {
    r = collapse_gh<Prob_plane, Mesh3>(mesh, ueRatio);
  } else if(policy == "PT") {
    r = collapse_gh<Prob_tri, Mesh3>(mesh, ueRatio);
  } else if(policy == "PL") {
    r = collapse_gh<Plane_and_line, Mesh3>(mesh, ueRatio);
  } else {
      Rcpp::stop("Wrong Garland-Heckbert policy.");
  }
  // update faces
  mesh.collect_garbage();

  // report
  if(verbose) {
      std::string msg;
      std::string policy_str = policy;
      msg = std::to_string(r) + " edges removed with method GH and policy " + policy_str + ".";
      rmessage(msg);
  }

  // surface mesh simplification does not guarantee that result has no self intersections
  if(repairMesh) {
      if(PMP::does_self_intersect(mesh)) {
          Mesh3 mesh_tmp = remove_selfint_mesh<K, Mesh3, Point3>(mesh, 1, verbose);
          mesh = std::move(mesh_tmp);
      }
  }
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}
