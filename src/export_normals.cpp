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

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

// ----------------------------------------------------------------------- //
// EPIC kernel only
// [[Rcpp::export]]
Rcpp::NumericMatrix normals_jet_pca_cpp(const Rcpp::NumericMatrix pts,
                                        const unsigned int nNeighbors,
                                        const Rcpp::String method) {
  const std::size_t nPts = pts.ncol();
  std::vector<P3V3> points_wn(nPts);
  for(std::size_t i = 0; i < nPts; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    points_wn[i] = std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)),
                                  Vector3(0.0, 0.0, 0.0));
  }

  if(method == "jet") {
      CGAL::jet_estimate_normals<PIA_TAG>(
          points_wn, nNeighbors,
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3V3>())
                          .normal_map(CGAL::Second_of_pair_property_map<P3V3>()));
  } else if(method == "pca") {
      CGAL::pca_estimate_normals<PIA_TAG>(
          points_wn, nNeighbors,
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3V3>())
                          .normal_map(CGAL::Second_of_pair_property_map<P3V3>()));
  } else {
      Rcpp::stop("Wrong method.");
  }

  CGAL::mst_orient_normals(
      points_wn, nNeighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3V3>())
                      .normal_map(CGAL::Second_of_pair_property_map<P3V3>()));

  Rcpp::NumericMatrix normals_mat(3, nPts);
  for(std::size_t i = 0; i < nPts; i++) {
    Rcpp::NumericVector normal_i(3);
    const Vector3 normal = points_wn[i].second;
    normal_i(0) = normal.x();
    normal_i(1) = normal.y();
    normal_i(2) = normal.z();
    normals_mat(Rcpp::_, i) = normal_i;
  }

  return normals_mat;
}
