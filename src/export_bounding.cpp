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

#include <CGAL/optimal_bounding_box.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List getBoundingBox_cpp(const Rcpp::List rmesh) {
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3, Vector3>(
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
  Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3, Vector3>(
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
// adapted from
// https://doc.cgal.org/latest/Bounding_volumes/Approximate_min_ellipsoid_d_2ellipsoid_8cpp-example.html
// [[Rcpp::export]]
Rcpp::List getBoundingEllipsoid_cpp(const Rcpp::NumericMatrix rpoints, const double eps) {
    using CK        = CGAL::Cartesian_d<double>;      // kernel
    using Traits    = CGAL::Approximate_min_ellipsoid_d_traits_d<CK, CGAL::MP_Float>;
    using Point_vec = std::vector<Traits::Point>;
    using AME       = CGAL::Approximate_min_ellipsoid_d<Traits>;

    std::vector<Point3> points = matrix_to_points3<Point3>(rpoints);
    Point_vec pts_ame;
    const int dim = 3;

    for(Point3 pt : points) {
        double px = CGAL::to_double<typename K::FT>(pt.x());
        double py = CGAL::to_double<typename K::FT>(pt.y());
        double pz = CGAL::to_double<typename K::FT>(pt.z());
        const double coords[] = {px, py, pz};
        pts_ame.emplace_back(dim, coords, coords + dim);
    }
    // compute approximation:
    Traits traits;
    AME ame(eps, pts_ame.begin(), pts_ame.end(), traits);
    AME::Center_coordinate_iterator cc_ib = ame.center_cartesian_begin();
    AME::Center_coordinate_iterator cc_ie = ame.center_cartesian_end();
    // dynamic growing performance ok
    std::vector<double> ctr;      // centroid
    std::vector<double> sa_lens;  // semi-axis lengths
    std::vector<Point3> sa_dirs;  // semi-axis direction
    // ellipsoid center
    for( ; cc_ib != cc_ie; ++cc_ib) {
        ctr.push_back(CGAL::to_double<K::FT>(*cc_ib));
    }

    // ellipsoid  axes
    AME::Axes_lengths_iterator al_ib = ame.axes_lengths_begin();
    for(unsigned int i=0; i < dim; ++i) {
        sa_lens.push_back(CGAL::to_double<K::FT>(*al_ib));
        al_ib++;
        AME::Axes_direction_coordinate_iterator adc_ib = ame.axis_direction_cartesian_begin(i);
        AME::Axes_direction_coordinate_iterator adc_ie = ame.axis_direction_cartesian_end(i);
        std::vector<double> dir;
        for( ; adc_ib != adc_ie; ++adc_ib) {
            dir.push_back(CGAL::to_double<K::FT>(*adc_ib));
        }
        sa_dirs.push_back(Point3(dir[0], dir[1], dir[2]));
    }

    Rcpp::NumericVector rctr     = Rcpp::wrap(ctr);      // no dynamic growing
    Rcpp::NumericVector rsa_lens = Rcpp::wrap(sa_lens);  // no dynamic growing
    Rcpp::NumericMatrix rsa_dirs = points3_to_matrix<K, Point3>(sa_dirs);
    return Rcpp::List::create(Rcpp::Named("center")     = rctr,
                              Rcpp::Named("lengths")    = rsa_lens,
                              Rcpp::Named("directions") = Rcpp::transpose(rsa_dirs));
}

// ----------------------------------------------------------------------- //
// adapted from
// https://doc.cgal.org/latest/Bounding_volumes/Min_sphere_d_2min_sphere_3_8cpp-example.html
// [[Rcpp::export]]
Rcpp::List getBoundingSphere_cpp(const Rcpp::NumericMatrix rpoints) {
    using Traits     = CGAL::Min_sphere_of_points_d_traits_3<K, double>;
    using Min_sphere = CGAL::Min_sphere_of_spheres_d<Traits>;
    std::vector<Point3> points = matrix_to_points3<Point3>(rpoints);

    Min_sphere  ms(points.begin(), points.end());       // smallest enclosing sphere
    Min_sphere::Cartesian_const_iterator ccib = ms.center_cartesian_begin();
    Min_sphere::Cartesian_const_iterator ccie = ms.center_cartesian_end();
    const double r = CGAL::to_double<typename K::FT>(ms.radius());
    std::vector<double> ctr;  // dynamic growing performance ok
    for( ; ccib != ccie; ++ccib) {
        ctr.push_back(CGAL::to_double<K::FT>(*ccib));
    }
    Rcpp::NumericVector rctr = Rcpp::wrap(ctr);  // no dynamic growing
    return Rcpp::List::create(Rcpp::Named("center") = rctr,
                              Rcpp::Named("radius") = r);
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
