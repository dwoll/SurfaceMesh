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

#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <locale>     // tolower()
#include <filesystem> // path()

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
Rcpp::List faces_to_list(std::vector<std::vector<std::size_t>> &faces) {
    const std::size_t nFaces = faces.size();
    Rcpp::List face_list(nFaces);
    for(std::size_t i = 0; i < nFaces; i++) {
      const std::vector<std::size_t> face_i = faces[i];
      Rcpp::IntegerVector col_i(face_i.begin(), face_i.end());
      face_list(i) = col_i + 1;  // vectorized + 1?
    }
    return face_list;
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List readFileSoup_cpp(const std::string filename, const bool verbose) {
    std::vector<Point3> points;
    std::vector<std::vector<std::size_t>> faces;
    const bool ok = CGAL::IO::read_polygon_soup(
        filename, points, faces, CGAL::parameters::verbose(verbose));
    if(!ok) {
      Rcpp::stop("Reading failure.");
    }
    Rcpp::NumericMatrix vertex_mat = points3_to_matrix<K, Point3>(points);
    Rcpp::List          face_list  = faces_to_list(faces);
    return Rcpp::List::create(
        Rcpp::Named("vertices") = Rcpp::transpose(vertex_mat),
        Rcpp::Named("faces")    = face_list);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List readFileMesh_cpp(
  const std::string filename, const bool normals, const bool verbose) {
  Mesh3 mesh;
  const bool ok = CGAL::IO::read_polygon_mesh(
      filename, mesh, CGAL::parameters::verbose(verbose));
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  if(verbose) {
    run_mesh_checks<Mesh3>(mesh);
  }
  return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, normals);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
void writeFile_cpp(const std::string filename,
                   const bool binary,
                   const int precision,
                   const Rcpp::NumericMatrix vertices,
                   const Rcpp::List faceList) {
  // using P3V3 = std::pair<Point3, Vector3>;  // Point3 with normal Vector3
  // std::vector<P3V3> points_nv;
  // if(rmesh.containsElementNamed("normals")) {
  //     const Rcpp::NumericMatrix rnormals = Rcpp::as<Rcpp::NumericMatrix>(rmesh["normals"]);
  //     const std::vector<VectorT> vnormals = matrix_to_points3<VectorT>(rnormals);
  //     set_vnormals<MeshT, VectorT>(mesh, vnormals);
  // }
  // points_nv.push_back(P3V3(Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0)));
  // CGAL::IO::write_PLY(filename, points_nv, faces.first,
  //                     CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3V3>())
  //                                     .normal_map(CGAL::Second_of_pair_property_map<P3V3>()));

  const std::vector<Point3> points = matrix_to_points3<Point3>(vertices);
  const std::pair<std::vector<std::vector<std::size_t>>, bool> faces =
      list_to_faces2(faceList);
  if(filename.length() < 5) {
      Rcpp::stop("`filename` needs at least 5 characters, including dot file-extension.");
  }
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  if(ext == ".ply") {
    ok = CGAL::IO::write_PLY(
      filename, points, faces.first,
      CGAL::parameters::use_binary_mode(binary).stream_precision(precision));
  } else if(ext == ".stl") {
    if(!faces.second) {
      Rcpp::stop("STL files only accept triangular faces.");
    }
    ok = CGAL::IO::write_STL(
      filename, points, faces.first,
      CGAL::parameters::use_binary_mode(binary).stream_precision(precision));
  } else if(ext == ".obj") {
    ok = CGAL::IO::write_OBJ(
      filename, points, faces.first,
      CGAL::parameters::stream_precision(precision));
  } else if(ext == ".off") {
    ok = CGAL::IO::write_OFF(
      filename, points, faces.first,
      CGAL::parameters::stream_precision(precision));
  } else {
    Rcpp::stop("Unknown file extension.");
  }
  if(!ok) {
    Rcpp::stop("Failed to write file.");
  }
}
