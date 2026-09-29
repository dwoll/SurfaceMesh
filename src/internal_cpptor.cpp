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

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for (E)Point3 and (E)Vector3
template <typename KernelT, typename PointT>
Rcpp::NumericMatrix points3_to_matrix(const std::vector<PointT> &points) {
    const std::size_t nPts = points.size();
    Rcpp::NumericMatrix pts_mat(3, nPts);
    for(std::size_t i = 0; i < nPts; i++) {
        Rcpp::NumericVector col_i(3);
        const PointT point = points[i];
        col_i(0) = CGAL::to_double<typename KernelT::FT>(point.x());
        col_i(1) = CGAL::to_double<typename KernelT::FT>(point.y());
        col_i(2) = CGAL::to_double<typename KernelT::FT>(point.z());
        pts_mat(Rcpp::_, i) = col_i;
    }
    return pts_mat;
}

template Rcpp::NumericMatrix points3_to_matrix<K,  Point3>(const   std::vector<Point3>&);
template Rcpp::NumericMatrix points3_to_matrix<EK, EPoint3>(const  std::vector<EPoint3>&);
template Rcpp::NumericMatrix points3_to_matrix<K,  Vector3>(const  std::vector<Vector3>&);
template Rcpp::NumericMatrix points3_to_matrix<EK, EVector3>(const std::vector<EVector3>&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::NumericMatrix get_vertices(const MeshT &mesh) {
    const std::size_t nVerts = mesh.number_of_vertices();
    Rcpp::NumericMatrix vertex_mat(3, nVerts);
    std::size_t i = 0;
    for(typename MeshT::Vertex_index vd : vertices(mesh)) {
        Rcpp::NumericVector col_i(3);
        const PointT vertex = mesh.point(vd);
        col_i(0) = CGAL::to_double<typename KernelT::FT>(vertex.x());
        col_i(1) = CGAL::to_double<typename KernelT::FT>(vertex.y());
        col_i(2) = CGAL::to_double<typename KernelT::FT>(vertex.z());
        vertex_mat(Rcpp::_, i) = col_i;
        i++;
    }
    return vertex_mat;
}

template Rcpp::NumericMatrix get_vertices<K,  Mesh3,  Point3>(const  Mesh3&);
template Rcpp::NumericMatrix get_vertices<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame get_edges(const MeshT &mesh) {
  const std::size_t nEdges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nEdges);
  Rcpp::IntegerVector I2(nEdges);
  Rcpp::NumericVector Length(nEdges);
  Rcpp::NumericVector Angle(nEdges);
  Rcpp::LogicalVector Exterior(nEdges);
  Rcpp::LogicalVector Coplanar(nEdges);
  std::size_t i = 0;
  for(typename MeshT::Edge_index ed : edges(mesh)) {
      typename MeshT::Vertex_index s    = source(ed, mesh);
      typename MeshT::Vertex_index t    = target(ed, mesh);
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
      I1(i) = static_cast<int>(s) + 1;
      I2(i) = static_cast<int>(t) + 1;
      std::vector<PointT> points(4);
      points[0] = mesh.point(s);
      points[1] = mesh.point(t);
      points[2] = mesh.point(mesh.target(mesh.next(h0)));
      points[3] = mesh.point(mesh.target(mesh.next(h1)));
      typename KernelT::FT len   = PMP::edge_length(h0, mesh);
      typename KernelT::FT angle = CGAL::abs(CGAL::approximate_dihedral_angle(
          points[0], points[1], points[2], points[3]));
      Length(i)   = CGAL::to_double<typename KernelT::FT>(len);
      Angle(i)    = CGAL::to_double<typename KernelT::FT>(angle);
      Exterior(i) = (angle < 179.0) || (angle > 181.0);
      Coplanar(i) = CGAL::coplanar(points[0], points[1], points[2], points[3]);
      i++;
  }
  Rcpp::DataFrame Edges = Rcpp::DataFrame::create(
      Rcpp::Named("i1")       = I1,
      Rcpp::Named("i2")       = I2,
      Rcpp::Named("length")   = Length,
      Rcpp::Named("angle")    = Angle,
      Rcpp::Named("exterior") = Exterior,
      Rcpp::Named("coplanar") = Coplanar);
  return Edges;
}

template Rcpp::DataFrame get_edges<K,  Mesh3,  Point3>(const  Mesh3&);
template Rcpp::DataFrame get_edges<EK, EMesh3, EPoint3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for possibly heterogeneous faces -> list
template <typename MeshT>
Rcpp::List get_faces1(const MeshT &mesh) {
  const std::size_t nFaces = mesh.number_of_faces();
  Rcpp::List face_list(nFaces);
  std::size_t i = 0;
  for(typename MeshT::Face_index fd : faces(mesh)) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      face_list(i) = col_i;
      i++;
  }
  return face_list;
}

template Rcpp::List get_faces1<Mesh3>(const  Mesh3&);
template Rcpp::List get_faces1<EMesh3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for homogeneous faces -> matrix
template <typename MeshT>
Rcpp::IntegerMatrix get_faces2(const MeshT &mesh, const std::size_t nSides) {
  const std::size_t nFaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix face_mat(nSides, nFaces);
  std::size_t i = 0;
  for(typename MeshT::Face_index fd : faces(mesh)) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      face_mat(Rcpp::_, i++) = col_i;
  }
  return face_mat;
}

template Rcpp::IntegerMatrix get_faces2<Mesh3>(const Mesh3&,   const std::size_t);
template Rcpp::IntegerMatrix get_faces2<EMesh3>(const EMesh3&, const std::size_t);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for possibly heterogeneous faces -> list
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh1(const MeshT &mesh, const bool normals) {
    using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
    Rcpp::DataFrame     Edges    = get_edges<KernelT, MeshT, PointT>(mesh);
    Rcpp::NumericMatrix Vertices = get_vertices<KernelT, MeshT, PointT>(mesh);
    Rcpp::List          Faces    = get_faces1<MeshT>(mesh);
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                        Rcpp::Named("edges")    = Edges,
                                        Rcpp::Named("faces")    = Faces);

    // check whether there are vertex normals
    std::optional<vertex_normals_map> vnormmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

    // compute new vertex normals?
    if(normals) {
        std::vector<VectorT> vnormals = compute_vnormals<KernelT, MeshT, VectorT>(mesh);
        out["normals"] = points3_to_matrix<KernelT, VectorT>(vnormals);
    } else if(vnormmap_.has_value()) {
        std::optional<std::vector<VectorT>> vnormals = get_vnormals<KernelT, MeshT, VectorT>(mesh);
        if(vnormals.has_value()) {
            out["normals"] = points3_to_matrix<KernelT, VectorT>(vnormals.value());;
        }
    }
    return out;
}

template Rcpp::List make_rmesh1<K,  Mesh3,  Point3,  Vector3>(const  Mesh3&,  const bool);
template Rcpp::List make_rmesh1<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// for homogeneous faces -> matrix
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List make_rmesh2(const MeshT &mesh, const bool normals, const std::size_t nSides) {
    using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
    Rcpp::DataFrame     Edges    = get_edges<KernelT, MeshT, PointT>(mesh);
    Rcpp::NumericMatrix Vertices = get_vertices<KernelT, MeshT, PointT>(mesh);
    Rcpp::IntegerMatrix Faces    = get_faces2<MeshT>(mesh, nSides);
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                        Rcpp::Named("edges") = Edges,
                                        Rcpp::Named("faces") = Faces);

    // check whether there are vertex normals
    std::optional<vertex_normals_map> vnormmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

    // compute new vertex normals?
    if(normals) {
        std::vector<VectorT> vnormals = compute_vnormals<KernelT, MeshT, VectorT>(mesh);
        out["normals"] = points3_to_matrix<KernelT, VectorT>(vnormals);
    } else if(vnormmap_.has_value()) {
        std::optional<std::vector<VectorT>> vnormals = get_vnormals<KernelT, MeshT, VectorT>(mesh);
        if(vnormals.has_value()) {
            out["normals"] = points3_to_matrix<KernelT, VectorT>(vnormals.value());;
        }
    }
    return out;
}

template Rcpp::List make_rmesh2<K,  Mesh3,  Point3,  Vector3>(const  Mesh3&,  const bool, const std::size_t);
template Rcpp::List make_rmesh2<EK, EMesh3, EPoint3, EVector3>(const EMesh3&, const bool, const std::size_t);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// PMP::triangulate_faces() modifies -> mesh not const
template <typename KernelT, typename MeshT, typename PointT, typename VectorT>
Rcpp::List get_rmesh(MeshT &mesh, const bool triangulate, const bool normals) {
  bool is_triangle = CGAL::is_triangle_mesh(mesh);
  if(triangulate && !is_triangle) {
    is_triangle = PMP::triangulate_faces(mesh);
  }
  Rcpp::List rmesh;
  if(is_triangle) {
    rmesh = make_rmesh2<KernelT, MeshT, PointT, VectorT>(mesh, normals, 3);
  } else if(CGAL::is_quad_mesh(mesh)) {
    rmesh = make_rmesh2<KernelT, MeshT, PointT, VectorT>(mesh, normals, 4);
  } else {
    rmesh = make_rmesh1<KernelT, MeshT, PointT, VectorT>(mesh, normals);
  }

  return rmesh;
}

template Rcpp::List get_rmesh<K,  Mesh3,  Point3,  Vector3>( Mesh3&,  const bool, const bool);
template Rcpp::List get_rmesh<EK, EMesh3, EPoint3, EVector3>(EMesh3&, const bool, const bool);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// CAVE: pass by reference, modifies mesh
template <typename MeshT, typename VectorT>
void remove_properties(MeshT &mesh, const std::vector<std::string> &props) {
  using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
  using face_descriptor    = typename boost::graph_traits<MeshT>::face_descriptor;
  using vertex_colors_map  = typename MeshT::template Property_map<vertex_descriptor, CGAL::IO::Color>;
  using face_colors_map    = typename MeshT::template Property_map<face_descriptor,   CGAL::IO::Color>;
  using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
  using vertex_scalars_map = typename MeshT::template Property_map<vertex_descriptor, double>;
  using face_scalars_map   = typename MeshT::template Property_map<face_descriptor,   double>;
  if((props.size() == 1) && (props[0] == "all:all")) {
      mesh.remove_all_property_maps();
      return;
  }
  for(std::size_t i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "v:color") {
      std::optional<vertex_colors_map> pmap_ =
        mesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");
      if(pmap_.has_value()) {
        mesh.remove_property_map(pmap_.value());
      }
    } else if(prop == "f:color") {
      std::optional<face_colors_map> pmap_ =
        mesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");
      if(pmap_.has_value()) {
        mesh.remove_property_map(pmap_.value());
      }
    } else if(prop == "v:normal") {
      std::optional<vertex_normals_map> pmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

      if(pmap_.has_value()) {
        mesh.remove_property_map(pmap_.value());
      }
    } else if(prop == "v:scalar") {
      std::optional<vertex_scalars_map> pmap_ =
        mesh.template property_map<vertex_descriptor, double>("v:scalar");
      if(pmap_.has_value()) {
        mesh.remove_property_map(pmap_.value());
      }
    } else if(prop == "f:scalar") {
      std::optional<face_scalars_map> pmap_ =
        mesh.template property_map<face_descriptor, double>("f:scalar");
      if(pmap_.has_value()) {
        mesh.remove_property_map(pmap_.value());
      }
    }
  }
}

template void remove_properties<Mesh3,  Vector3>(Mesh3&,   const std::vector<std::string>&);
template void remove_properties<EMesh3, EVector3>(EMesh3&, const std::vector<std::string>&);
