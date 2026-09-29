// ----------------------------------------------------------------------- //
// CAVE
// the code in this file is not tested
// it may or may not work
// ----------------------------------------------------------------------- //

template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::IntegerMatrix> getFColors(const MeshT&);

template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::IntegerMatrix> getVColors(const MeshT&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// return existing f:color property map as R matrix
template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::IntegerMatrix> getFColors(const MeshT &mesh) {
    using face_descriptor = typename boost::graph_traits<MeshT>::face_descriptor;
    using face_colors_map = typename MeshT::template Property_map<face_descriptor, CGAL::IO::Color>;
    std::optional<Rcpp::IntegerMatrix> colors_mat;
    std::optional<face_colors_map> pmap_ =
        mesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");

    if(pmap_.has_value()) {
        rmessage("pmap_ has value.");
        std::size_t i = 0;
        Rcpp::IntegerMatrix cm(4, mesh.number_of_faces());
        face_colors_map fcolmap = pmap_.value();
        for(face_descriptor fd : faces(mesh)) {
            Rcpp::IntegerVector col_i(4);
            CGAL::IO::Color color = fcolmap[fd];   // CGAL::SM_Face_index(i)
            col_i(0) = color.red();
            col_i(1) = color.green();
            col_i(2) = color.blue();
            col_i(3) = color.alpha();
            cm(Rcpp::_, i) = col_i;
            i++;
        }
        colors_mat = std::move(cm);
    }
    return colors_mat;
}

template std::optional<Rcpp::IntegerMatrix> getFColors<K,  Mesh3,  Vector3>(const  Mesh3&);
template std::optional<Rcpp::IntegerMatrix> getFColors<EK, EMesh3, EVector3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// return existing v:color property map as R matrix
template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::IntegerMatrix> getVColors(const MeshT &mesh) {
    using vertex_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_colors_map = typename MeshT::template Property_map<vertex_descriptor, CGAL::IO::Color>;
    std::optional<Rcpp::IntegerMatrix> colors_mat;
    std::optional<vertex_colors_map> pmap_ =
        mesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");

    if(pmap_.has_value()) {
        std::size_t i = 0;
        Rcpp::IntegerMatrix cm(4, mesh.number_of_vertices());
        vertex_colors_map vcolmap = pmap_.value();
        for(vertex_descriptor vd : vertices(mesh)) {
            Rcpp::IntegerVector col_i(4);
            CGAL::IO::Color color = vcolmap[vd];    // CGAL::SM_Vertex_index(i)
            col_i(0) = color.red();
            col_i(1) = color.green();
            col_i(2) = color.blue();
            col_i(3) = color.alpha();
            cm(Rcpp::_, i) = col_i;
            i++;
        }
        colors_mat = std::move(cm);
    }
    return colors_mat;
}

template std::optional<Rcpp::IntegerMatrix> getVColors<K,  Mesh3,  Vector3>(const Mesh3&);
template std::optional<Rcpp::IntegerMatrix> getVColors<EK, EMesh3, EVector3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// return existing v:normal property map as R matrix
template <typename KernelT, typename MeshT, typename VectorT>
std::optional<Rcpp::NumericMatrix> getVNormals(const MeshT &mesh) {
    using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
    std::optional<Rcpp::NumericMatrix> normals_mat;
    std::optional<vertex_normals_map> pmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

    if(pmap_.has_value()) {
        std::size_t i = 0;
        Rcpp::NumericMatrix nm(3, mesh.number_of_vertices());
        vertex_normals_map vnormmap = pmap_.value();
        for(vertex_descriptor vd : vertices(mesh)) {
            Rcpp::NumericVector col_i(3);
            VectorT normal = vnormmap[vd];    // CGAL::SM_Vertex_index(i)
            col_i(0) = CGAL::to_double<typename KernelT::FT>(normal.x());
            col_i(1) = CGAL::to_double<typename KernelT::FT>(normal.y());
            col_i(2) = CGAL::to_double<typename KernelT::FT>(normal.z());
            nm(Rcpp::_, i) = col_i;
            i++;
        }
        normals_mat = std::move(nm);
    }
    return normals_mat;
}

template std::optional<Rcpp::NumericMatrix> getVNormals<K,  Mesh3,  Vector3>(const Mesh3&);
template std::optional<Rcpp::NumericMatrix> getVNormals<EK, EMesh3, EVector3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::IntegerMatrix> getFaceColors_cpp(const Rcpp::List rmesh) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    const std::optional<Rcpp::IntegerMatrix> colors_mat = getFColors<K, Mesh3, Vector3>(mesh);
    if(colors_mat.has_value()) {
        return Rcpp::transpose(colors_mat.value());
    } else {
        return R_NilValue;
    }
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::IntegerMatrix> getVertexColors_cpp(const Rcpp::List rmesh) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    const std::optional<Rcpp::IntegerMatrix> colors_mat = getVColors<K, Mesh3, Vector3>(mesh);
    if(colors_mat.has_value()) {
        return Rcpp::transpose(colors_mat.value());
    } else {
        return R_NilValue;
    }
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List setFaceColors_cpp(const Rcpp::List rmesh,
                             const Rcpp::IntegerMatrix colors) {
    using face_colors_map = Mesh3::Property_map<fc_dscrptr, CGAL::IO::Color>;
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        false,       // triangulate
        false,       // repair_soup
        false);      // verbose
    if((colors.ncol() != 1) && (colors.ncol() != mesh.number_of_faces())) {
        Rcpp::stop("The number of colors does not match the number of faces.");
    }
    if(colors.nrow() != 4) {
        Rcpp::stop("colors must have 4 channels for R, G, B, A.");
    }
    remove_properties<Mesh3, Vector3>(mesh, {"f:color"});
    face_colors_map fcolmap =
        mesh.add_property_map<fc_dscrptr, CGAL::IO::Color>("f:color", CGAL::IO::white()).first;
    if(colors.ncol() == 1) {
        for(Mesh3::Face_index fi : mesh.faces()) {
            const Rcpp::IntegerVector v = colors(Rcpp::_, 0);
            fcolmap[fi] = CGAL::IO::Color(v(0), v(1), v(2), v(3));
        }
    } else {
        std::size_t i = 0;
        for(Mesh3::Face_index fi : mesh.faces()) {
            const Rcpp::IntegerVector v = colors(Rcpp::_, i);
            fcolmap[fi] = CGAL::IO::Color(v(0), v(1), v(2), v(3));
            i++;
        }
    }
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, false);
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List setVertexColors_cpp(const Rcpp::List rmesh,
                               const Rcpp::IntegerMatrix colors) {
    using vertex_colors_map = Mesh3::Property_map<vrtx_dscrptr, CGAL::IO::Color>;
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        false,       // triangulate
        false,       // repair_soup
        false);      // verbose
    if(colors.ncol() != mesh.number_of_vertices()) {
        Rcpp::stop("The number of colors does not match the number of vertices.");
    }
    remove_properties<Mesh3, Vector3>(mesh, {"v:color"});
    vertex_colors_map vcolmap =
        mesh.add_property_map<vrtx_dscrptr, CGAL::IO::Color>("v:color").first;
    std::size_t i = 0;
    for(Mesh3::Vertex_index vi : mesh.vertices()) {
        const Rcpp::IntegerVector v = colors(Rcpp::_, i);
        vcolmap[vi] = CGAL::IO::Color(v(0), v(1), v(2), v(3));
        i++;
    }
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, false);
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

// in make_rmesh*()
using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
using face_descriptor    = typename boost::graph_traits<MeshT>::face_descriptor;
using vertex_colors_map  = typename MeshT::template Property_map<vertex_descriptor, CGAL::IO::Color>;
using face_colors_map    = typename MeshT::template Property_map<face_descriptor,   CGAL::IO::Color>;
using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;

std::optional<face_colors_map> fcolmap_ =
  mesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");

std::optional<vertex_colors_map> vcolmap_ =
  mesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");

std::optional<vertex_normals_map> vnormmap_ =
  mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

if () {
    ...
} else if(vnormmap_.has_value()) {
    std::optional<Rcpp::NumericMatrix> normals_mat = getVNormals<KernelT, MeshT, VectorT>(mesh);
    if(normals_mat.has_value()) {
        out["normals"] = normals_mat.value();
    }
}

// vertex colors or face colors? either - or
if(fcolmap_.has_value() || vcolmap_.has_value()) {
    if(vcolmap_.has_value()) {
        std::optional<Rcpp::IntegerMatrix> colors_mat = getVColors<KernelT, MeshT, VectorT>(mesh);
        if(colors_mat.has_value()) {
            out["colors"] = colors_mat.value();
        }
    } else {
        std::optional<Rcpp::IntegerMatrix> colors_mat = getFColors<KernelT, MeshT, VectorT>(mesh);
        if(colors_mat.has_value()) {
            out["colors"] = colors_mat.value();
        }
    }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

typedef boost::graph_traits<Mesh3>::vertex_descriptor                vrtx_dscrptr;
typedef Mesh3::Property_map<vrtx_dscrptr, Rcpp::NumericVector>       nrmls_map_r;

typedef boost::graph_traits<EMesh3>::vertex_descriptor               vrtx_descriptor;
typedef EMesh3::Property_map<vrtx_descriptor, Rcpp::NumericVector>   normals_map_r;

typedef boost::graph_traits<EMesh3>::edge_descriptor                 dg_descriptor;
typedef boost::graph_traits<EMesh3>::halfedge_descriptor             hlfdg_descriptor;

// EPoint3 with normal EVector3
typedef std::pair<EPoint3, EVector3>                                 EP3EV3;
typedef boost::graph_traits<EMesh3>::face_descriptor                 fc_descriptor;

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

template <typename MeshT, typename PointT>
bool is_small_hole(typename boost::graph_traits<MeshT>::halfedge_descriptor h,
                   const MeshT &mesh,
                   const double max_hole_diam,
                   const int max_num_hole_edges) {
  using halfedge_descriptor = typename boost::graph_traits<MeshT>::halfedge_descriptor;
  int num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh)) {
    const PointT& p = mesh.point(target(hc, mesh));

    hole_bbox += p.bbox();
    ++num_hole_edges;

    // exit early, to avoid unnecessary traversal of large holes
    if (num_hole_edges > max_num_hole_edges) return false;
    if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
    if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
    if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
  }

  return true;
}

template bool is_small_hole<Mesh3,  Point3>(typename  boost::graph_traits<Mesh3>::halfedge_descriptor,  const Mesh3&,  const double, const int);
template bool is_small_hole<EMesh3, EPoint3>(typename boost::graph_traits<EMesh3>::halfedge_descriptor, const EMesh3&, const double, const int);

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

Rcpp::NumericVector defaultNormal() {
  Rcpp::NumericVector def =
    {
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na(),
      Rcpp::NumericVector::get_na()
    };
  return def;
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// compatibility wrapper for CGAL property_map(std::string) API changes:
// older returned std::pair<Property_map, bool>
// newer returns  std::optional<Property_map>
// property_map_pair returns a std::pair<Property_map, bool> in both cases

template <typename KeyT, typename T, typename MeshT>
std::pair<typename MeshT::template Property_map<KeyT,T>, bool>
property_map_pair(MeshT &mesh, const std::string name) {
  using RetType = decltype(mesh.template property_map<KeyT,T>(name));
  using Pmap    = typename MeshT::template Property_map<KeyT,T>;
  if constexpr (std::is_same_v<RetType, std::pair<Pmap, bool>>) {
    return mesh.template property_map<KeyT,T>(name);
  } else {
    auto opt = mesh.template property_map<KeyT,T>(name);
    if(opt) {
        return std::make_pair(*opt, true);
    }
    return std::make_pair(Pmap(), false);
  }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

(const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_)
using norm_map_r   = typename MeshT::template Property_map<v_descriptor, Rcpp::NumericVector>;
using vertex_descriptor = typename boost::graph_traits<MeshT>::vertex_descriptor;
if(normals_.isNotNull()) {
  Rcpp::NumericMatrix normals_mat(normals_);
  const unsigned int nNormals = static_cast<unsigned int>(normals_mat.ncol());
  if(mesh.number_of_vertices() != nNormals) {
    Rcpp::stop(
      "The number of normals does not match the number of vertices.");
  }
  Rcpp::NumericVector def = defaultNormal();
  norm_map_r normalsmap =
    mesh.template add_property_map<v_descriptor, Rcpp::NumericVector>(
      "v:normal", def).first;
  for(std::size_t j = 0; j < nNormals; j++) {
    Rcpp::NumericVector normal = normals_mat(Rcpp::_, j);
    normalsmap[CGAL::SM_Vertex_index(j)] = normal;
  }
}
