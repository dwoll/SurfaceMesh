#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

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

template std::optional<Rcpp::IntegerMatrix> getFColors<K,  Mesh3,  Vector3>(const Mesh3&);
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
Rcpp::Nullable<Rcpp::NumericMatrix> getVertexNormals_cpp(const Rcpp::List rmesh) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    const std::optional<Rcpp::NumericMatrix> normals_mat = getVNormals<K, Mesh3, Vector3>(mesh);
    if(normals_mat.has_value()) {
        return Rcpp::transpose(normals_mat.value());
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
    remove_properties<Mesh3, Vector3>(mesh, {"f:color"});
    face_colors_map fcolmap =
        mesh.add_property_map<fc_dscrptr, CGAL::IO::Color>("f:color").first;
        // ("f:color", CGAL::IO::white()).first for default
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
// [[Rcpp::export]]
Rcpp::List setVertexNormals(const Rcpp::List rmesh,
                            const Rcpp::NumericMatrix normals) {
    using vertex_normals_map = Mesh3::Property_map<vrtx_dscrptr, Vector3>;
    rmessage("Before make_surf_mesh_valid");
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        false,       // triangulate
        false,       // repair_soup
        false);      // verbose
    rmessage("After make_surf_mesh_valid");
    if(normals.ncol() != mesh.number_of_vertices()) {
        Rcpp::stop("The number of normals does not match the number of vertices.");
    }
    rmessage("Before remove_properties");
    remove_properties<Mesh3, Vector3>(mesh, {"v:normal"});
    rmessage("After remove_properties");
    vertex_normals_map vnormmap =
        mesh.add_property_map<vrtx_dscrptr, Vector3>(
            "v:normal", CGAL::NULL_VECTOR).first;
    std::size_t i = 0;
    rmessage("Before for() loop");
    for(Mesh3::Vertex_index vi : mesh.vertices()) {
        const Rcpp::NumericVector v = normals(Rcpp::_, i);
        vnormmap[vi] = Vector3(v(0), v(1), v(2));
        i++;
    }
    rmessage("After for() loop");
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, false);
}
