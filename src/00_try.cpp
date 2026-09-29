#ifndef _CGALMESHHEADER_
#include "SurfaceMesh.h"
#endif

// ----------------------------------------------------------------------- //
// return existing v:normal property map as R matrix
template <typename KernelT, typename MeshT, typename VectorT>
std::optional<std::vector<VectorT>> get_vnormals(const MeshT &mesh) {
    using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
    std::optional<std::vector<VectorT>> vnormals;
    std::optional<vertex_normals_map> vnormmap_ =
        mesh.template property_map<vertex_descriptor, VectorT>("v:normal");

    if(vnormmap_.has_value()) {
        std::vector<VectorT> vnv;
        vnv.reserve(mesh.number_of_vertices());
        vertex_normals_map vnormmap = vnormmap_.value();
        for(vertex_descriptor vd : vertices(mesh)) {
            vnv.emplace_back(vnormmap[vd]);
        }
        vnormals = std::move(vnv);
    }
    return vnormals;
}

template std::optional<std::vector<Vector3>>  get_vnormals<K,  Mesh3,  Vector3>(const  Mesh3&);
template std::optional<std::vector<EVector3>> get_vnormals<EK, EMesh3, EVector3>(const EMesh3&);

// ----------------------------------------------------------------------- //
// set vertex normal property map to given input
// mesh is changed -> not const
template <typename MeshT, typename VectorT>
void set_vnormals(
    MeshT &mesh, const std::vector<VectorT> &vnormals) {
    using vertex_descriptor  = typename boost::graph_traits<MeshT>::vertex_descriptor;
    using vertex_normals_map = typename MeshT::template Property_map<vertex_descriptor, VectorT>;
    if(vnormals.size() != mesh.number_of_vertices()) {
        Rcpp::stop("The number of normals does not match the number of vertices.");
    }
    remove_properties<MeshT, VectorT>(mesh, {"v:normal"});
    vertex_normals_map vnormmap =
        mesh.template add_property_map<vertex_descriptor, VectorT>(
            "v:normal", CGAL::NULL_VECTOR).first;
    std::size_t i = 0;
    for(vertex_descriptor vd : vertices(mesh)) {
        vnormmap[vd] = vnormals[i];
        i++;
    }
}

template void set_vnormals<Mesh3,  Vector3>( Mesh3&,  const std::vector<Vector3>&);
template void set_vnormals<EMesh3, EVector3>(EMesh3&, const std::vector<EVector3>&);

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::NumericMatrix> getVertexNormals_cpp(const Rcpp::List rmesh) {
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        true,        // triangulate - must be triangle
        false,       // repair_soup
        false);      // verbose
    const std::optional<std::vector<Vector3>> vnormals = get_vnormals<K, Mesh3, Vector3>(mesh);
    if(vnormals.has_value()) {
        const Rcpp::NumericMatrix normals_mat = points3_to_matrix<K, Vector3>(vnormals.value());
        return Rcpp::transpose(normals_mat);
    } else {
        return R_NilValue;
    }
}

// ----------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List setVertexNormals_cpp(const Rcpp::List rmesh,
                                const Rcpp::NumericMatrix rnormals) {
    using vertex_normals_map = Mesh3::Property_map<vrtx_dscrptr, Vector3>;
    Mesh3 mesh = make_surf_mesh_valid<Mesh3, Point3>(
        rmesh,
        false,       // soup
        false,       // triangulate
        false,       // repair_soup
        false);      // verbose
    std::vector<Vector3> normals = matrix_to_points3<Vector3>(rnormals);
    set_vnormals<Mesh3, Vector3>(mesh, normals);
    return get_rmesh<K, Mesh3, Point3, Vector3>(mesh, false, false);
}
