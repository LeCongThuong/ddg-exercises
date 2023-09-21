// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <Eigen/Sparse>


using namespace geometrycentral;
using namespace geometrycentral::surface;


typedef Eigen::Triplet<size_t> T;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    std::vector<T> tripletList;
    size_t nVertices = mesh->nVertices();
    size_t nEdges = mesh->nEdges();
    
    for(Edge e: mesh->edges()){
        size_t idx_row = e.getIndex();
        size_t idx_col_1 = e.firstVertex().getIndex();
        size_t idx_col_2 = e.secondVertex().getIndex();
        tripletList.push_back(T(idx_row, idx_col_1, 1));
        tripletList.push_back(T(idx_row, idx_col_2, 1)); 
    }

    SparseMatrix<size_t> a_0(nEdges,nVertices);
    a_0.setFromTriplets(tripletList.begin(),  tripletList.end());
    return a_0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    std::vector<T> tripletList;
    size_t nFaces = mesh->nFaces();
    size_t nEdges = mesh->nEdges();

    for(Face f: mesh->faces()){
        for(Vertex v : f.adjacentVertices()){
            size_t row_idx = f.getIndex();
            size_t idx_col = v.getIndex();
            tripletList.push_back(T(row_idx, idx_col, 1));
        }
    }
    SparseMatrix<size_t> a_1(nFaces, nEdges);
    a_1.setFromTriplets(tripletList.begin(),  tripletList.end());
    return a_1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    size_t num_vertices = subset.vertices.size();

    Vector<size_t> subset_simplices_vertices;
    size_t count = 0;
    for(auto v: subset.vertices){
        subset_simplices_vertices(count, 0) = v;
        count += 1;
    }
    return subset_simplices_vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    size_t num_edges = subset.edges.size();

    Vector<size_t> subset_simplices_edges;
    size_t count = 0;
    for(auto e: subset.edges){
        subset_simplices_edges(count, 0) = e;
        count += 1;
    }
    return subset_simplices_edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    size_t num_edges = subset.faces.size();

    Vector<size_t> subset_simplices_faces;
    size_t count = 0;
    for(auto e: subset.faces){
        subset_simplices_faces(count, 0) = e;
        count += 1;
    }
    return subset_simplices_faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> face_simplices;
    std::set<size_t> edge_simplices;

    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    for(size_t v: vertices){
        for(Edge e: mesh->edges()){
            size_t idx_row = e.getIndex();
            size_t idx_col_1 = e.firstVertex().getIndex();
            size_t idx_col_2 = e.secondVertex().getIndex();
            if(v == idx_col_1 || v == idx_col_2){
                edge_simplices.insert(idx_row);
            }
        }
    }


    for(size_t e: edges){
        edge_simplices.insert(e);
    }

    for(size_t e: edge_simplices){
        for(Face f: mesh->faces()){
            for(Vertex v : f.adjacentVertices()){
                size_t row_idx = f.getIndex();
                size_t idx_col = v.getIndex();
                if(e == idx_col){
                    face_simplices.insert(row_idx);
                    break;
                }
            }
        }
    }
 

    for(size_t f: faces){
        face_simplices.insert(f);
    }

    return MeshSubset(vertices, edge_simplices, face_simplices);
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    //Add all vertices, edges and faces of subset into closure. Beside that, it also add the edges, vertices to make sure Closure is the SubComplex.
    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    std::set<size_t> face_simplices;
    std::set<size_t> edge_simplices;
    std::set<size_t> vertice_simplices;

    std::vector<Edge> mesh_edge; 
    for (Edge e : mesh->edges()) {
        mesh_edge.push_back(e);
    }

    std::vector<Vertex> mesh_vertex; 
    for (Vertex v: mesh->vertices()) {
        mesh_vertex.push_back(v);
    }

    std::vector<Face> mesh_face; 
    for (Face f: mesh->faces()) {
        mesh_face.push_back(f);
    }

    for(size_t f_index: faces){
        face_simplices.insert(f_index);
        Face f = mesh_face[f_index];
        for(Edge e : f.adjacentEdges()){
            size_t idx_col = e.getIndex();
            edge_simplices.insert(idx_col);
              
        }
    }

    for(size_t e_index: edges){
        edge_simplices.insert(e_index);
    }

    for(size_t e_index: edge_simplices){
        Edge e = mesh_edge[e_index];
        size_t idx_col_1 = e.firstVertex().getIndex();
        size_t idx_col_2 = e.secondVertex().getIndex();
        vertice_simplices.insert(idx_col_1);
        vertice_simplices.insert(idx_col_2);
    }

    for(size_t v_index: vertices){
        vertice_simplices.insert(v_index);
    }

    return MeshSubset(vertices, edge_simplices, face_simplices);

    
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset cl_st = closure(star(subset));
    MeshSubset st_cl = star(closure(subset));

    std::set<size_t> cl_st_vertices = cl_st.vertices;
    std::set<size_t> cl_st_edges = cl_st.edges;
    std::set<size_t> cl_st_faces = cl_st.faces;

    std::set<size_t> st_cl_vertices = st_cl.vertices;
    std::set<size_t> st_cl_edges = st_cl.edges;
    std::set<size_t> st_cl_faces = st_cl.faces;

    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;
    
    std::set_difference(cl_st_vertices.begin(), cl_st_vertices.end(), st_cl_vertices.begin(), st_cl_vertices.end(), std::inserter(vertices, vertices.begin()));
    std::set_difference(cl_st_edges.begin(), cl_st_edges.end(), st_cl_edges.begin(), st_cl_edges.end(), std::inserter(edges, edges.begin()));
    std::set_difference(cl_st_faces.begin(), cl_st_faces.end(), st_cl_faces.begin(), st_cl_faces.end(), std::inserter(faces, faces.begin()));

    return MeshSubset(vertices, edges, faces);
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */

bool isSubset(std::set<size_t> s1, std::set<size_t> s2){
    std::set<size_t> res;
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(res, res.begin()));
    if(res.size() == 0){
        return true;
    }
    return false;
}


bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // Iterate the each faces, list the all edges and vertices -> check the the set is the subset of original edges and vertices 

    // Iterate the each edges, list all vertices -> check the set is the subset of the original of vertices 
    std::set<size_t> face_simplices;
    std::set<size_t> edge_simplices;
    std::set<size_t> vertice_simplices;

    std::vector<Edge> mesh_edge; 
    for (Edge e : mesh->edges()) {
        mesh_edge.push_back(e);
    }

    std::vector<Vertex> mesh_vertex; 
    for (Vertex v: mesh->vertices()) {
        mesh_vertex.push_back(v);
    }

    std::vector<Face> mesh_face; 
    for (Face f: mesh->faces()) {
        mesh_face.push_back(f);
    }

    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    for(size_t f_index: faces){
        Face f = mesh_face[f_index];
        for(Edge e : f.adjacentEdges()){
            size_t idx_col = e.getIndex();
            edge_simplices.insert(idx_col);
            size_t idx_col_1 = e.firstVertex().getIndex();
            size_t idx_col_2 = e.secondVertex().getIndex();
            vertice_simplices.insert(idx_col_1);
            vertice_simplices.insert(idx_col_2); 
        }

    }
    if(!isSubset(edge_simplices, edges) || !isSubset(vertice_simplices, vertices)){
        return false;
    }

    for(size_t e_index: edges){
        Edge e = mesh_edge[e_index];
        size_t idx_col_1 = e.firstVertex().getIndex();
        size_t idx_col_2 = e.secondVertex().getIndex();
        vertice_simplices.insert(idx_col_1);
        vertice_simplices.insert(idx_col_2);
    }

     if(!isSubset(vertice_simplices, vertices)){
        return false;
    }

    return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}