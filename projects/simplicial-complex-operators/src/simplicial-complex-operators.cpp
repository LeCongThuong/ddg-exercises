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
        for(Edge e: f.adjacentEdges()){
            size_t row_idx = f.getIndex();
            size_t idx_col = e.getIndex();
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
    // Init zero encoding vertex vector
    Vector<size_t> encoding_vertex_vec;
    for(size_t i=0; i < mesh->nVertices(); i++){
        encoding_vertex_vec(i) = 0;
    }

    for(size_t v_index: subset.vertices){
        encoding_vertex_vec(v_index) = 1;
    }
    return encoding_vertex_vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> encoding_edge_vec;
    for(size_t i=0; i < mesh->nEdges(); i++){
        encoding_edge_vec(i) = 0;
    }

    for(size_t e_index: subset.edges){
        encoding_edge_vec(e_index) = 1;
    }
    return encoding_edge_vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> encoding_face_vec;
    for(size_t i=0; i < mesh->nFaces(); i++){
        encoding_face_vec(i) = 0;
    }

    for(size_t f_index: subset.faces){
        encoding_face_vec(f_index) = 1;
    }
    return encoding_face_vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    /*
        1. Init star with all faces, edges, vertices from subset
        3. Iterate each complex (vertices, edges) of the subset, check the complex belong what complices of mesh based on these adj matrices 
    */

    MeshSubset st = subset.deepCopy();

    for (size_t v_index: st.vertices){
        for (SparseMatrix<size_t>::InnerIterator it(A0, v_index); it; ++it)
        {
            st.addEdge(it.row());
        }
    }
    
    for(size_t e_index: st.edges){
        for (SparseMatrix<size_t>::InnerIterator it(A1, e_index); it; ++it){
            st.addFace(it.row());
        }
    }
    return st;

}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    MeshSubset lk = subset.deepCopy();
    std::set<size_t> faces = subset.faces;
    std::set<size_t> edges = subset.edges;

    for (int k=0; k<A1.outerSize(); ++k){
        for (SparseMatrix<size_t>::InnerIterator it(A1,k); it; ++it){
            size_t r_index = it.row();   // row index
            size_t c_index = it.col();   // col index (here it is equal to k)
            if(faces.find(r_index) != faces.end()){
                lk.addEdge(c_index);
            }
        }
    }

    for (int k=0; k<A0.outerSize(); ++k){
        for (SparseMatrix<size_t>::InnerIterator it(A0,k); it; ++it){
            size_t r_index = it.row();   // row index
            size_t c_index = it.col();   // col index (here it is equal to k)
            if(edges.find(r_index) != edges.end()){
                lk.addVertex(c_index);
            }
        }
    }
    return lk;
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

    MeshSubset lk = st_cl.deepCopy();
    lk.deleteSubset(cl_st);

    return lk;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */


bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    MeshSubset cl = closure(subset);
    if(subset.vertices == cl.vertices && subset.edges == cl.edges){
        return true;
    }
    else{
        return false;
    }
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */

int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    
    if(isComplex(subset)){
        if(A1.sum() > 0){
            for (int k=0; k<A1.outerSize(); ++k){
                bool flag = false;
                for (SparseMatrix<size_t>::InnerIterator it(A1,k); it; ++it){
                    if(it.value() != 0){
                        flag = true;
                        break;
                    }
                }
                if(!flag){
                    return -1;
                } 
            }
        }
       

        for (int k=0; k<A0.outerSize(); ++k){
            bool flag = false;
            for (SparseMatrix<size_t>::InnerIterator it(A0,k); it; ++it){
                if(it.value() != 0){
                    flag = true;
                    break;
                }
            }
            if(!flag){
                return -1;
            }
        }
        return 2;
        
    }else{
        return -1; // placeholder
    }
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