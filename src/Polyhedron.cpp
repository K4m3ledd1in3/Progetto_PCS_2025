#include "Eigen/Eigen"
#include "Polygon.hpp"
#include <iostream>
#include "Polyhedron.hpp"
#include <limits>
#include <vector>
#include <set>
#include <fstream>
using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;
namespace PolyhedronLibrary{

    //______________________
    


   
    //__________________
    void _Polyhedron::WriteTXT() const
    {
        // Scrittura di Cell0Ds.txt: informazioni sui vertici
        std::ofstream ofile0("Cell0Ds.txt");
        if (ofile0.is_open()) {
            for (size_t i = 0; i < vertices.size(); ++i) {
                ofile0 << vertices[i].id << ", " << vertices[i].x << ", " << vertices[i].y << ", " << vertices[i].z << ", " << vertices[i].ShortPath << "\n";
            }
            ofile0.close();
        }

        // Scrittura di Cell1Ds.txt: informazioni sui lati
        std::ofstream ofile1("Cell1Ds.txt");
        if (ofile1.is_open()) {
            for (size_t i = 0; i < edges.size(); ++i) {
                ofile1 << edges[i].id << ", " << edges[i].origin.id << ", " << edges[i].end.id << ", " << edges[i].ShortPath << ", " << edges[i].length() << "\n";
            }
            ofile1.close();
        }

        // Scrittura di Cell2Ds.txt: informazioni sulle facce
        std::ofstream ofile2("Cell2Ds.txt");
        if (ofile2.is_open()) {
            for (size_t i = 0; i < faces.size(); ++i) {
                ofile2 << faces[i].id << ", " << faces[i].vertices.size() << ", " << faces[i].edges.size() << ", ";
                for (size_t j = 0; j < faces[i].vertices.size(); ++j) {
                    ofile2 << faces[i].vertices[j].id << " ";
                }
                ofile2 << ", ";
                for (size_t j = 0; j < faces[i].edges.size(); ++j) {
                    ofile2 << faces[i].edges[j].id << " ";
                }
                ofile2 << "\n";
            }
            ofile2.close();
        }

        // Scrittura di Cell3Ds.txt: informazioni generali del poliedro
        std::ofstream ofile3("Cell3Ds.txt");
        if (ofile3.is_open()) {
            ofile3 << "0, " << vertices.size() << ", " << edges.size() << ", " << faces.size() << ", ";
            for (size_t i = 0; i < vertices.size(); ++i) {
                ofile3 << vertices[i].id << " ";
            }
            ofile3 << ", ";
            for (size_t i = 0; i < edges.size(); ++i) {
                ofile3 << edges[i].id << " ";
            }
            ofile3 << ", ";
            for (size_t i = 0; i < faces.size(); ++i) {
                ofile3 << faces[i].id << " ";
            }
            ofile3 << "\n";
            ofile3.close();
        }
}
// Funzione per calcolare il cammino più breve tra due vertici
    std::vector<PolygonalLibrary::vertex> ShortestPath(
    std::vector<PolygonalLibrary::vertex>& vertices,
    std::vector<PolygonalLibrary::Edge>& edges,
    unsigned int id1, unsigned int id2) {
    const int N = vertices.size();
    std::vector<double> dist(N, std::numeric_limits<double>::infinity());
    std::vector<int> prev(N, -1);
    std::vector<bool> visited(N, false);
    dist[id1] = 0.0;
    while (true) {
        int u = -1;
        double min_dist = numeric_limits<double>::infinity();

        // Trova il vertice non visitato con distanza minima
        for (int i = 0; i < N; i++) {
            if (!visited[i] && dist[i] < min_dist) {
                min_dist = dist[i];
                u = i;
            }
        }
        if (u == -1 || u == id2) break; // Nessun vertice raggiungibile o raggiunto il target
        visited[u] = true;
        // Aggiorna le distanze dei vicini tramite gli archi
        for (size_t i = 0; i < edges.size(); i++) {
            unsigned int from = edges[i].origin.id;
            unsigned int to = edges[i].end.id;
            double length = edges[i].length();

            if (from == u && dist[to] > dist[u] + length) {
                dist[to] = dist[u] + length;
                prev[to] = u;
            }
            if (to == u && dist[from] > dist[u] + length) {
                dist[from] = dist[u] + length;
                prev[from] = u;
            }
        }
    }

    // Ricostruzione del cammino
    std::vector<unsigned int> path_ids;
    int at = id2;
    if (prev[at] == -1 && at != (int)id1) {
        std::cout << "Nessun cammino trovato.\n";
        return std::vector<PolygonalLibrary::vertex>();
    }

    while (at != -1) {
        path_ids.push_back(at);
        at = prev[at];
    }

    // Inverti il cammino
    std::vector<unsigned int> correct_path;
    for (int i = path_ids.size() - 1; i >= 0; i--) {
        correct_path.push_back(path_ids[i]);
    }
    // Marcare ShortPath
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i].ShortPath = 0;
    }
    for (size_t i = 0; i < correct_path.size(); i++) {
        vertices[correct_path[i]].ShortPath = 1;
    }

    for (size_t i = 0; i < edges.size(); i++) {
        unsigned int from = edges[i].origin.id;
        unsigned int to = edges[i].end.id;
       for (int i = 0; i < correct_path.size() - 1; ++i) {
        unsigned int from = correct_path[i];
        unsigned int to = correct_path[i + 1];

        for (auto& edge : edges) {
            if ((edge.origin.id == from && edge.end.id == to) ||
                (edge.origin.id == to && edge.end.id == from)) {
                edge.ShortPath = 1;
                break;
            }
        }
}

    }
    // Calcolo lunghezza e numero lati
    double total_length = 0.0;
    int num_edges = correct_path.size() - 1;
    for (int i = 0; i < num_edges; i++) {
        unsigned int from = correct_path[i];
        unsigned int to = correct_path[i+1];
        for (size_t j = 0; j < edges.size(); j++) {
            if ((edges[j].origin.id == from && edges[j].end.id == to) ||
                (edges[j].origin.id == to && edges[j].end.id == from)) {
                total_length += edges[j].length();
                break;
            }
        }
    }
    cout << "Cammino minimo tra " << id1 << " e " << id2 << ": "
              << num_edges << " lati, lunghezza totale " << total_length << "\n";

    // Ritorna il cammino come vector di vertex
    vector<PolygonalLibrary::vertex> path_vertices;
    for (size_t i = 0; i < correct_path.size(); i++) {
        path_vertices.push_back(vertices[correct_path[i]]);
    }

    return path_vertices;
}

	void Icosahedron::display() const{
			for(size_t i=0; i<vertices.size(); i++){
				cout << "(" << this->vertices[i].x << "," << this->vertices[i].y << "," << this->vertices[i].z << ")" << endl;}
		}		
	void Icosahedron::printFaces() {
    		for(auto face: faces){
		printFace(face);
			}	
			}
	void _Polyhedron::printFaces() {
    			for(auto face: faces){
				printFace(face);
				}
				}
	void Dodecahedron::display() const{
			for(size_t i=0; i<vertices.size(); i++){
				cout << "(" << this->vertices[i].x << "," << this->vertices[i].y << "," << this->vertices[i].z << ")" << endl;}
		}	
	void Dodecahedron::printFaces() {
    		for(auto face: faces){
				printFace(face);
			}	
			}			
	void Octahedron::display() const{
			for(size_t i=0; i<vertices.size(); i++){
				cout << "(" << this->vertices[i].x << "," << this->vertices[i].y << "," << this->vertices[i].z << ")" << endl;}
		}
	void Octahedron::printFaces() {
    		for(auto face: faces){
		printFace(face);
		}	
		}
	void Cube::display() const{
			for(size_t i=0; i<vertices.size(); i++){
				cout << "(" << this->vertices[i].x << "," << this->vertices[i].y << "," << this->vertices[i].z << ")" << endl;}
		}
	void Cube::printFaces() {
    for(auto face: faces){
		printFace(face);
	}	
	}
	void Tetrahedron::display() const{
			for(size_t i=0; i<vertices.size(); i++){
				cout << "(" << vertices[i].x << "," << vertices[i].y << "," << vertices[i].z << ")" << endl;}
		}
	void Tetrahedron::printFaces() {
    for(auto face: faces){
		printFace(face);
	}	
	}
    void _Polyhedron::Triangulation() {
    vector<Edge> edges_1;
    vector<Face> faces_1;
    vector<vertex> ver_1;
    if ((b == 0 && c != 0) || (b != 0 && c == 0)) {
        if (p == 3) {
            unsigned int T = b * b + c * b + c * c;

            switch (q) {
                case 3:
                    NumVer = 2 * T + 2; NumEdg = 6 * T; NumFcs = 4 * T;
                    break;
                case 4:
                    NumVer = 4 * T + 2; NumEdg = 12 * T; NumFcs = 8 * T;
                    break;
                case 5:
                    NumVer = 10 * T; NumEdg = 30 * T; NumFcs = 20 * T;
                    break;
            }

            faces_1.reserve(NumFcs);
            edges_1.reserve(NumEdg);
            edges.reserve(NumEdg);
            vertices.reserve(NumVer);
            ver_1.reserve(NumVer);
            unsigned int id_edge = 0;
            unsigned int id_faces = 0;
            vertex v;
            vertex v1;
            map<unsigned int, vector<vertex>> Triangle;
            vector<size_t> index_edges(3);
            size_t cv = 0;
            if(b==0){
                b=c;
            }else{
                c=b;
            }
            unsigned int count=0;
            map<tuple<double, double, double>, unsigned int> vertex_map;
            unsigned int next_id = 0;

            for (size_t i = 0; i < faces.size(); i++) {
                for (size_t k = 0; k < b ; k++) {
                    Triangle[k].resize(k + 2);
                    vertex v0;
                    vertex v_origin;

                    if (faces[i].edges[1].end != faces[i].edges[0].end) {
                        if (faces[i].edges[1].end == faces[i].edges[0].origin) {
                            v = faces[i].edges[0].origin;
                            v0 = vertex(
                                (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (k + 1) / float(b) + faces[i].edges[0].origin.x,
                                (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (k + 1) / double(b) + faces[i].edges[0].origin.y,
                                (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (k + 1) / double(b) + faces[i].edges[0].origin.z,
                                cv++
                            );
                            v1 = vertex(
                                (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (b - k - 1) / double(b) + faces[i].edges[1].origin.x,
                                (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (b - k - 1) / double(b) + faces[i].edges[1].origin.y,
                                (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (b - k - 1) / double(b) + faces[i].edges[1].origin.z,
                                cv++
                            );
                        } else if (faces[i].edges[1].origin == faces[i].edges[0].origin) {
                            v = faces[i].edges[0].origin;
                            v0 = vertex(
                                (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (k + 1) / float(b) + faces[i].edges[0].origin.x,
                                (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (k + 1) / double(b) + faces[i].edges[0].origin.y,
                                (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (k + 1) / double(b) + faces[i].edges[0].origin.z,
                                cv++
                            );
                            v1 = vertex(
                                (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (k + 1) / double(b) + faces[i].edges[1].origin.x,
                                (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (k + 1) / double(b) + faces[i].edges[1].origin.y,
                                (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (k + 1) / double(b) + faces[i].edges[1].origin.z,
                                cv++
                            );
                        } else if (faces[i].edges[0].end == faces[i].edges[1].origin) {
                            v = faces[i].edges[0].end;
                            v0 = vertex(
                                (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (b - k - 1) / float(b) + faces[i].edges[0].origin.x,
                                (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (b - k - 1) / double(b) + faces[i].edges[0].origin.y,
                                (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (b - k - 1) / double(b) + faces[i].edges[0].origin.z,
                                cv++
                            );
                            v1 = vertex(
                                (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (k + 1) / double(b) + faces[i].edges[1].origin.x,
                                (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (k + 1) / double(b) + faces[i].edges[1].origin.y,
                                (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (k + 1) / double(b) + faces[i].edges[1].origin.z,
                                cv++
                            );
                        }
                    } else if (faces[i].edges[0].end == faces[i].edges[1].end) {
                        v = faces[i].edges[0].end;
                        v0 = vertex(
                            (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (b - k - 1) / double(b) + faces[i].edges[0].origin.x,
                            (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (b - k - 1) / double(b) + faces[i].edges[0].origin.y,
                            (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (b - k - 1) / double(b) + faces[i].edges[0].origin.z,
                            cv++
                        );
                        v1 = vertex(
                            (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (b - k - 1) / double(b) + faces[i].edges[1].origin.x,
                            (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (b - k - 1) / double(b) + faces[i].edges[1].origin.y,
                            (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (b - k - 1) / double(b) + faces[i].edges[1].origin.z,
                            cv++
                        );
                    }

                    // Filtro duplicati per v0 e v1
                    tuple<double, double, double> key_v0(v0.x, v0.y, v0.z);
                    if (vertex_map.count(key_v0)) {
                        v0.id = vertex_map[key_v0];
                    } else {
                        v0.id = next_id++;
                        vertex_map[key_v0] = v0.id;
                        ver_1.push_back(v0);
                    }
                    tuple<double, double, double> key_v1(v1.x, v1.y, v1.z);
                    if (vertex_map.count(key_v1)) {
                        v1.id = vertex_map[key_v1];
                    } else {
                        v1.id = next_id++;
                        vertex_map[key_v1] = v1.id;
                        ver_1.push_back(v1);
                    }

                    Triangle[k][0] = v0;
                    Triangle[k][k+1] = v1;

                    for (size_t j = 0; j < k; j++) {
                        double t = double(j+1) / (k + 1);
                        vertex v2(v0.x + t*(v1.x-v0.x), v0.y + t*(v1.y-v0.y), v0.z + t*(v1.z-v0.z), cv++);
                        tuple<double, double, double> key_v2(v2.x, v2.y, v2.z);
                        if (vertex_map.count(key_v2)) {
                            v2.id = vertex_map[key_v2];
                        } else {
                            v2.id = next_id++;
                            vertex_map[key_v2] = v2.id;
                            ver_1.push_back(v2);
                        }
                        Triangle[k][j+1] = v2;
                    }

                    if (k > 0) {
                        for (size_t j = 0; j < k + 1 ; j++) {
                            if (j + 1 >= Triangle[k].size() || j >= Triangle[k-1].size()) continue;
                            Edge e0(Triangle[k-1][j], Triangle[k][j], id_edge++, count+1);
                            Edge e1(Triangle[k][j], Triangle[k][j+1], id_edge++, count+1);
                            Edge e2(Triangle[k][j+1], Triangle[k-1][j], id_edge++, count+1);
                            count+=3;
                            vector<Edge> e = {e0, e1, e2};
                            vector<vertex> ver = {Triangle[k-1][j], Triangle[k][j], Triangle[k][j+1]};
                            if (j>0) {
                                vector<vertex> hidden = {Triangle[k-1][j], Triangle[k-1][j-1], Triangle[k][j]};
                                Edge _e1(Triangle[k][j], Triangle[k-1][j-1], -1);
                                Edge _e2(Triangle[k-1][j-1], Triangle[k-1][j], -1);
                                vector<Edge> _e = {_e1, _e2, e0};
                                faces_1.push_back(Face(hidden, _e, id_faces++, 0));
                            }
                            faces_1.push_back(Face(ver, e, id_faces++, 0));
                            edges_1.push_back(e0);
                            edges_1.push_back(e1);
                            edges_1.push_back(e2);
                        }
                    } else if(k==0) {
                        tuple<double, double, double> key_v(v.x, v.y, v.z);
                        if (vertex_map.count(key_v)) {
                            v.id = vertex_map[key_v];
                        } else {
                            v.id = next_id++;
                            vertex_map[key_v] = v.id;
                            ver_1.push_back(v);
                        }

                        Edge e0(v, v0, id_edge++, count+1);
                        Edge e1(v0, v1, id_edge++, count+1);
                        Edge e2(v1, v, id_edge++, count+1);
                        count+=3;
                        vector<Edge> e = {e0, e1, e2};
                        vector<vertex> ver = {v, v0, v1};
                        faces_1.push_back(Face(ver, e, id_faces++, 0));
                        edges_1.push_back(e0);
                        edges_1.push_back(e1);
                        edges_1.push_back(e2);
                    }
                }
            }
        }
    }
    edges = edges_1;
    faces = faces_1;
    vertices = ver_1;
}

    void _Polyhedron::First_Triangulation(
   	vector<vertex>& vertices,
    vector<Edge>& edges,
    Face& originalFace,
    unsigned int& cv,
    unsigned int& ce,
    unsigned int& cf,
    vector<Face>& fill,unsigned int code, 
	size_t k, size_t j, 
   map<size_t, map<size_t,vertex>>& _Triangle,
   map<size_t, map<size_t,vertex>>& _TriangleI
)
{
    const double heightRatio = sqrt(3.0) / 4.0;
    
    Edge v_0 = originalFace.edges[0];
    Edge v_1 = originalFace.edges[1];
    Edge v_2 = originalFace.edges[2];

    vertex shared01, shared12, shared02;
    if (v_1.end != v_0.end) {
        if (v_1.end == v_0.origin || v_1.origin == v_0.origin)
            shared01 = v_0.origin;
        else if (v_0.end == v_1.origin)
            shared01 = v_0.end;
    } else {
        shared01 = v_0.end;
    }

    if (v_2.end != v_1.end) {
        if (v_2.end == v_1.origin || v_2.origin == v_1.origin)
            shared12 = v_1.origin;
        else if (v_1.end == v_2.origin)
            shared12 = v_1.end;
    } else {
        shared12 = v_1.end;
    }

    if (v_2.end != v_0.end) {
        if (v_2.end == v_0.origin || v_2.origin == v_0.origin)
            shared02 = v_0.origin;
        else if (v_0.end == v_2.origin)
            shared02 = v_0.end;
    } else {
        shared02 = v_0.end;
    }
    vertex center = vertex(
        0.5 * (v_2.origin.x + v_2.end.x),
        0.5 * (v_2.origin.y + v_2.end.y),
        0.5 * (v_2.origin.z + v_2.end.z),
        cv++
    );
    vertex apex = vertex(
        (shared01.x - center.x) * heightRatio + center.x,
        (shared01.y - center.y) * heightRatio + center.y,
        (shared01.z - center.z) * heightRatio + center.z,
        cv++
    );
    vertex midEdge0 = vertex(
        0.5 * (v_0.origin.x + v_0.end.x),
        0.5 * (v_0.origin.y + v_0.end.y),
        0.5 * (v_0.origin.z + v_0.end.z),
        cv++
    );
    vertex midEdge1 = vertex(
        0.5 * (v_1.origin.x + v_1.end.x),
        0.5 * (v_1.origin.y + v_1.end.y),
        0.5 * (v_1.origin.z + v_1.end.z),
        cv++
    );
	 vertices.insert(vertices.end(), {center, apex, midEdge0, midEdge1, shared02, shared12, shared01});

 if(code==0){

    Edge edge_midEdge0_to_shared02(midEdge0, shared02, ce++);
	Edge edge_shared02_to_apex(shared02, apex, ce++);
	Edge edge_apex_to_midEdge0(apex, midEdge0, ce++);
	Edge edge_shared01_to_midEdge0(shared01, midEdge0, ce++);
	Edge edge_apex_to_shared01(apex, shared01, ce++);
	Edge edge_shared01_to_midEdge1(shared01, midEdge1, ce++);
	Edge edge_midEdge1_to_apex(midEdge1, apex, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_midEdge1_to_shared12(midEdge1, shared12, ce++);
    Face face2({midEdge0, shared02, apex}, {edge_midEdge0_to_shared02, edge_shared02_to_apex, edge_apex_to_midEdge0}, cf++, 0);
    Face face3({shared01, midEdge0, apex}, {edge_shared01_to_midEdge0, reverseEdge(edge_apex_to_midEdge0), edge_apex_to_shared01}, cf++, 0);
    Face face4({shared01, midEdge1, apex}, {edge_shared01_to_midEdge1, edge_midEdge1_to_apex, edge_apex_to_shared01}, cf++, 0);
    Face face5({midEdge1, apex, shared12}, {edge_midEdge1_to_apex, edge_apex_to_shared12, reverseEdge(edge_midEdge1_to_shared12)}, cf++, 0);
	fill.insert(fill.end(),{face2,face3,face4,face5});
	    _Triangle[k][j]=apex;
}
	if(code==11){
	Edge edge_center_to_apex(center, apex, ce++);
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_midEdge0(apex, midEdge0, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_midEdge1_to_apex(midEdge1, apex, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    Edge edge_shared01_to_midEdge0(shared01, midEdge0, ce++);
    Edge edge_midEdge0_to_shared02(midEdge0, shared02, ce++);
    Edge edge_shared01_to_midEdge1(shared01, midEdge1, ce++);
    Edge edge_midEdge1_to_shared12(midEdge1, shared12, ce++);
    Edge edge_shared12_to_center(shared12, center, ce++);
    Edge edge_shared02_to_center(shared02, center, ce++);
    edges.insert(edges.end(), {
        edge_center_to_apex, edge_apex_to_shared01, edge_apex_to_midEdge0,
        edge_apex_to_shared12, edge_midEdge1_to_apex, edge_shared02_to_apex,
        edge_shared01_to_midEdge0, edge_midEdge0_to_shared02,
        edge_shared01_to_midEdge1, edge_midEdge1_to_shared12,
        edge_shared12_to_center, edge_shared02_to_center
    });
    Face face2({midEdge0, shared02, apex}, {edge_midEdge0_to_shared02, edge_shared02_to_apex, edge_apex_to_midEdge0}, cf++, 0);
    Face face3({shared01, midEdge0, apex}, {edge_shared01_to_midEdge0, reverseEdge(edge_apex_to_midEdge0), edge_apex_to_shared01}, cf++, 0);
     fill.insert(fill.end(), {face2, face3});
         _Triangle[k][j]=apex;
	}
	if(code==13){
	
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    edges.insert(edges.end(), {
		edge_apex_to_shared01,edge_apex_to_shared12,edge_shared02_to_apex
    });
        _Triangle[k][j]=apex;
	}
	if(code==12){	
	Edge edge_center_to_apex(center, apex, ce++);
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_midEdge0(apex, midEdge0, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_midEdge1_to_apex(midEdge1, apex, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    Edge edge_shared01_to_midEdge0(shared01, midEdge0, ce++);
    Edge edge_midEdge0_to_shared02(midEdge0, shared02, ce++);
    Edge edge_shared01_to_midEdge1(shared01, midEdge1, ce++);
    Edge edge_midEdge1_to_shared12(midEdge1, shared12, ce++);
    Edge edge_shared12_to_center(shared12, center, ce++);
    Edge edge_shared02_to_center(shared02, center, ce++);
    edges.insert(edges.end(), {
        edge_center_to_apex, edge_apex_to_shared01, edge_apex_to_midEdge0,
        edge_apex_to_shared12, edge_midEdge1_to_apex, edge_shared02_to_apex,
        edge_shared01_to_midEdge0, edge_midEdge0_to_shared02,
        edge_shared01_to_midEdge1, edge_midEdge1_to_shared12,
        edge_shared12_to_center, edge_shared02_to_center
    });
    Face face4({shared01, midEdge1, apex}, {edge_shared01_to_midEdge1, edge_midEdge1_to_apex, edge_apex_to_shared01}, cf++, 0);
    Face face5({midEdge1, apex, shared12}, {edge_midEdge1_to_apex, edge_apex_to_shared12, reverseEdge(edge_midEdge1_to_shared12)}, cf++, 0);
     fill.insert(fill.end(), {face4, face5});
	     _Triangle[k][j]=apex; 	
	}
	if(code==21){
	Edge edge_center_to_apex(center, apex, ce++);
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_midEdge0(apex, midEdge0, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    Edge edge_shared01_to_midEdge0(shared01, midEdge0, ce++);
    Edge edge_midEdge0_to_shared02(midEdge0, shared02, ce++);
    Edge edge_shared12_to_center(shared12, center, ce++);
    Edge edge_shared02_to_center(shared02, center, ce++);
    edges.insert(edges.end(), {
edge_center_to_apex,edge_apex_to_shared01,edge_apex_to_midEdge0,
edge_apex_to_shared12,edge_shared02_to_apex,edge_shared01_to_midEdge0,
edge_midEdge0_to_shared02,edge_shared12_to_center,edge_shared02_to_center
    });
    Face face0({center, shared12, apex}, {edge_apex_to_shared12, edge_shared12_to_center, edge_center_to_apex}, cf++, 0);
    printFace(face0);
    Face face1({shared02, center, apex}, {edge_shared02_to_center, edge_center_to_apex, edge_shared02_to_apex}, cf++, 0);
    Face face2({midEdge0, shared02, apex}, {edge_midEdge0_to_shared02, edge_shared02_to_apex, edge_apex_to_midEdge0}, cf++, 0);
    Face face3({shared01, midEdge0, apex}, {edge_shared01_to_midEdge0, reverseEdge(edge_apex_to_midEdge0), edge_apex_to_shared01}, cf++, 0);
    fill.insert(fill.end(), {face0, face1, face2, face3});
        _Triangle[k][j]=apex;
	}
	if(code==22){
	    Edge edge_center_to_apex(center, apex, ce++);
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_midEdge1_to_apex(midEdge1, apex, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    Edge edge_midEdge0_to_shared02(midEdge0, shared02, ce++);
    Edge edge_shared01_to_midEdge1(shared01, midEdge1, ce++);
    Edge edge_midEdge1_to_shared12(midEdge1, shared12, ce++);
    Edge edge_shared12_to_center(shared12, center, ce++);
    Edge edge_shared02_to_center(shared02, center, ce++);
    edges.insert(edges.end(), {
    edge_apex_to_shared01,edge_apex_to_shared01,edge_apex_to_shared12
    ,edge_midEdge1_to_apex	,edge_shared02_to_apex	,edge_midEdge0_to_shared02 	
    ,edge_shared01_to_midEdge1	,edge_midEdge1_to_shared12,edge_shared12_to_center,edge_shared02_to_center
    });
    Face face0({center, shared12, apex}, {edge_apex_to_shared12, edge_shared12_to_center, edge_center_to_apex}, cf++, 0);
    printFace(face0);
    Face face1({shared02, center, apex}, {edge_shared02_to_center, edge_center_to_apex, edge_shared02_to_apex}, cf++, 0);
    Face face4({shared01, midEdge1, apex}, {edge_shared01_to_midEdge1, edge_midEdge1_to_apex, edge_apex_to_shared01}, cf++, 0);
    Face face5({midEdge1, apex, shared12}, {edge_midEdge1_to_apex, edge_apex_to_shared12, reverseEdge(edge_midEdge1_to_shared12)}, cf++, 0);
    fill.insert(fill.end(), {face0, face1,face4,face5});
        _Triangle[k][j]=apex;
	}
	if(code==23){
	Edge edge_center_to_apex(center, apex, ce++);
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    Edge edge_shared12_to_center(shared12, center, ce++);
    Edge edge_shared02_to_center(shared02, center, ce++);
    edges.insert(edges.end(), {
    	edge_center_to_apex,
    	edge_apex_to_shared01,
    	edge_apex_to_shared12,
    	edge_shared02_to_apex,
		edge_shared12_to_center
		,edge_shared02_to_center
    });
    Face face0({center, shared12, apex}, {edge_apex_to_shared12, edge_shared12_to_center, edge_center_to_apex}, cf++, 0);
    Face face1({shared02, center, apex}, {edge_shared02_to_center, edge_center_to_apex, edge_shared02_to_apex}, cf++, 0);
    fill.insert(fill.end(), {face0, face1});
        _Triangle[k][j]=apex;
	}
	if(code==14){
    Edge edge_apex_to_shared01(apex, shared01, ce++);
    Edge edge_apex_to_shared12(apex, shared12, ce++);
    Edge edge_shared02_to_apex(shared02, apex, ce++);
    edges.insert(edges.end(), {
	edge_apex_to_shared01,
	edge_apex_to_shared12,
	edge_shared02_to_apex});
	_TriangleI[k][j]=apex;
	
	}
}


                            
void _Polyhedron::Assembler(
					//	map<size_t, map<size_t,vertex>>& Triangle,
					//	map<size_t, map<size_t,vertex>>& TC,
						map<size_t,vector<vertex>>& Triangle,
					    map<size_t, map<size_t,vertex>>& TC,
						map<size_t, map<size_t,vertex>>& TCII,
						size_t k, 
						size_t j , 
						unsigned int code,
						vector<Face>& fill,
						vector<Edge>& edges,
						unsigned int& cf,
						unsigned int& ce)
		{
		
		if(code==0){
			Edge e0(TC[k][j], TC[k+1][j+1], ce++);
			Edge e1(TC[k+1][j+1], Triangle[k][j], ce++ );
			Edge e2(TC[k+1][j+1], Triangle[k][j+1], ce++);
			Edge e3(TC[k][j], Triangle[k][j], ce++);
			Edge e4(TC[k][j], Triangle[k][j+1], ce++);
			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);
			Face f0({TC[k][j], TC[k+1][j+1], Triangle[k][j+1]}, {e0,e2,reverseEdge(e4)}, cf++,0);
			Face f1({TC[k][j], TC[k+1][j+1],Triangle[k][j]} , {e0,e1,reverseEdge(e3)}, cf++,0);
			fill.push_back(f0);
			fill.push_back(f1);
			}				
		if(code==11){
			Edge e0 (TCII[k][j+1],Triangle[k][j],ce++);
			Edge e1 (Triangle[k][j],TC[k][j],ce++);
			Edge e2 (TC[k][j], TCII[k][j+1],ce++);
			
			Edge e3 (TCII[k+1][j+1],TC[k][j],ce++);
			Edge e4 (TC[k][j],Triangle[k][j+1],ce++);
			Edge e5 (Triangle[k][j+1],TCII[k+1][j+1],ce++);
			
			Face f0({TCII[k][j+1],Triangle[k][j],TC[k][j]},{e0,e1,e2},cf++,0);
			Face f1({TCII[k+1][j+1],TC[k][j],Triangle[k][j+1]},{e3,e4,e5},cf++,0);
			
			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);				
			fill.push_back(f0);fill.push_back(f1);
		}
		if(code==12){
			Edge e0 (TCII[k][j],Triangle[k][j],ce++);
			Edge e1 (Triangle[k][j],TC[k][j],ce++);
			Edge e2 (TC[k][j], TCII[k][j],ce++);
	
			Edge e3 (TCII[k+1][j],TC[k][j],ce++);
			Edge e4 (TC[k][j],Triangle[k][j+1],ce++);
			Edge e5 (Triangle[k][j+1],TCII[k+1][j],ce++);

			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);
			
			Face f0({TCII[k][j],Triangle[k][j],TC[k][j]},{e0,e1,e2},cf++,0);
			Face f1({TCII[k+1][j],TC[k][j],Triangle[k][j+1]},{e3,e4,e5},cf++,0);
			fill.push_back(f0);fill.push_back(f1);
		}		
		if(code==13){
			Edge e0(TCII[k][j],TC[k][j],ce++);
			Edge e1(TC[k][j],Triangle[k][j],ce++);
			Edge e2(Triangle[k][j],TCII[k][j],ce++);
			Face f0({TCII[k][j],TC[k][j],Triangle[k][j]},{e0,e1,e2},cf++,0);
			
			Edge e3(TC[k][j],TCII[k+1][j+1],ce++);
			Edge e4(TCII[k+1][j+1],Triangle[k][j],ce++);
			Edge e5(Triangle[k][j],TC[k][j],ce++);
			Face f1({TC[k][j],TCII[k+1][j+1],Triangle[k][j]},{e3,e4,e5},cf++,0);
			
			Edge e6(TC[k][j],TCII[k+1][j+1],ce++);
			Edge e7(TCII[k+1][j+1],Triangle[k][j+1],ce++);
			Edge e8(Triangle[k][j+1],TC[k][j],ce++);
			Face f2({TC[k][j],TCII[k+1][j+1],Triangle[k][j+1]},{e6,e7,e8},cf++,0);
			
			Edge e9(TC[k][j],Triangle[k][j+1],ce++);
			Edge e10(Triangle[k][j+1],TCII[k][j+1],ce++);
			Edge e11(TCII[k][j+1],TC[k][j],ce++);
			Face f3({TC[k][j],Triangle[k][j+1],TCII[k][j+1]},{e9,e10,e11},cf++,0);

			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);
			add_unique_edge(edges,e6);
			add_unique_edge(edges,e7);
			add_unique_edge(edges,e8);
			add_unique_edge(edges,e9);
			add_unique_edge(edges,e10);			
			add_unique_edge(edges,e11);
			fill.push_back(f0);	fill.push_back(f1);fill.push_back(f2);fill.push_back(f3);
		}
		if(code==21){
			Edge e0(TCII[k][j+1],TC[k][j],ce++);
			Edge e1(TC[k][j],Triangle[k][j+1],ce++);
			Edge e2(Triangle[k][j+1],TCII[k][j+1],ce++);
			Face f0({TCII[k][j+1],TC[k][j],Triangle[k][j+1]},{e0,e1,e2},cf++,0);
			Edge e3(Triangle[k-1][j],TC[k][j],ce++);
			Edge e4(TC[k][j],TCII[k][j+1],ce++);
			Edge e5(TCII[k][j+1],Triangle[k-1][j]);
			Face f1({Triangle[k-1][j],TC[k][j],TCII[k][j+1]},{e3,e4,e5},cf++,0);
			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);	
			fill.push_back(f0);fill.push_back(f1);
			
		}
		if(code==22){
			Edge e0(TC[k][j],TCII[k][j+1],ce++);
			Edge e1(TCII[k][j+1],Triangle[k-1][j],ce++);
			Edge e2(Triangle[k-1][j],TC[k][j],ce++);
			Face f0({TC[k][j],TCII[k][j+1],Triangle[k-1][j]},{e0,e1,e2},cf++,0);
			
			
			Edge e3(TC[k][j],TCII[k][j+1],ce++);
			Edge e4(TCII[k][j+1],Triangle[k][j+1],ce++);
			Edge e5(Triangle[k][j+1],TC[k][j],ce++);
			Face f1({TC[k][j],TCII[k][j+1],Triangle[k][j+1]},{e3,e4,e5},cf++,0);
				
			
			Edge e6(TC[k][j],TCII[k-1][j],ce++);
			Edge e7(TCII[k-1][j],Triangle[k-1][j],ce++);
			Edge e8(Triangle[k-1][j],TC[k][j],ce++);
			Face f2({TC[k][j],TCII[k-1][j],Triangle[k-1][j]},{e6,e7,e8},cf++,0);	
	
			
			Edge e9(TC[k][j],TCII[k][j],ce++);
			Edge e10(TCII[k][j],Triangle[k][j],ce++);
			Edge e11(Triangle[k][j],TC[k][j],ce++);
			Face f3({TC[k][j],TCII[k][j],Triangle[k][j]},{e9,e10,e11},cf++,0);

			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);
			add_unique_edge(edges,e6);
			add_unique_edge(edges,e7);
			add_unique_edge(edges,e8);
			add_unique_edge(edges,e9);
			add_unique_edge(edges,e10);			
			add_unique_edge(edges,e11);	
			fill.push_back(f0);fill.push_back(f1);fill.push_back(f2);fill.push_back(f3);
		
		}
		if(code==23){
			Edge e0(TC[k][j],TCII[k][j],ce++);
			Edge e1(TCII[k][j],Triangle[k-1][j],ce++);
			Edge e2(Triangle[k-1][j],TC[k][j],ce++);
			Face f0({TC[k][j],TCII[k][j],Triangle[k-1][j]},{e0,e1,e2},cf++,0);
			
			Edge e3(TC[k][j],TCII[k][j],ce++);
			Edge e4(TCII[k][j],Triangle[k][j],ce++);
			Edge e5(Triangle[k][j],TC[k][j],ce++);
			Face f1({TC[k][j],TCII[k][j],Triangle[k][j]},{e3,e4,e5},cf++,0);
			add_unique_edge(edges,e0);
			add_unique_edge(edges,e1);
			add_unique_edge(edges,e2);
			add_unique_edge(edges,e3);
			add_unique_edge(edges,e4);			
			add_unique_edge(edges,e5);
			fill.push_back(f0); fill.push_back(f1);
				
		}
		}
    void _Polyhedron::Triangulation_2() {
    vector<Edge> edges_1;
    vector<Face> faces_1;
    vector<vertex> ver_1;
    vector<vertex> ver_2;
	vector<Face> faces_2;
	vector<Edge> edges_2;	
	unsigned int id_f1=0;
	unsigned int id_v1=0;
	unsigned int id_e1=0;	
        if (p == 3) {
            unsigned int T = b * b + c * b + c * c;
            switch (q) {
                case 3:
                    NumVer = 2 * T + 2; NumEdg = 6 * T; NumFcs = 4 * T;
                    break;
                case 4:
                    NumVer = 4 * T + 2; NumEdg = 12 * T; NumFcs = 8 * T;
                    break;
                case 5:
                    NumVer = 10 * T; NumEdg = 30 * T; NumFcs = 20 * T;
                    break;
            }
            faces_1.reserve(NumFcs);
            edges_1.reserve(NumEdg);
            edges.reserve(NumEdg);
            vertices.reserve(NumVer);
            ver_1.reserve(NumVer);
            ver_2.reserve(NumVer);
            edges_2.reserve(NumEdg);
            faces_2.reserve(NumFcs);
            unsigned int id_edge = 0;
            unsigned int id_faces = 0;
            size_t j=0;
            vertex v;
            vertex v1;
            map<size_t, map<size_t,vertex>> _Triangle_;
            vector<size_t> index_edges(3);
            size_t cv = 0;
        	if(b==0){
        		b=c;
			}else{
				c=b;
			}					
  			for (size_t i = 0; i < faces.size(); i++) {
  			map<size_t, vector<vertex>> Triangle;
  			map<size_t, map<size_t,vertex>> Triangle_ClassII;
            map<size_t, map<size_t,vertex>> Triangle_ClassII_I;
               for (size_t k = 0; k < b ; k++) {
                    Triangle[k].resize(k + 2);
                    
					vertex v0;
					vertex v_origin;				  		
			if (faces[i].edges[1].end != faces[i].edges[0].end) 
			{
    		if (faces[i].edges[1].end == faces[i].edges[0].origin) {
        	v = faces[i].edges[0].origin;
        	v0 = vertex(
            (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (k + 1) / float(b) + faces[i].edges[0].origin.x,
            (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (k + 1) / double(b) + faces[i].edges[0].origin.y,
            (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (k + 1) / double(b) + faces[i].edges[0].origin.z,
            cv++
        	)	;
        	v1 = vertex(
            (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (b - k - 1) / double(b) + faces[i].edges[1].origin.x,
            (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (b - k - 1) / double(b) + faces[i].edges[1].origin.y,
            (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (b - k - 1) / double(b) + faces[i].edges[1].origin.z,
            cv++
        	);
    		} else if (faces[i].edges[1].origin == faces[i].edges[0].origin) {
        	v = faces[i].edges[0].origin;
        	v0 = vertex(
            (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (k + 1) / float(b) + faces[i].edges[0].origin.x,
            (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (k + 1) / double(b) + faces[i].edges[0].origin.y,
            (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (k + 1) / double(b) + faces[i].edges[0].origin.z,
            cv++
        	);
        	v1 = vertex(
            (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (k + 1) / double(b) + faces[i].edges[1].origin.x,
            (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (k + 1) / double(b) + faces[i].edges[1].origin.y,
            (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (k + 1) / double(b) + faces[i].edges[1].origin.z,
            cv++
        	);
    		} else if (faces[i].edges[0].end == faces[i].edges[1].origin) {
        	v = faces[i].edges[0].end;
        	v0 = vertex(
            (faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (b - k - 1) / float(b) + faces[i].edges[0].origin.x,
            (faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (b - k - 1) / double(b) + faces[i].edges[0].origin.y,
            (faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (b - k - 1) / double(b) + faces[i].edges[0].origin.z,
            cv++
        	);
     	  	v1 = vertex(
            (faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (k + 1) / double(b) + faces[i].edges[1].origin.x,
            (faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (k + 1) / double(b) + faces[i].edges[1].origin.y,
            (faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (k + 1) / double(b) + faces[i].edges[1].origin.z,
            cv++
        	);
    		}
			}
			else if (faces[i].edges[0].end == faces[i].edges[1].end)
			{
    		v = faces[i].edges[0].end;
    		v0 = vertex(
        	(faces[i].edges[0].end.x - faces[i].edges[0].origin.x) * (b - k - 1) / double(b) + faces[i].edges[0].origin.x,
        	(faces[i].edges[0].end.y - faces[i].edges[0].origin.y) * (b - k - 1) / double(b) + faces[i].edges[0].origin.y,
        	(faces[i].edges[0].end.z - faces[i].edges[0].origin.z) * (b - k - 1) / double(b) + faces[i].edges[0].origin.z,
        	cv++
    		);
    		v1 = vertex(
        	(faces[i].edges[1].end.x - faces[i].edges[1].origin.x) * (b - k - 1) / double(b) + faces[i].edges[1].origin.x,
        	(faces[i].edges[1].end.y - faces[i].edges[1].origin.y) * (b - k - 1) / double(b) + faces[i].edges[1].origin.y,
        	(faces[i].edges[1].end.z - faces[i].edges[1].origin.z) * (b - k - 1) / double(b) + faces[i].edges[1].origin.z,
        	cv++
    		);	
			}
            Triangle[k][0] = v0;
            Triangle[k][k+1] = v1;
                for (size_t j = 0; j < k; j++) {				
                double t = double(j+1)/(k + 1);
					vertex v2(
    				v0.x + t * (v1.x - v0.x),
 				   	v0.y + t * (v1.y - v0.y),
    				v0.z + t * (v1.z - v0.z),
    				cv++
				);
                Triangle[k][j+1] = v2; 
                }          
                if (k > 0) {
                    for (j = 0; j < k + 1 ; j++) {   
                         if (j + 1 >= Triangle[k].size() || j >= Triangle[k - 1].size()) continue;
						    Edge e1(Triangle[k][j], Triangle[k][j + 1], id_edge++);	
							Edge e2(Triangle[k][j + 1], Triangle[k - 1][j], id_edge++);
                            Edge e0(Triangle[k - 1][j], Triangle[k][j], id_edge++);
                            vector<Edge> _e = {e0, e1, e2};
                            vector<vertex> ver = {Triangle[k - 1][j], Triangle[k][j], Triangle[k][j + 1]};
                            if(j>0){		
							vector<vertex> hidden = {Triangle[k - 1][j], Triangle[k-1][j-1], Triangle[k][j]};
						    Edge _e1(Triangle[k][j], Triangle[k - 1][j-1], -1);
							Edge _e2(Triangle[k-1][j-1], Triangle[k-1][j], -1);
                            vector<Edge> _e = {e0,_e2,_e1};
                            Face fb_1(hidden, _e, id_faces++, 0);
                            faces_1.push_back(fb_1);
                            First_Triangulation(ver_2,edges_2,fb_1,id_v1,id_e1,id_f1,faces_2,14,k,j,Triangle_ClassII, Triangle_ClassII_I);
						    }
						    Face fb(ver, _e, id_faces++, 0);
                            faces_1.push_back(fb);
							if(k<b-1){
									if(j==0){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,11,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}else if(j>0 && j<k){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,12,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}else if(j==k){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,13,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}			
							}
							if(k==b-1){
									if(j==0){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,21,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}else if(j>0 && j<k){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,23,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}else if(j==k){
								First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,22,k,j,Triangle_ClassII, Triangle_ClassII_I);										
									}								
							}	
                            edges_1.push_back(e1); 
							edges_1.push_back(e2);
							edges_1.push_back(e0);
                            ver_1.push_back(Triangle[k - 1][j]);
                            ver_1.push_back(Triangle[k][j]);
                            ver_1.push_back(Triangle[k][j + 1]);
                        }
                    } 
					else if(k==0) {
                    	vertex c;
                        Edge e0(v, v0, id_edge++);
                        Edge e1(v0, v1, id_edge++);
                        Edge e2(v1, v, id_edge++);
                        vector<Edge> e = {e0, e1, e2};
                        vector<vertex> ver = {v, v0, v1};
                        Face fb(ver, e, id_faces++, 0);
                        faces_1.push_back(fb);                    
                        First_Triangulation(ver_2,edges_2,fb,id_v1,id_e1,id_f1,faces_2,22,k,j,Triangle_ClassII, Triangle_ClassII_I);
                        edges_1.push_back(e0);
                        edges_1.push_back(e1);
                        edges_1.push_back(e2);
                        ver_1.push_back(v);
                        ver_1.push_back(v0);
                        ver_1.push_back(v1);
                    }
                }
			for(size_t k=0; k<b; k++){
				for(size_t j=0; j<k;j++){
					if(k==0){
						Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,0,faces_2,edges_2,id_f1,id_e1);
					}
					if(k>0 && k<b-1){
						if(j==0){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,11,faces_2,edges_2,id_f1,id_e1);
						}
						if(j>0 && j<k){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,13,faces_2,edges_2,id_f1,id_e1);
						}
						if(j==k){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,12,faces_2,edges_2,id_f1,id_e1);
						}
					}
					if(k==b-1 && k>0){
						if(j==0){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,21,faces_2,edges_2,id_f1,id_e1);
						}
						if(j>0 && j<k-1){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,23,faces_2,edges_2,id_f1,id_e1);	
						}
						if(j==k-1){
							Assembler(Triangle,Triangle_ClassII,Triangle_ClassII_I,k,j,22,faces_2,edges_2,id_f1,id_e1);
						}
					}
				}
			}		
            }//loop faces.   
        }
    edges = edges_2;
    faces = faces_2;
    vertices = ver_2;
        }
	void _Polyhedron::OverAll_Triangulation(){
			if ((b == 0 && c > 0) || (b >0 && c == 0)) {
				_Polyhedron::Triangulation();
			}
			else if(b>0 && c>0 && b==c){
				_Polyhedron::Triangulation_2();
			}
            for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i].id = i;
            std::vector<Edge> unique_edges;
            std::set<std::pair<unsigned int, unsigned int>> edge_set;

            for (const auto& e : edges) {
                unsigned int a = e.origin.id;
                unsigned int b = e.end.id;
                if (a > b) std::swap(a, b);
                std::pair<unsigned int, unsigned int> edge_pair(a, b);
                
                if (edge_set.find(edge_pair) == edge_set.end()) {
                    edge_set.insert(edge_pair);
                    unique_edges.push_back(e);
                }
}

edges = unique_edges;
		}	};
	void _Polyhedron::GenerateDual() {
    	map<unsigned int, vector<unsigned int>> M;  
    	map<unsigned int, vertex> centroid;    
    	size_t n = faces[0].vertices.size();     
    	for (const auto& face : faces) {
        double cx = 0, cy = 0, cz = 0;
        for (const auto& v : face.vertices) {
            cx += v.x;
            cy += v.y;
            cz += v.z;
        }
        centroid[face.id] = vertex(cx / n, cy / n, cz / n);
        centroid[face.id].normalize();
    	}
   	 	for (size_t i = 0; i < faces.size(); ++i) {
        	for (size_t j = i + 1; j < faces.size(); ++j) {
            	int shared = 0;
            		for (const auto& vi : faces[i].vertices) {
                		for (const auto& vj : faces[j].vertices) {
                    		if (vi == vj) shared++;
                		}
            			}
            			if (shared == 2) {
                			M[faces[i].id].push_back(faces[j].id);
                			M[faces[j].id].push_back(faces[i].id);
            			}
        		}
    		}
    		vector<Edge> dual_edges;
    		unsigned int id_edg = 0;
 	   for (size_t i = 0; i < faces.size(); ++i) {
    	    unsigned int fid = faces[i].id;
        		for (unsigned int adj_id : M[fid]) {
           		 bool already_added = false;
            		for (const auto& e : dual_edges) {
                		if ((e.origin == centroid[fid] && e.end == centroid[adj_id]) ||
                    		(e.end == centroid[fid] && e.origin == centroid[adj_id])) {
                    			already_added = true;
                    			break;
                }
            }
            if (!already_added) {
                Edge e(centroid[fid], centroid[adj_id], id_edg++);
                dual_edges.push_back(e);
            }
        }
    }
    for (const auto& e : dual_edges) {
        cout << "Edge " << e.id << ": (" << e.origin.x << ", " << e.origin.y << ", " << e.origin.z
             << ") -> (" << e.end.x << ", " << e.end.y << ", " << e.end.z
             << ") length: " << e.length() << "\n";
    }
	//bulk of the face.
	map<unsigned int, vector<unsigned int>> vertex_to_faces;
	for (const auto& face : faces) {
    for (const auto& v : face.vertices) {
        vertex_to_faces[v.id].push_back(face.id);}
	}
	}
    	MatrixXd ConvertVerticesToEigen(const vector<vertex>& vertices) {
        MatrixXd points(3, vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i) {
            points(0, i) = vertices[i].x;
            points(1, i) = vertices[i].y;
            points(2, i) = vertices[i].z;
        }
        return points;
    }
	MatrixXi ConvertEdgesToEigen(const std::vector<Edge>& edges) {
        MatrixXi segments(2, edges.size());
        for (size_t i = 0; i < edges.size(); ++i) {
            segments(0, i) = edges[i].origin.id;
            segments(1, i) = edges[i].end.id;
        }
        return segments;
    }

	}	
	
	

