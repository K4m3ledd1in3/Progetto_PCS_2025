#include "Eigen/Eigen"
#include "Polygon.hpp"
#include <iostream>
#include "Polyhedron.hpp"
#include <limits>
#include <vector>
#include <set>
#include <fstream>
#include <cmath>
#include <queue>
using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;
namespace PolyhedronLibrary{

    
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


    vector<PolygonalLibrary::vertex> ShortestPath(
    vector<PolygonalLibrary::vertex>& vertices,
    vector<PolygonalLibrary::Edge>& edges,
    unsigned int id1, unsigned int id2) {

    const int N = vertices.size();

    // Verifica validità ID
    if (id1 >= N || id2 >= N) {
        std::cerr << "Errore: ID vertici non valido.\n";
        return {};
    }

    // Costruzione lista di adiacenza non orientata
    std::unordered_map<int, std::vector<std::pair<int, double>>> adj;
    for (const auto& edge : edges) {
        int u = edge.origin.id;
        int v = edge.end.id;
        double w = edge.length();
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w); // grafo non orientato
    }

    // Inizializzazione strutture per Dijkstra
    std::vector<double> dist(N, std::numeric_limits<double>::infinity());
    std::vector<int> prev(N, -1);
    std::vector<bool> visited(N, false);
    using P = std::pair<double, int>; // (distanza, vertice)
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;

    dist[id1] = 0.0;
    pq.push({0.0, id1});

    // Algoritmo di Dijkstra
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (visited[u]) continue;
        visited[u] = true;

        for (auto [v, w] : adj[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    // Ricostruzione del cammino minimo
    std::vector<unsigned int> path_ids;
    for (int at = id2; at != -1; at = prev[at]) {
        path_ids.push_back(at);
    }
    if (path_ids.back() != id1) {
        std::cout << "Nessun cammino trovato.\n";
        return {};
    }
    std::reverse(path_ids.begin(), path_ids.end());

    // Marcare ShortPath = 1 sui vertici
    for (auto& v : vertices) v.ShortPath = 0;
    for (int id : path_ids) vertices[id].ShortPath = 1;

    // Marcare ShortPath = 1 sugli archi
    for (auto& e : edges) e.ShortPath = 0;
    for (size_t i = 0; i + 1 < path_ids.size(); ++i) {
        int a = path_ids[i], b = path_ids[i + 1];
        bool found = false;
        for (auto& e : edges) {
            if ((e.origin.id == a && e.end.id == b) || 
                (e.origin.id == b && e.end.id == a)) {
                e.ShortPath = 1;
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "⚠ Arco mancante tra " << a << " e " << b << ", lo aggiungo forzatamente.\n";
            Edge e(vertices[a], vertices[b], edges.size());
            e.ShortPath = 1;
            edges.push_back(e);
        }

    }

    // Calcolo numero di lati e lunghezza totale
    double total_length = 0.0;
    for (size_t i = 0; i + 1 < path_ids.size(); ++i) {
        int a = path_ids[i], b = path_ids[i + 1];
        for (const auto& e : edges) {
            if ((e.origin.id == a && e.end.id == b) || 
                (e.origin.id == b && e.end.id == a)) {
                total_length += e.length();
                break;
            }
        }
    }

    std::cout << "Cammino minimo tra " << id1 << " e " << id2
              << ": " << path_ids.size() - 1
              << " lati, lunghezza totale " << total_length << "\n";

    // Costruzione del vettore di vertex da restituire
    std::vector<PolygonalLibrary::vertex> path_vertices;
    for (int id : path_ids) {
        path_vertices.push_back(vertices[id]);
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
    //______________________
        // Trova il vertice comune tra due archi
    vertex find_common_vertex(const Edge& e1, const Edge& e2) {
        if (e1.origin == e2.origin || e1.origin == e2.end)
            return e1.origin;
        else
            return e1.end;
    }

    // Dato un arco e il vertice comune, restituisce l'altro vertice
    vertex get_other_vertex(const Edge& e, const vertex& common) {
        return (e.origin == common) ? e.end : e.origin;
    }

    


    double _Polyhedron::round6(double x) {
        return round(x * 1e6) / 1e6;
    }

    void _Polyhedron::Triangulation() {
    faces_original = faces;
    vector<Edge> edges_1;
    vector<Face> faces_1;
    vector<vertex> ver_1;
    
    if(b==0){
            b=c;
        }else{
            c=b;
        }
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
        
        unsigned int count=0;
        map<tuple<double, double, double>, unsigned int> vertex_map;
        unsigned int next_id = 0;

        for (size_t i = 0; i < faces.size(); i++) {
            for (size_t k = 0; k < b ; k++) {
                Triangle[k].resize(k + 2);
                            // Prendi i primi due vertici della faccia
                vertex v = faces[i].vertices[0];
                vertex v_end_0 = faces[i].vertices[1];
                vertex v_end_1 = faces[i].vertices[2];

                // Interpolazione su edge[0]
                vertex v0(
                    (v_end_0.x - v.x) * (k + 1) / double(b) + v.x,
                    (v_end_0.y - v.y) * (k + 1) / double(b) + v.y,
                    (v_end_0.z - v.z) * (k + 1) / double(b) + v.z,
                    0 // id temporaneo
                );

                // Interpolazione su edge[1]
                vertex v1(
                    (v_end_1.x - v.x) * (k + 1) / double(b) + v.x,
                    (v_end_1.y - v.y) * (k + 1) / double(b) + v.y,
                    (v_end_1.z - v.z) * (k + 1) / double(b) + v.z,
                    0 // id temporaneo
                );
                

                // Filtro duplicati per v0 
                tuple<double, double, double> key_v0(round6(v0.x), round6(v0.y), round6(v0.z));
                if (vertex_map.count(key_v0)) {
                    v0.id = vertex_map[key_v0];
                } else {
                    v0.id = next_id++;
                    vertex_map[key_v0] = v0.id;
                    ver_1.push_back(v0);
                }
                //Filtro duplicati per v1
                tuple<double, double, double> key_v1(round6(v1.x), round6(v1.y), round6(v1.z));
                if (vertex_map.count(key_v1)) {
                    v1.id = vertex_map[key_v1];
                } else {
                    v1.id = next_id++;
                    vertex_map[key_v1] = v1.id;
                    ver_1.push_back(v1);
                }
                // Inserisco estremi nella riga
                Triangle[k][0] = v0;
                Triangle[k][k+1] = v1;
                // Inserisco punti interni
                for (size_t j = 0; j < k; j++) {
                    double t = double(j+1) / (k + 1);
                    //vertex v2(v0.x + t*(v1.x-v0.x), v0.y + t*(v1.y-v0.y), v0.z + t*(v1.z-v0.z), cv++);
                    vertex v2(v0.x + t * (v1.x - v0.x),v0.y + t * (v1.y - v0.y),v0.z + t * (v1.z - v0.z),0 );
                    tuple<double, double, double> key_v2(round6(v2.x), round6(v2.y), round6(v2.z));
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
                    
                        //Primo triangolo
                        Edge e0(Triangle[k-1][j], Triangle[k][j], id_edge++, 0);
                        Edge e1(Triangle[k][j], Triangle[k][j+1], id_edge++, 0);
                        Edge e2(Triangle[k][j+1], Triangle[k-1][j], id_edge++, 0);
                        count+=3;
                        vector<Edge> e = {e0, e1, e2};
                        vector<vertex> ver = {Triangle[k-1][j], Triangle[k][j], Triangle[k][j+1]};
                        faces_1.push_back(Face(ver, e, id_faces++, 0));
                        if (j > 0 && j - 1 < Triangle[k - 1].size()) {
                            vector<vertex> hidden = {Triangle[k - 1][j], Triangle[k - 1][j - 1], Triangle[k][j]};

                            Edge _e1(Triangle[k][j], Triangle[k - 1][j - 1], id_edge++, 0);
                            Edge _e2(Triangle[k - 1][j - 1], Triangle[k - 1][j], id_edge++, 0);
                            Edge _e3(Triangle[k - 1][j], Triangle[k][j], id_edge++, 0);

                            vector<Edge> _e = {_e1, _e2, _e3};
                            faces_1.push_back(Face(hidden, _e, id_faces++, 0));

                            add_unique_edge(edges_1, _e1);
                            add_unique_edge(edges_1, _e2);
                            add_unique_edge(edges_1, _e3);
                            }

                                                    
                        add_unique_edge(edges_1, e0);
                        add_unique_edge(edges_1, e1);
                        add_unique_edge(edges_1, e2);
                    }
                } else if(k==0) {
                    tuple<double, double, double> key_v(round6(v.x), round6(v.y), round6(v.z));
                    if (vertex_map.count(key_v)) {
                        v.id = vertex_map[key_v];
                    } else {
                        v.id = next_id++;
                        vertex_map[key_v] = v.id;
                        ver_1.push_back(v);
                    }

                    Edge e0(v, v0, id_edge++, 0);
                    Edge e1(v0, v1, id_edge++, 0);
                    Edge e2(v1, v, id_edge++, 0);
                    count+=3;
                    vector<Edge> e = {e0, e1, e2};
                    vector<vertex> ver = {v, v0, v1};
                    faces_1.push_back(Face(ver, e, id_faces++, 0));
                    add_unique_edge(edges_1, e0);
                    add_unique_edge(edges_1, e1);
                    add_unique_edge(edges_1, e2);
                }
            }
        }
    }
    edges = edges_1;
    faces = faces_1;
    vertices = ver_1;
    }

    
        void _Polyhedron::Triangulation_2() {
            // Backup dei dati di Classe 1
            vector<vertex> vertices_class1 = vertices;
        
            vector<Face> faces_class1 = faces;
            
            vector<Edge> edges_class1 = edges;

            // Contenitori per Classe 2
            vector<vertex> ver_2;
            vector<Face> faces_2;
            vector<Edge> edges_2;
            
            // Mappa per evitare duplicati nei nuovi vertici di Classe 2
            map<tuple<double, double, double>, unsigned int> vertex_map_2;
            unsigned int next_id_2 = 0;
            vector<vertex> centroidi;
            //Inserisco i vertici  della Classe 1 in Classe 2
            for (const vertex& v : vertices_class1){
                tuple<double, double, double> key(v.x, v.y, v.z);
                if (!vertex_map_2.count(key)) {
                    vertex new_v = v;
                    new_v.id = next_id_2++;
                    vertex_map_2[key] = new_v.id;
                    ver_2.push_back(new_v);
                }
                
            }

            //Calcolo i baricentri delle facce della Classe 1
            for (const Face& faccia : faces_class1){
                const vertex& v1 = faccia.vertices[0];
                const vertex& v2 = faccia.vertices[1];
                const vertex& v3 = faccia.vertices[2];

                double centroid_x = (v1.x + v2.x + v3.x) / 3.0;
                double centroid_y = (v1.y + v2.y + v3.y) / 3.0;
                double centroid_z = (v1.z + v2.z + v3.z) / 3.0;
                vertex centroide(centroid_x,centroid_y,centroid_z, 0);
                tuple<double, double, double> key(centroid_x, centroid_y, centroid_z);
                
                if (vertex_map_2.count(key)){
                    centroide.id = vertex_map_2[key];
                }else{
                    centroide.id = next_id_2++;
                    vertex_map_2[key] = centroide.id;
                    ver_2.push_back(centroide);
                }
            
                // Aggiungo il baricentro come nuovo vertice
                centroidi.push_back(centroide);
            }
            
            //collego i baricentri ai vertici della Classe 1
            
            unsigned int id_edge_2 = 0;
            unsigned int id_face_2 = 0;
            for (size_t i=0; i < faces_class1.size();i++){
                const Face& faccia = faces_class1[i]; // prendo la faccia corrente
                const vertex& v1 = faccia.vertices[0]; //Estraggo i tre vertici della faccia triangolare
                const vertex& v2 = faccia.vertices[1];
                const vertex& v3 = faccia.vertices[2];
                const vertex& centroide = centroidi[i]; // prendo il baricentro della faccia corrente
                
                //Creo i 3 spigoli
                Edge e1(v1, centroide, id_edge_2++, 0);
                Edge e2(v2, centroide, id_edge_2++, 0);
                Edge e3(v3, centroide, id_edge_2++, 0);
                // Aggiungo gli spigoli alla Classe 2
                add_unique_edge(edges_2, e1);
                add_unique_edge(edges_2, e2);
                add_unique_edge(edges_2, e3);

            }
            
            map<unsigned int, vector<size_t>> edge_to_faces;

            // Per ogni faccia originale, associa i lati
            for (size_t id_face = 0; id_face < faces_original.size(); ++id_face) {
                const Face& face = faces_class1[id_face];

                // Considera ogni coppia di vertici della faccia
                for (size_t i = 0; i < face.vertices.size(); ++i) {
                    for (size_t j = i + 1; j < face.vertices.size(); ++j) {
                        unsigned int id1 = face.vertices[i].id;
                        unsigned int id2 = face.vertices[j].id;

                        // Trova l’edge nella triangolazione che collega questi vertici
                        for (const Edge& e : edges_class1) {
                            if ((e.origin.id == id1 && e.end.id == id2) ||
                                (e.origin.id == id2 && e.end.id == id1)) {
                                edge_to_faces[e.id].push_back(id_face);
                                break;
                            }
                        }
                    }
                }
            }
        
            for (const Edge& e : edges_class1) {
                auto it = edge_to_faces.find(e.id);
                if (it == edge_to_faces.end()) continue;

                const auto& faces_list = it->second;

                // **CASO 1: bordo della faccia del poliedro (solo 1 faccia)**
                if (faces_list.size() == 1) {
                    size_t face_id = faces_list[0];  // ID della faccia originale
                    const vertex& baricentro = centroidi[face_id];

                    // Calcolo del punto medio dello spigolo
                    double mx = (e.origin.x + e.end.x) / 2.0;
                    double my = (e.origin.y + e.end.y) / 2.0;
                    double mz = (e.origin.z + e.end.z) / 2.0;
                    vertex midpoint(mx, my, mz, 0);
                    tuple<double,double,double> key(round6(mx), round6(my), round6(mz));
                    
                    if (vertex_map_2.count(key)) {
                        midpoint.id = vertex_map_2[key];
                    } else {
                        midpoint.id = next_id_2++;
                        vertex_map_2[key] = midpoint.id;
                        ver_2.push_back(midpoint);

                    }
                    


                    // Collego il baricentro al punto medio
                    Edge eb(baricentro, midpoint, id_edge_2++, 0);
                    add_unique_edge(edges_2, eb);

                    //  Divido lo spigolo originario in due: origin–midpoint e midpoint–end
                    Edge e1(midpoint, e.origin, id_edge_2++, 0);
                    Edge e2(midpoint, e.end,    id_edge_2++, 0);
                    Edge e3(e.origin, baricentro, id_edge_2++, 0) ;
                    add_unique_edge(edges_2, e1);
                    add_unique_edge(edges_2, e2);
                    add_unique_edge(edges_2, e3);
                    vector<vertex> vertici_faccia = {baricentro, midpoint, e.origin};
                    vector<Edge> spigoli_faccia = {e1, e2, e3};
                    faces_2.push_back(Face(vertici_faccia, spigoli_faccia, id_face_2++, 0));

                }

                // CASO 2: spigolo interno alla faccia del poliedro (2 facce)
                else if (faces_list.size() == 2) {
                    const vertex& b1 = centroidi[faces_list[0]];
                    const vertex& b2 = centroidi[faces_list[1]];
                    // Scelgo un vertice dello spigolo originale
                    const vertex& v_orig = e.origin;

                    //  Collego i due baricentri tra di loro
                    Edge arco_baricentri(b1, b2, id_edge_2++, 0);
                    add_unique_edge(edges_2, arco_baricentri);

                    // Creo gli spigoli per la faccia
                    Edge arco_baricentro1_vertice(b1, v_orig, id_edge_2++, 0);
                    Edge arco_baricentro2_vertice(v_orig, b2, id_edge_2++, 0);
                    add_unique_edge(edges_2, arco_baricentro1_vertice);
                    add_unique_edge(edges_2, arco_baricentro2_vertice);
                    // Creo la faccia triangolare: (b1, v_orig, b2)
                    vector<vertex> vertici_face = {b1, v_orig, b2};
                    vector<Edge> spigoli_face = {arco_baricentri, arco_baricentro1_vertice, arco_baricentro2_vertice};
                    faces_2.push_back(Face(vertici_face, spigoli_face, id_face_2++, 0));

                }

            
            }


            faces = faces_2;        
            vertices = ver_2;
            edges = edges_2;


        
        }

	void _Polyhedron::OverAll_Triangulation(){
            
			if ((b == 0 && c > 0) || (b >0 && c == 0)) {
				_Polyhedron::Triangulation();
			}
			else if(b>0 && c>0 && b==c){
                _Polyhedron::Triangulation();
				_Polyhedron::Triangulation_2();
			}

            for (auto& v : vertices) {
                v.normalize(); // Porta ogni vertice sulla sfera di raggio 1
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
		}	
    };
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
    /*for (const auto& e : dual_edges) {
        cout << "Edge " << e.id << ": (" << e.origin.x << ", " << e.origin.y << ", " << e.origin.z
             << ") -> (" << e.end.x << ", " << e.end.y << ", " << e.end.z
             << ") length: " << e.length() << "\n";
    }*/
	//bulk delle facce.
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
	
	