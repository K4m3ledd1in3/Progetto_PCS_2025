#include "Eigen/Eigen"
#include <iostream>
#include "Polygon.hpp"
#include <vector>

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary{

void printFace(Face f) {
    const unsigned int INVALID_ID = numeric_limits<unsigned int>::max();

    cout << "Face ID: " << f.id << endl;

    for (auto e : f.edges) {
    if(e.id!=INVALID_ID){
        cout << "Edge ID: " << e.id 
             << " | From: (" << e.origin.x << ", " << e.origin.y << ", " << e.origin.z << ")"
             << " To: (" << e.end.x << ", " << e.end.y << ", " << e.end.z << ")"
             << " Length: " << e.length() << endl;
    }else{
		cout << "Edge ID: " << e.type
             << " | From: (" << e.origin.x << ", " << e.origin.y << ", " << e.origin.z << ")"
             << " To: (" << e.end.x << ", " << e.end.y << ", " << e.end.z << ")"
             << " Length: " << e.length() << endl;    	
	}
    }
    cout << endl;
}
Edge reverseEdge(Edge e){
	return Edge(e.end, e.origin, e.id);
}
void add_unique_edge(vector<Edge>& edges, const Edge& e) {
    for (const auto& existing : edges)
        if (existing == e) return; 
    edges.push_back(e); 
    return;
}
Edge add_unique_edge(vector<Edge>& edges, const vertex& v1, const vertex& v2, unsigned int& counter) {
    Edge new_edge(v1, v2, counter);
    for (const auto& existing : edges)
        if (existing == new_edge)
            return existing; // Se esiste, lo restituisco
    edges.push_back(new_edge);
    ++counter;
    return new_edge;
}



void add_unique_vertex(vector<vertex>& vertices, const vertex& v, unsigned int& counter) {
    for (const auto& existing : vertices)
        if (existing == v) return;

    vertex new_v = v;
    new_v.id = counter++;
    vertices.push_back(new_v);
}



}