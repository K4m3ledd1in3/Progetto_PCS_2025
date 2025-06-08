#include "Eigen/Eigen"
#include <iostream>
#include "Polygon.hpp"
#include <vector>
using namespace std;
using namespace Eigen;


namespace PolygonalLibrary{

void printFace(Face f) {
    cout << "Face ID: " << f.id << endl;

    for (auto e : f.edges) {
    if(e.id!=4294967295){
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
void add_unique_vertex(vector<vertex>& vertices, const vertex& v, int& counter) {

    for (const auto& existing : vertices)
        if (existing == v) return;
    vertices.push_back(v);
        cout << counter++ << endl;

    return;
}

}