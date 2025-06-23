#pragma once


#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
using namespace Eigen;
namespace PolygonalLibrary{

class vertex {
	public:
    double x, y, z;
    unsigned int id;
	int ShortPath = 0;
    unsigned int marker;
    	vertex() : x(0), y(0), z(0), marker(0) {}
    	vertex(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    	vertex(double _x, double _y, double _z, unsigned int id) : x(_x), y(_y), z(_z), id(id){}
    	vertex(double _x, double _y, double _z, unsigned int id, unsigned int _marker) : x(_x), y(_y), z(_z), id(id), marker(_marker){}
    double length() const {
        return sqrt(x * x + y * y + z * z);
    }
    void normalize(){
    	double buff = length();
    	x/=buff;
    	y/=buff;
    	z/=buff;	
	}
		
		bool operator==(const vertex& other) const {
    		const double EPS = 1e-16;
    		return (fabs(x - other.x) < EPS &&
            		fabs(y - other.y) < EPS &&
           	 		fabs(z - other.z) < EPS);
}

		
		bool operator!=(const vertex& other) const {
    		return !(*this == other);
		}

};
class Edge {
	public:
    	vertex origin;
    	vertex end;
		int ShortPath = 0;
    unsigned int id;
    unsigned int type;

    Edge() = default;
    Edge(const 	vertex& start, const 	vertex& final, unsigned int code, unsigned int _type) : origin(start), end(final), id(code), type(_type) {}
    Edge(const 	vertex& start, const 	vertex& final, unsigned int code) : origin(start), end(final), id(code) {}
    Edge(const 	vertex& start, const 	vertex& final) : origin(start), end(final) {}	
    double length() const {
        return sqrt((origin.x - end.x) * (origin.x - end.x) +
                    (origin.y - end.y) * (origin.y - end.y) +
                    (origin.z - end.z) * (origin.z - end.z));
    }
    double length_2d_xy() const{
    	        return sqrt((origin.x - end.x) * (origin.x - end.x) +
                    (origin.y - end.y) * (origin.y - end.y));
	}
	
	bool operator==(const Edge& other) const {
    return (origin == other.origin && end == other.end) ||
           (origin == other.end && end == other.origin);
		}

		
	bool operator!=(const Edge& other) const {
    		return !(*this == other);
		}
};

class Face {
	public:
    vector<Edge> edges;
    vector<	vertex> vertices;
    unsigned int id, type; 
    
   	Face()=default; 
   	
    Face(const vector<	vertex>& vec, const vector<Edge>& edg, unsigned int code, unsigned int type) {
	vertices.reserve(vec.size()) ;
	edges.reserve(edg.size());
	 vertices=vec; edges=edg;
	  id=code; 
	  this->type=type;}
	
};

Edge add_unique_edge(vector<Edge>& edges, const vertex& v1, const vertex& v2, unsigned int& counter);


void add_unique_vertex(vector<vertex>& vertices, const vertex& v, unsigned int& counter);
void add_unique_edge(vector<Edge>& , const Edge& );
void printFace(Face f) ;
Edge reverseEdge(Edge e);
 }
 