
#include <iostream>
#include "UCDUtilities.hpp"
#include "Polyhedron.hpp"
#include "Polygon.hpp"
#include <unordered_set>

using namespace Eigen;
using namespace PolygonalLibrary;
using namespace PolyhedronLibrary;
using namespace std;
using namespace Gedim;

int main() {
    unsigned int p, q, b, c;
    cout << "Inserisci p, q, b, c : ";
    cin >> p >> q >> b >> c;
    _Polyhedron pp;
    // Costruisci il poliedro in base a p, q, b, c
    if (p == 3 && q == 3) {
        Tetrahedron t;
        pp = _Polyhedron(t.faces, t.edge, t.vertices, b, c, p, q);
    } else if (p == 4 && q == 3) {
        Cube C;
        pp = _Polyhedron(C.faces, C.edge, C.vertices, b, c, p, q);
    } else if (p == 3 && q == 4) {
        Octahedron o;
        pp = _Polyhedron(o.faces, o.edge, o.vertices, b, c, p, q);
    } else if (p == 5 && q == 3) {
        Dodecahedron d;
        pp = _Polyhedron(d.faces, d.edge, d.vertices, b, c, p, q);
    } else if (p == 3 && q == 5) {
        Icosahedron i;
        pp = _Polyhedron(i.faces, i.edge, i.vertices, b, c, p, q);
    } else {
        cout << "Parametri non validi. Termino.\n";
        return 1;
    }

    pp.OverAll_Triangulation();
    pp.WriteTXT(); // Scrive i dati in file TXT
    pp.printFaces();  // Stampa facce generate (facoltativo)
    pp.GenerateDual(); // Genera duale (facoltativo)
    
    // Mostra i vertici generati (per scegliere ID)
    cout << "\nVertici generati:\n";
    for (auto v : pp.vertices) {
        cout << "ID: " << v.id
             << " - (" << v.x << ", " << v.y << ", " << v.z << ")\n";
    }

    // Input ID vertici per il cammino minimo
    unsigned int id_start, id_end;
    cout << "Inserisci id vertice di partenza: ";
    cin >> id_start;
    cout << "Inserisci id vertice di arrivo: ";
    cin >> id_end;

    // Calcola cammino minimo
    vector<vertex> path = ShortestPath(pp.vertices, pp.edges, id_start, id_end);

    if (!path.empty()) {
        cout << "Cammino minimo trovato:\n";
        for (size_t i = 0; i < path.size(); i++) {
            cout << "ID: " << path[i].id
                 << " (" << path[i].x << ", " << path[i].y << ", " << path[i].z << ")\n";
        }
    } else {
        cout << "Nessun cammino trovato.\n";
    }

    vector<double> punti(pp.vertices.size(), 0.0);
    vector<double> segmenti(pp.edges.size(), 0.0);

    for (const auto& v:path)
    {   
        punti[v.id]=1.0;
        
    }
   

    for(unsigned int i = 0;i < path.size() -1;i++)
    {
        unsigned int a = path[i].id;
        unsigned int b = path[i+1].id;
        for (const auto& arco: pp.edges){
            unsigned int c = arco.id;
            unsigned int o = arco.origin.id;
            unsigned int e = arco.end.id;
            unordered_set<unsigned int> insieme_controllo = {o, e};
            if ((insieme_controllo.count(a)+insieme_controllo.count(b))==2)
            {
                segmenti[c]=1.0;
                break;
            }

        }
    }

    


    // Export completo (UCD e txt)
    UCDUtilities utilities;
    vector<Gedim::UCDProperty<double>> points_properties;
    vector<Gedim::UCDProperty<double>> segments_properties;
    Gedim::UCDProperty<double> point_proprietà;
    Gedim::UCDProperty<double> segmenti_proprietà;
    point_proprietà.Label = "ShortPath";
    point_proprietà.NumComponents = 1;
    point_proprietà.Data = punti.data(); // data = puntatore ai dati del vettore
    points_properties.push_back(point_proprietà);
    segmenti_proprietà.Label = "ShortPath";
    segmenti_proprietà.NumComponents = 1;
    segmenti_proprietà.Data = segmenti.data(); // data = puntatore ai dati del vettore
    segments_properties.push_back(segmenti_proprietà);


   

    // 🔸 Export UCD
    utilities.ExportPoints("./Cell0Ds.inp", ConvertVerticesToEigen(pp.vertices), points_properties);
    utilities.ExportSegments("./Cell1Ds.inp", ConvertVerticesToEigen(pp.vertices), ConvertEdgesToEigen(pp.edges), points_properties, segments_properties);
    //utilities.ExportPolygons("./Cell2Ds.inp", P.cell0Ds_coordinates, P.cell2Ds_vertices, points_properties, polygons_properties, materials);

    
    return 0;
}
