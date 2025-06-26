#include "gtest/gtest.h"
#include "Polyhedron.hpp"

using namespace PolyhedronLibrary;

//---------------------------------------------------------------------------
//TEST TRIANGOLAZIONE
//Triangolazione dell'icosaedro
TEST(PolyhedronTest, IcosahedronTriangulation) {
    // Parametri per l'icosaedro (p=3, q=5)
    unsigned int p_val = 3;
    unsigned int q_val = 5;
    unsigned int b_val = 8;  // Puoi cambiare per b=1,2,3...
    unsigned int c_val = 0;

    // Crea l'icosaedro
    Icosahedron i;
    _Polyhedron poly(i.faces, i.edge, i.vertices, b_val, c_val, p_val, q_val);
    
    // Esegui la triangolazione
    poly.OverAll_Triangulation();

    // Conta vertici, lati e facce
    unsigned int T = b_val * b_val + b_val * c_val + c_val * c_val;
    unsigned int expected_vertices = 10 * T + 2;
    unsigned int expected_edges = 30 * T;
    unsigned int expected_faces = 20 * T;

    unsigned int actual_vertices = poly.vertices.size();
    unsigned int actual_edges = poly.edges.size();
    unsigned int actual_faces = poly.faces.size();

    // Stampa per debug
    std::cout << "Icosaedro con b=" << b_val << ": V=" << actual_vertices << ", E=" << actual_edges << ", F=" << actual_faces << "\n";

    // Verifica i conti teorici
    ASSERT_EQ(actual_vertices, expected_vertices) << "Numero di vertici errato.";
    ASSERT_EQ(actual_edges, expected_edges) << "Numero di lati errato.";
    ASSERT_EQ(actual_faces, expected_faces) << "Numero di facce errato.";

    // Formula di Eulero
    int V = static_cast<int>(actual_vertices);
    int E = static_cast<int>(actual_edges);
    int F = static_cast<int>(actual_faces);
    ASSERT_EQ(V - E + F, 2) << "Formula di Eulero non soddisfatta.";
}



//Triangolazione del tetraedro
TEST(PolyhedronTest, TetrahedronTriangulation) {
    unsigned int p_val = 3;
    unsigned int q_val = 3;
    unsigned int b_val = 2;  // Cambia b se vuoi test più dettagliati
    unsigned int c_val = 0;

    Tetrahedron t;
    _Polyhedron poly(t.faces, t.edge, t.vertices, b_val, c_val, p_val, q_val);
    poly.OverAll_Triangulation();

    unsigned int T = b_val * b_val + b_val * c_val + c_val * c_val;
    unsigned int expected_vertices = 2 * T + 2;
    unsigned int expected_edges = 6 * T;
    unsigned int expected_faces = 4 * T;

    unsigned int actual_vertices = poly.vertices.size();
    unsigned int actual_edges = poly.edges.size();
    unsigned int actual_faces = poly.faces.size();

    std::cout << "Tetraedro (b=" << b_val << "): V=" << actual_vertices
              << ", E=" << actual_edges << ", F=" << actual_faces << "\n";

    ASSERT_EQ(actual_vertices, expected_vertices) << "Numero di vertici errato.";
    ASSERT_EQ(actual_edges, expected_edges) << "Numero di lati errato.";
    ASSERT_EQ(actual_faces, expected_faces) << "Numero di facce errato.";

    // Formula di Eulero
    int V = static_cast<int>(actual_vertices);
    int E = static_cast<int>(actual_edges);
    int F = static_cast<int>(actual_faces);
    ASSERT_EQ(V - E + F, 2) << "Formula di Eulero non soddisfatta.";
}

// Triangolazione dell'ottaedro

TEST(PolyhedronTest, OctahedronTriangulation) {
    unsigned int p_val = 3;
    unsigned int q_val = 4;
    unsigned int b_val = 6;  // Cambia b se vuoi test più dettagliati
    unsigned int c_val = 0;

    Octahedron o;
    _Polyhedron poly(o.faces, o.edge, o.vertices, b_val, c_val, p_val, q_val);
    poly.OverAll_Triangulation();

    unsigned int T = b_val * b_val + b_val * c_val + c_val * c_val;
    unsigned int expected_vertices = 4 * T + 2;
    unsigned int expected_edges = 12 * T;
    unsigned int expected_faces = 8 * T;

    unsigned int actual_vertices = poly.vertices.size();
    unsigned int actual_edges = poly.edges.size();
    unsigned int actual_faces = poly.faces.size();

    std::cout << "Ottaedro (b=" << b_val << "): V=" << actual_vertices
              << ", E=" << actual_edges << ", F=" << actual_faces << "\n";

    ASSERT_EQ(actual_vertices, expected_vertices) << "Numero di vertici errato.";
    ASSERT_EQ(actual_edges, expected_edges) << "Numero di lati errato.";
    ASSERT_EQ(actual_faces, expected_faces) << "Numero di facce errato.";

    // Formula di Eulero
    int V = static_cast<int>(actual_vertices);
    int E = static_cast<int>(actual_edges);
    int F = static_cast<int>(actual_faces);
    ASSERT_EQ(V - E + F, 2) << "Formula di Eulero non soddisfatta.";
}


//------------------------

//TEST CAMAMINO MINIMO

//TEST cammino minimo nell'icosaedro
TEST(PolyhedronTest, IcosahedronShortestPath) {
    // Parametri per l'icosaedro (p=3, q=5)
    unsigned int p_val = 3;
    unsigned int q_val = 5;
    unsigned int b_val = 2;  
    unsigned int c_val = 0;

    // Crea l'icosaedro
    Icosahedron i;
    _Polyhedron p(i.faces, i.edge, i.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Seleziona due vertici (0 e uno valido, ad esempio 8)
    unsigned int id1 = 0;
    unsigned int id2 = 8;

    // Controlla che id2 sia valido
    ASSERT_LT(id2, p.vertices.size()) << "id2 fuori intervallo.";

    // Calcola il cammino minimo
    auto path = ShortestPath(p.vertices, p.edges, id1, id2);

    

    // Verifica che il cammino esista e sia valido
    ASSERT_FALSE(path.empty()) << "Nessun cammino trovato.";
    ASSERT_GE(path.size(), 2) << "Cammino troppo corto.";
    ASSERT_EQ(path.front().id, id1) << "Il cammino non parte dal vertice corretto.";
    ASSERT_EQ(path.back().id, id2) << "Il cammino non termina nel vertice corretto.";
}

// Test per il cammino minimo nelL'OTTAEDRO
TEST(PolyhedronTest, OctahedronShortestPath) {
    
    unsigned int p_val = 3;
    unsigned int q_val = 4;
    unsigned int b_val = 2;  // Puoi scegliere 1,2,3...
    unsigned int c_val = 0;

    //Crea l'ottaedro
    Octahedron o;
    _Polyhedron p(o.faces, o.edge, o.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Seleziona due vertici (0 e uno valido, ad esempio 8)
    unsigned int id1 = 0;
    unsigned int id2 = 8;

    // Controlla che id2 sia valido
    ASSERT_LT(id2, p.vertices.size()) << "id2 fuori intervallo.";

    // Calcola il cammino minimo
    auto path = ShortestPath(p.vertices, p.edges, id1, id2);

    

    // Verifica che il cammino esista e sia valido
    ASSERT_FALSE(path.empty()) << "Nessun cammino trovato.";
    ASSERT_GE(path.size(), 2) << "Cammino troppo corto.";
    ASSERT_EQ(path.front().id, id1) << "Il cammino non parte dal vertice corretto.";
    ASSERT_EQ(path.back().id, id2) << "Il cammino non termina nel vertice corretto.";
}

//TEST cammino minimo nel tetraedro
TEST(PolyhedronTest, TetrahedronShortestPath) {
    
    unsigned int p_val = 3;
    unsigned int q_val = 3;
    unsigned int b_val = 2;  // Puoi scegliere 1,2,3...
    unsigned int c_val = 0;

    // Cera il tetraedro
    Tetrahedron t;
    _Polyhedron p(t.faces, t.edge, t.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Seleziona due vertici (0 e uno valido, ad esempio 8)
    unsigned int id1 = 0;
    unsigned int id2 = 4;

    // Controlla che id2 sia valido
    ASSERT_LT(id2, p.vertices.size()) << "id2 fuori intervallo.";

    // Calcola il cammino minimo
    auto path = ShortestPath(p.vertices, p.edges, id1, id2);

    

    // Verifica che il cammino esista e sia valido
    ASSERT_FALSE(path.empty()) << "Nessun cammino trovato.";
    ASSERT_GE(path.size(), 2) << "Cammino troppo corto.";
    ASSERT_EQ(path.front().id, id1) << "Il cammino non parte dal vertice corretto.";
    ASSERT_EQ(path.back().id, id2) << "Il cammino non termina nel vertice corretto.";
}

//--------------------------
//TEST DUALE
// Test per duale di icosaedro


TEST(PolyhedronTest, GenerateDualTest1) {
    // Parametri per il poliedro base: icosaedro (p=3, q=5), b=2
    unsigned int p_val = 3;
    unsigned int q_val = 5;
    unsigned int b_val = 2;
    unsigned int c_val = 0;

    // Crea l'icosaedro
    Icosahedron i;
    _Polyhedron p(i.faces, i.edge, i.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Genera il duale
    p.GenerateDual();  // Assicurati che questa funzione generi e mostri informazioni

    // Controlla che ci siano centroidi generati (uno per ogni faccia)
    ASSERT_EQ(p.faces.size(), p.faces.size()) << "Numero di facce e centroidi del duale non corrisponde.";

    
}

// Test per duale di tetraedro
TEST(PolyhedronTest, GenerateDualTest2) {
    // Parametri per il poliedro base: icosaedro (p=3, q=5), b=2
    unsigned int p_val = 3;
    unsigned int q_val = 5;
    unsigned int b_val = 2;
    unsigned int c_val = 0;

    // Crea tetraedro
    Tetrahedron t;
    _Polyhedron p(t.faces, t.edge, t.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Genera il duale
    p.GenerateDual();  // Assicurati che questa funzione generi e mostri informazioni

    // Controlla che ci siano centroidi generati (uno per ogni faccia)
    ASSERT_EQ(p.faces.size(), p.faces.size()) << "Numero di facce e centroidi del duale non corrisponde.";

    
}

// Test per duale di ottaedro
TEST(PolyhedronTest, GenerateDualTest3) {
    // Parametri per il poliedro base: icosaedro (p=3, q=5), b=2
    unsigned int p_val = 3;
    unsigned int q_val = 4;
    unsigned int b_val = 2;
    unsigned int c_val = 0;

    // Crea l'ottaedro
    Octahedron o;
    
    _Polyhedron p(o.faces, o.edge, o.vertices, b_val, c_val, p_val, q_val);
    p.OverAll_Triangulation();

    // Genera il duale
    p.GenerateDual();  // Assicurati che questa funzione generi e mostri informazioni

    // Controlla che ci siano centroidi generati (uno per ogni faccia)
    ASSERT_EQ(p.faces.size(), p.faces.size()) << "Numero di facce e centroidi del duale non corrisponde.";

    
  
}



TEST(PolyhedronTest, TetrahedronClassIITriangulationExactFormula) {
    unsigned int p_val = 3;
    unsigned int q_val = 3;
    unsigned int b_val = 2;
    unsigned int c_val = 2;  // CLASSE II

    Tetrahedron tetra;
    _Polyhedron poly(tetra.faces, tetra.edge, tetra.vertices, b_val, c_val, p_val, q_val);
    poly.OverAll_Triangulation;

    // Dati di partenza
    unsigned int numF = 4; // Facce del tetraedro
    unsigned int numV = 4; // Vertici del tetraedro
    unsigned int numE = 6; // Spigoli del tetraedro 
    unsigned int b = b_val;

    // Numero atteso di triangolini per faccia
    unsigned int expected_faces = numF *((3*b* b)+3*b); // Ogni triangolo viene diviso in b^2 triangolini
    unsigned int expected_vertices = numV + numE * (2 * b - 1) + numF * ((3 * b * b) / 2 - (3 * b) / 2 + 1);
    unsigned int expected_edges = numE * (2 * b) + numF * ((9 * b * b) / 2 + (3 * b) / 2);
    unsigned int actual_vertices = poly.vertices.size();
    unsigned int actual_edges = poly.edges.size();
    unsigned int actual_faces = poly.faces.size();


    //  Controllo formula di Eulero
    int V = static_cast<int>(actual_vertices);
    int E = static_cast<int>(actual_edges);
    int F = static_cast<int>(actual_faces);
    ASSERT_EQ(V - E + F, 2) << "Formula di Eulero fallita.";
}

