// graph.h <Starter Code>
// < Your name >
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_map>

using namespace std;



template<typename VertexT, typename WeightT>
class graph {
 private:
  unordered_map< VertexT, unordered_map<VertexT, WeightT > > Vertices;
 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() {
  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() {
    return this->Vertices.size();
  }
  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int checkEdges(unordered_map<VertexT, WeightT> edgeN) {
    int count = 0;
    for (auto& edge : edgeN) {
      count++;
    }
    return count;
  }

  int NumEdges() {
    int count = 0;
    for (auto& myPair : this->Vertices) {
      count += checkEdges(myPair.second);
    }
    return count;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool vertexExists(VertexT& v) {
    if (this->Vertices.find(v) == this->Vertices.end()) {
      return false;
    }
    return true;
  }

  bool addVertex(VertexT v) {
    // If set find returns end iterator, vertex doesn't exist yet
    // if it doesn't, then it exists and return false. Prevents duplicates
    if (vertexExists(v)) {
      return false;
    }
    // Create temp vertice to add into set. 
    this->Vertices[v];
    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // Checks if both vertexs exist in list; If either of them don't exist
    // return false
    if (!vertexExists(from) || !vertexExists(to)) {
      return false;
    }
    this->Vertices[from][to] = weight;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) {
    if (!vertexExists(from)) {
      return false;
    }
    unordered_map<VertexT, WeightT> temp = this->Vertices[from];
    if (temp.count(to) == 0) {
      return false;
    }
    weight = this->Vertices[from][to];
    return true;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) {
    set<VertexT>  S;
    if (!vertexExists(v)) {
      return S;
    }
    unordered_map<VertexT, WeightT> temp = this->Vertices[v];
    for (auto& vertexPair : temp) {
      S.insert(vertexPair.first);
    }
    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() {
    vector<VertexT> V;
    for (auto& vertex : this->Vertices) {
      V.push_back(vertex.first);
    }
    return V;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    int i = 0;
    for (auto& vertex : this->Vertices) {
      output << " " << i << ". " << vertex.first << endl;
      i++;
    }
    output << endl;
    output << "**Edges:" << endl;
    for (auto& vertex : this->Vertices) {
      output << " row " << vertex.first << ": ";
      for (auto& vertex2 : this->Vertices) {
        unordered_map<VertexT, WeightT> vertexNeighBors = vertex.second;
        if(vertexNeighBors.find(vertex2.first) == vertexNeighBors.end())
        {
          output << "F ";
        } else {
          output << "(T, "
            << vertexNeighBors[vertex2.first]
            << ") ";
        }
      }
      output << endl;
    }
    output << "**************************************************" << endl;
  }
};
