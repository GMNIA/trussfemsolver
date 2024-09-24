#pragma once
#include <vector>
#include <string>

// Node class to represent a mesh node
class Node {
public:
    int id;
    std::vector<double> coordinates;

    Node(const int& id, const std::vector<double>& coords)
        : id(id), coordinates(coords) {}
};

// Edge class to represent a connection between two nodes
class Edge {
public:
    int id;
    Node* startNode;
    Node* endNode;

    Edge(const int& id, Node* start, Node* end)
        : id(id), startNode(start), endNode(end) {}
};
