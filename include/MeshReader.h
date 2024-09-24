#pragma once
#include "MeshStructure.h"
#include <string>
#include <list>
#include <map>
#include <fstream>


class MeshReader {
public:
    std::list<Node*> nodes;  // List to store nodes
    std::map<int, Node*> nodeMap;  // Map to easily access nodes by their ID
    std::list<Edge*> edges;  // List to store edges
    std::map<int, Edge*> edgeMap;  // Map to easily access edges by their ID
    std::map<std::string, std::list<Edge*>> edgeGroups;
    std::map<std::string, std::list<Node*>> nodeGroups;

    MeshReader(const std::string &filename);
    ~MeshReader(); // Destructor to handle cleanup
    void readFile();
    void parseNodes(std::ifstream& file);
    void parseEdges(std::ifstream& file);
    void parseGroups(std::ifstream& file, const std::string& groupKeyWord);
    int extractInteger(const std::string& input);
    std::vector<int> getboundaryNodeIds();
    std::vector<int> getForceNodeIds();
private:
    std::string filename;
};
