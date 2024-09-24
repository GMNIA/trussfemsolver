#include "../include/MeshReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <algorithm>

MeshReader::MeshReader(const std::string &filename) : filename(filename) {}

MeshReader::~MeshReader() {
    for (Node* node : nodes) {
        delete node;  // Delete each node
    }
    for (Edge* edge : edges) {
        delete edge;  // Delete each edge
    }
}

void MeshReader::readFile() {
    std::ifstream file(this->filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << this->filename << std::endl;
        return;
    }

    // Parse nodes
    file.seekg(0, std::ios::beg);
    parseNodes(file);
    file.clear();

    // Parse linear elements (edges)
    file.seekg(0, std::ios::beg);
    parseEdges(file);
    file.clear();

    // Parse group of nodes
    file.seekg(0, std::ios::beg);
    parseGroups(file, "GROUP_NO");
    file.clear();

    // Parse group of elements
    file.seekg(0, std::ios::beg);
    parseGroups(file, "GROUP_MA");
    file.clear();
}

void MeshReader::parseNodes(std::ifstream& file) {
    std::string line;
    bool readData = false;
    while (getline(file, line)) {
        if (line == "COOR_3D") {
            readData = true;
            continue;
        }
        if (readData && line == "FINSF") break;
        
        if (readData) {
            std::istringstream iss(line);
            std::string idString;
            double x, y, z;
            if (iss >> idString >> x >> y) {
                int id = extractInteger(idString);
                if (iss >> z) {
                    // Three values present (x, y, z)
                    Node* node = new Node(id, {x, y, z});
                    nodes.push_back(node);
                    nodeMap[id] = node;
                } else {
                    // If only two coordinates are present, save only two and forget about z
                    Node* node = new Node(id, {x, y});
                    nodes.push_back(node);
                    nodeMap[id] = node;
                }
            } else {
                std::cerr << "Failed to parse line: " << line << std::endl;
            }
        }
    }
}

void MeshReader::parseEdges(std::ifstream& file) {
    std::string line;
    bool readData = false;

    while (getline(file, line)) {
        if (line == "SEG2") {
            readData = true;
            continue;
        }
        if (readData && line == "FINSF") break;
        
        if (readData) {
            std::istringstream iss(line);
            std::string idString, startIdString, endIdString;
            if (iss >> idString >> startIdString >> endIdString) {
                int id = extractInteger(idString);
                int startId = extractInteger(startIdString);
                int endId = extractInteger(endIdString);
                if (nodeMap.find(startId) != nodeMap.end() && nodeMap.find(endId) != nodeMap.end()) {
                    Node* startNode = nodeMap[startId];
                    Node* endNode = nodeMap[endId];
                    Edge* edge = new Edge(id, startNode, endNode);
                    edges.push_back(edge);
                    edgeMap[id] = edge;
                } else {
                    std::cerr << "Node ID not found for edge: " << id << std::endl;
                    std::cerr << "Line " << line << std::endl;
                    std::cerr << "Read element " << idString << " " <<
                    startIdString << " "<< endIdString << std::endl;
                }
            } else {
                std::cerr << "Failed to parse line: " << line << std::endl;
            }
        }
    }
}

int MeshReader::extractInteger(const std::string& input) {
    std::istringstream stream(input);
    std::string token;
    int number;

    while (stream >> token) {
        if (token[0] == 'N' || token[0] == 'B') {
            number = std::stoi(token.substr(1));
        }
    }
    return number;
}



void MeshReader::parseGroups(std::ifstream& file, const std::string& groupKeyWord) {
    std::string line;
    std::string groupName;
    std::list<Node*> nodes;
    std::list<Edge*> edges;

    while (getline(file, line)) {
        // Start processing when "GROUP_NO" is encountered
        if (line == groupKeyWord) {
            // If we've already read a group, save it before starting a new one
            if (!groupName.empty()) {
                if (groupKeyWord == "GROUP_NO") {
                    nodes.clear();  // Clear the list for the new group
                } else if (groupKeyWord == "GROUP_MA") {
                    edges.clear();  // Clear the list for the new group
                }
            }
            
            // Read the group name
            getline(file, groupName);
            continue;
        }

        // Stop processing when "FINSF" is encountered and save the current group
        if (line == "FINSF") {
            if (!groupName.empty()) {
                if (groupKeyWord == "GROUP_NO") {
                    nodeGroups[groupName] = nodes;  // Store the nodes for the current group
                    nodes.clear();  // Clear nodes for any future groups
                    groupName.clear();  // Clear the group name for future use
                } else if (groupKeyWord == "GROUP_MA") {
                    edgeGroups[groupName] = edges;  // Store the nodes for the current group
                    edges.clear();  // Clear nodes for any future groups
                    groupName.clear();  // Clear the group name for future use
                }
            }
            continue;
        }
        
        // If we are reading data, extract integers from the line
        if (!groupName.empty()) {  // Only process lines if we're inside a group
            std::istringstream stream(line);
            std::string token;
            while (stream >> token) {
                if (groupKeyWord == "GROUP_NO") {
                    int nodeId = extractInteger(token);
                    if (nodeMap.find(nodeId) != nodeMap.end()) {
                        Node* node = nodeMap[nodeId];
                        nodes.push_back(node);
                    }
                } else if (groupKeyWord == "GROUP_MA") {
                    int edgeId = extractInteger(token);
                    if (edgeMap.find(edgeId) != edgeMap.end()) {
                        Edge* edge = edgeMap[edgeId];
                        edges.push_back(edge);
                    }
                }
            }
        }   
    }
}


std::vector<int> MeshReader::getboundaryNodeIds() {
    std::vector<int> blockedNodeIds;

    // Iterate through each group in nodeGroups
    for (const auto& groupPair : nodeGroups) {
        std::string groupName = groupPair.first;
        const std::list<Node*>& nodes = groupPair.second;

        // Convert groupName to lowercase
        std::transform(groupName.begin(), groupName.end(), groupName.begin(), ::tolower);

        // Check if "boundary" or "constraint" is in the group name (case-insensitive)
        if (groupName.find("boundary") != std::string::npos
            || groupName.find("constraint") != std::string::npos
            || groupName.find("blocked") != std::string::npos
        ) {
            for (const Node* node : nodes) {
                blockedNodeIds.push_back(node->id);
            }
        }
    }
    return blockedNodeIds;
}

std::vector<int> MeshReader::getForceNodeIds() {
    std::vector<int> forceNodeIds;

    // Iterate through each group in nodeGroups
    for (const auto& groupPair : nodeGroups) {
        std::string groupName = groupPair.first;
        const std::list<Node*>& nodes = groupPair.second;

        // Convert groupName to lowercase
        std::transform(groupName.begin(), groupName.end(), groupName.begin(), ::tolower);

        // Check if "boundary" or "constraint" is in the group name (case-insensitive)
        if (groupName.find("force") != std::string::npos) {
            for (const Node* node : nodes) {
                forceNodeIds.push_back(node->id);
            }
        }
    }
    return forceNodeIds;
}