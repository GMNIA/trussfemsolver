#include "../include/Simulation.h"
#include "../include/MeshReader.h"
#include <vector>
#include <map>
#include <iostream>
#include "../lib/eigen-3.4.0/Eigen/Dense"
#include "../lib/eigen-3.4.0/Eigen/Sparse"
#include <fstream>
#include <string>
#include <array>
#include <cmath>


Simulation::Simulation(
    MeshReader& reader,
    std::vector<std::vector<int>>& blockedVectors,
    std::vector<std::vector<double>>& forceVectors,
    double E,
    std::vector<double>& AVector
) : meshReader(reader), blockedVectors(blockedVectors), forceVectors(forceVectors), E(E), AVector(AVector) {}

void Simulation::run() {
    // Initialize nodes and linear FEM elements
    std::list<Node*>* pNodes = &meshReader.nodes;
    std::list<Edge*>* pEdges = &meshReader.edges;
    std::list<Node*> nodes = *pNodes;
    std::list<Edge*> edges = *pEdges;

    // Initialise variables to check if 3d or 2d simulation
    bool allThree = true;
    bool allTwo = true;

    // Check if 3d or 2d simulation if all nodes have 2 or 3 degrees of freedom (dof)
    for (Node* node : nodes) {
        if (node->coordinates.size() == 3) {
            allTwo = false;
        } else if (node->coordinates.size() == 2) {
            allThree = false;
        } else {
            allTwo = false;
            allThree = false;
            break;
        }
    }

    // Create a map of the degrees of freedom (dof) to their id
    std::map<int, std::vector<int>> mapDof;
    for (Node* node : nodes) {
        if (allThree) {
            mapDof[node->id] = { 3 * node->id - 2, 3 * node->id - 1,  3 * node->id};
        } else if (allTwo) {
            mapDof[node->id] = { 2 * node->id - 1, 2 * node->id};
        }
    }    

    // Discriminate between 2 dof case (2D space example) and 3 dof case (3D-space example)
    int inputDofEdge;
    int inputDofNode;
    if (allThree) {
        inputDofEdge = 6;
        inputDofNode = 3;
    } else if (allTwo) {
        inputDofEdge = 4;
        inputDofNode = 2;
    }
    const int dofEdge = inputDofEdge;
    const int dofNode = inputDofNode;

    // Initialise structured data for stiffness
    std::list<std::vector<std::vector<double>>> Kels;
    std::vector<std::vector<double>> Q(dofEdge, std::vector<double>(dofEdge));
    std::vector<std::vector<double>> Q_T(dofEdge, std::vector<double>(dofEdge));
    std::vector<std::vector<double>> Kel(dofEdge, std::vector<double>(dofEdge));

    // Temporary variable helping calculation
    double x, y, l;
    double z = 0;
    double s = 0;

    // Input materials
    double memberProp;
    std::vector<double> startV;
    std::vector<double> endV;
    int j = 0;
    for (Edge* edge : edges) {
        if (edge != nullptr && edge->startNode != nullptr && edge->endNode != nullptr) {
            // Read coordinates of each node and calculate geometry parameters
            startV = edge->startNode->coordinates;
            endV = edge->endNode->coordinates;
            x = endV[0] - startV[0];
            y = endV[1] - startV[1];
            if (allThree) {
                z = endV[2] - startV[2];
            } else if (allTwo) {
                z = 0;
            }

            // If z not existing it is zero
            l = sqrt(pow(x, 2) + pow(y, 2) + pow(z,2));
            if (allThree) {
                // Global stiffness matrix in 3d
                s = E * AVector[j] / pow(l, 3);
                Kel = {
                    {s*x*x, s*x*y, s*x*z, -s*x*x, -s*x*y, -s*x*z},
                    {s*x*y, s*y*y, s*y*z, -s*x*y, -s*y*y, -s*y*z},
                    {s*x*z, s*y*z, s*z*z, -s*x*z, -s*y*z, -s*z*z},
                    {-s*x*x, -s*x*y, -s*x*z, s*x*x, s*x*y, s*x*z},
                    {-s*x*y, -s*y*y, -s*y*z, s*x*y, s*y*y, s*y*z},
                    {-s*x*z, -s*y*z, -s*z*z, s*x*z, s*y*z, s*z*z}
                };
                Kels.push_back(Kel);
            } else if (allTwo) {
                // Rotation matrix in 2d
                Q = {{
                    {x/l, y/l, 0, 0},
                    {-y/l, x/l, 0, 0}, 
                    {0, 0, x/l, y/l},
                    {0, 0, -y/l, x/l}
                }};

                // Local stiffness matrix in 2d
                Kel = {{
                    {1,  0, -1, 0},
                    {0,  0, 0,  0},
                    {-1, 0, 1,  0},
                    {0,  0, 0,  0},
                }};

                // Calculate the transpose of Q
                for (int i = 0; i < dofEdge; ++i) {
                    for (int j = 0; j < dofEdge; ++j) {
                        Q_T[j][i] = Q[i][j];
                    }
                }

                // Calculate and apply linear element properties: in this case member is 1 linear element
                memberProp = E * AVector[j] / l;
                for (int i = 0; i < dofEdge; ++i) {
                    for (int j = 0; j < dofEdge; ++j) {
                        Kel[i][j] = Kel[i][j] * memberProp;
                    }
                }

                // Calculate element stiffness matrix in global coordinates
                Kel = multiplyM(Q_T, Kel);
                Kel = multiplyM(Kel, Q);
                Kels.push_back(Kel);
            }

        }
        j++;
    }

    // Id matrix to define global model connectivity
    const int KelRows = dofEdge;
    const int nEdges = edges.size();
    int idMatrix[KelRows][nEdges];
    j = 0;
    for (Edge* edge : edges) {
        if (edge != nullptr && edge->startNode != nullptr && edge->endNode != nullptr) {
            auto startDof = mapDof.find(edge->startNode->id)->second;
            auto endDof = mapDof.find(edge->endNode->id)->second;
            // Explanation of id matrix: the columns are the edge (linear elements), the rows show the 
            if (allThree) {
                idMatrix[0][j] = startDof[0];
                idMatrix[1][j] = startDof[1];
                idMatrix[2][j] = startDof[2];
                idMatrix[3][j] = endDof[0];
                idMatrix[4][j] = endDof[1];
                idMatrix[5][j] = endDof[2];
            } else if (allTwo) {
                idMatrix[0][j] = startDof[0];
                idMatrix[1][j] = startDof[1];
                idMatrix[2][j] = endDof[0];
                idMatrix[3][j] = endDof[1];
            }
            j++;
        }
    }

    // Build global stiffness matrix using id matrix and local stiffness matrices in global coordinates
    Eigen::MatrixXd Kglobal(dofNode * nodes.size(), dofNode * nodes.size());
    int k = 0;
    for (const auto& Kel : Kels) {
        for (int i = 0; i < dofEdge; i++) {
            for (int j = 0; j < dofEdge; j++) {
                int globalRow = idMatrix[i][k] - 1;
                int globalCol = idMatrix[j][k] - 1;
                Kglobal(globalRow, globalCol) += Kel[i][j];
            }
        }
        k += 1;
    }

    // Initialise displacement d and Force F vectors
    Eigen::VectorXd d = Eigen::VectorXd::Zero(dofNode * nodes.size());
    Eigen::VectorXd F = Eigen::VectorXd::Zero(dofNode * nodes.size());

    // Impose boundary condtions (blocked nodes)
    std::map<int, double> imposed;
    std::vector<int> blockedNodes = meshReader.getboundaryNodeIds();
     if (allThree) {
        for (int i = 0; i < blockedNodes.size(); ++i) {
            if (blockedVectors[i][0] == 0) {
                imposed[blockedNodes[i] * dofNode - 3] = 0.0;
            }
            if (blockedVectors[i][1] == 0) {
                imposed[blockedNodes[i] * dofNode - 2] = 0.0;
            }
            if (blockedVectors[i][2] == 0) {
                imposed[blockedNodes[i] * dofNode - 1] = 0.0;
            }
        }
     } else if (allTwo) {
        for (int i = 0; i < blockedNodes.size(); ++i) {
            if (blockedVectors[i][0] == 0) {
                imposed[blockedNodes[i] * dofNode - 2] = 0.0;
            }
            if (blockedVectors[i][1] == 0) {
                imposed[blockedNodes[i] * dofNode - 1] = 0.0;
            }
        }
     }
    

    // Apply the input
    std::map<int, double> imposeF;
    std::vector<int> forceNodes = meshReader.getForceNodeIds();
    for (int i = 0; i < forceNodes.size(); ++i) {
        // Apply force Vectors which must match forceNodes size
        for (int j = 0; j < dofNode; ++j) {
            imposeF[(forceNodes[i] - 1) * dofNode + j] = forceVectors[i][j];
        }
    }

    // Apply imposed dof in force and deformations
    for (const auto& pair : imposed) {
        d(pair.first) = pair.second;
    }
    for (const auto& pair : imposeF) {
        F(pair.first) = pair.second;
    }

    // Collect indices for equationreduction
    std::vector<int> freeIndices;
    for (int i = 0; i < dofNode * nodes.size(); ++i) {
        if (imposed.find(i) == imposed.end()) {
            freeIndices.push_back(i);
        }
    }

    // Build reduced equation system
    int reducedSize = freeIndices.size();
    Eigen::MatrixXd Kglobal_red(reducedSize, reducedSize);
    Eigen::VectorXd F_red(reducedSize);
    for (int i = 0; i < reducedSize; ++i) {
        for (int j = 0; j < reducedSize; ++j) {
            Kglobal_red(i, j) = Kglobal(freeIndices[i], freeIndices[j]);
        }
        F_red(i) = F(freeIndices[i]);
    }

    // Solve the reduced system
    Eigen::VectorXd d_red = Kglobal_red.colPivHouseholderQr().solve(F_red);
    for (int i = 0; i < d_red.size(); ++i) {
        d(freeIndices[i]) = d_red(i);
    }
    
    // Recalculate F using the original Kglobal and solved d
    F = Kglobal * d;

    // Store results in the internal class vectors
    deformations.resize(d.size());
    internalForces.resize(F.size());

    for (int i = 0; i < d.size(); ++i) {
        deformations[i] = d(i);
    }
    for (int i = 0; i < F.size(); ++i) {
        internalForces[i] = F(i);
    }
}

// Getter for internal forces
const std::vector<double>& Simulation::getInternalForces() const {
    return internalForces;
}

// Getter for deformations
const std::vector<double>& Simulation::getDeformations() const {
    return deformations;
}

// Custom function to multiply two matrices
std::vector<std::vector<double>> Simulation::multiplyM(
    const std::vector<std::vector<double>>& mat1,
    const std::vector<std::vector<double>>& mat2
    ) {
    // Assume mat1 and mat2 are square matrices of the same size
    int size = mat1.size();
    std::vector<std::vector<double>> result(size, std::vector<double>(size, 0.0));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}
