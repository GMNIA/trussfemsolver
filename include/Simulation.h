#pragma once
#include "MeshReader.h"
#include <vector>
#include <map>
#include <iostream>

class Simulation {
public:
    // Constructor
    Simulation(
        MeshReader& reader,
        std::vector<std::vector<int>>& blockedVectors,
        std::vector<std::vector<double>>& forceVectors,
        double E,
        std::vector<double>& AVector
    );

    // Method to run the simulation
    void run();

    // Getter methods for results
    const std::vector<double>& getInternalForces() const;
    const std::vector<double>& getDeformations() const;

    // Custom multiplication method to test basic functions instead of using directly Eigen
    std::vector<std::vector<double>> multiplyM(
        const std::vector<std::vector<double>>& mat1,
        const std::vector<std::vector<double>>& mat2
    );

private:
    // Member variables
    MeshReader& meshReader;  // Reference to the MeshReader object
    std::vector<std::vector<int>> blockedVectors;  // List of blocked vectors
    std::vector<std::vector<double>> forceVectors; // List of force vectors
    double E;  // Modulus of elasticity
    std::vector<double> AVector;  // Cross-sectional areas
    std::vector<double> internalForces;  // To store F (internal forces)
    std::vector<double> deformations;    // To store d (deformations)
};