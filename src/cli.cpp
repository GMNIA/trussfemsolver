#include "../include/Simulation.h"
#include "../include/MeshReader.h"
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    // Default file path to the example mail file
    std::string filePath = "./testfiles/example.mail";

    // If the user provides a custom file path, use it
    if (argc > 1) {
        filePath = argv[1];
    }

    std::cout << "Running simulation using file: " << filePath << std::endl;

    // Define material properties and cross-section areas
    double E = 2.1E8;  // Young's modulus for steel
    std::vector<double> AVector = {0.008770};  // Cross-sectional area

    // Define boundary conditions (blocked nodes) and force vectors
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0}  // Example: Fully fixed node
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, 10, 0}  // Example: Applied force
    };

    // Initialize MeshReader and Simulation
    MeshReader meshReader(filePath);
    meshReader.readFile();
    
    Simulation simulation(meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Output results
    std::cout << "Simulation completed.\nDeformations:" << std::endl;
    for (double deformation : simulation.getDeformations()) {
        std::cout << deformation << " ";
    }
    std::cout << "\nInternal Forces:" << std::endl;
    for (double force : simulation.getInternalForces()) {
        std::cout << force << " ";
    }
    std::cout << std::endl;

    return 0;
}
