#include <gtest/gtest.h>
#include "../include/MeshReader.h"
#include "../include/Simulation.h"
#include <cmath>
#include <vector>

class SimulationTest : public ::testing::Test {
protected:
    MeshReader* meshReader;  // Use a pointer to MeshReader
    std::vector<int> blockedNodes;  // To store boundary node IDs
    std::vector<int> forceNodes;  // To store for node IDs
    std::vector<std::vector<int>> blockedVectors;
    std::vector<std::vector<double>> forceVectors;
    double E;
    std::vector<double> AVector;

    // Constructor
    SimulationTest() : meshReader(nullptr), E(0.0) {}

    // Set up resources before running each test
    void SetUpSimulation(
        const std::string& filePath, double EValue,
        const std::vector<double>& areaVector,
        const std::vector<std::vector<int>>& blockedVecs,
        const std::vector<std::vector<double>>& forceVecs
    ) {
        // Set material properties
        E = EValue;
        AVector = areaVector;

        // Initialize MeshReader dynamically
        meshReader = new MeshReader(filePath);
        meshReader->readFile();

        // Retrieve boundary nodes and assign vectors
        blockedNodes = meshReader->getboundaryNodeIds();
        forceNodes = meshReader->getForceNodeIds();
        blockedVectors = blockedVecs;
        forceVectors = forceVecs;
    }

    // Clean up after each test
    virtual void TearDown() override {
        delete meshReader;
        meshReader = nullptr;
    }
};


TEST_F(SimulationTest, OneBeam2DTest) {
    // Single steel beam, 1 Area for cross section
    double E = 2.1E8;
    std::vector<double> AVector(1, 0.008770);

    // Node ids for force and boundaries are defined in mail file: here define specifically dofs
    // size of blockedNodes in mail file must comply with vector size blockedVectors (similarly forceVectors)
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0},
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, 10, 0},
    };

    // Set up the simulation with a specific file path and parameters
    SetUpSimulation("../testfiles/1trussin2d.mail", E, AVector, blockedVectors, forceVectors);

    // Create and run the simulation
    Simulation simulation(*meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Check assertion about size of deformation and check deformation values
    std::vector<double> expectedDeformation = {0, 0, 0, 0.000108595};
    ASSERT_EQ(simulation.getDeformations().size(), expectedDeformation.size());
    for (size_t i = 0; i < expectedDeformation.size(); ++i) {
        EXPECT_NEAR(simulation.getDeformations()[i], expectedDeformation[i], 1e-6);
    }
}

TEST_F(SimulationTest, OneBeam3DTest) {
    // Single steel beam, 1 Area for cross section
    double E = 2.1E8;
    std::vector<double> AVector(1, 0.008770);

    // Node ids for force and boundaries are defined in mail file: here define specifically dofs
    // size of blockedNodes in mail file must comply with vector size blockedVectors (similarly forceVectors)
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0},
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, 10, 0},
    };

    // Set up the simulation with a specific file path and parameters
    SetUpSimulation("../testfiles/1trussin3d.mail", E, AVector, blockedVectors, forceVectors);

    // Create and run the simulation
    Simulation simulation(*meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Check assertion about size of deformation and check deformation values
    std::vector<double> expectedDeformation = {0, 0, 0, 0, 0.000108595, 0};
    ASSERT_EQ(simulation.getDeformations().size(), expectedDeformation.size());
    for (size_t i = 0; i < expectedDeformation.size(); ++i) {
        EXPECT_NEAR(simulation.getDeformations()[i], expectedDeformation[i], 1e-6);
    }
}


TEST_F(SimulationTest, Truss2din3d) {
    // 4 Beams in steel with same cross section
    double E = 2.1E8;
    std::vector<double> AVector = {0.008770, 0.008770, 0.008770, 0.008770};

    // Node ids for force and boundaries are defined in mail file: here define specifically dofs
    // size of blockedNodes in mail file must comply with vector size blockedVectors (similarly forceVectors)
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0},  // Corresponding to blockedNode 1 (fixed)
        {0, 0, 0},   // Corresponding to blockedNode 2 (fixed)
        {1, 0, 1},  // Corresponding to blockedNode 3 (in plane boundary)
        {1, 0, 1},   // Corresponding to blockedNode 4 (in plane boundary)
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, 0, -1000},
    };

    // Set up the simulation with a specific file pat for the meshh and parameters
    SetUpSimulation("../testfiles/truss2din3d.mail", E, AVector, blockedVectors, forceVectors);

    // Create and run the simulation
    Simulation simulation(*meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Check assertion about size of deformation and check deformation values for N4
    std::vector<double> checkedDeformation = {0.0108595, 0, -0.0415749};
    ASSERT_EQ(simulation.getDeformations().size(), 4 * 3);
    for (size_t i = 0; i < checkedDeformation.size(); ++i) {
        EXPECT_NEAR(simulation.getDeformations()[i + 3 * 3], checkedDeformation[i], 1e-6);
    }
}


TEST_F(SimulationTest, BridgeTruss3d) {
    // Define custom materials and sections
    double E = 1000;
    double Abottom = 2;
    double Atop = 10;
    double Avertical = 3;
    double Adiagonal = 1;
    std::vector<double> AVector = {
        Abottom, Abottom, Abottom, Abottom, Abottom, Abottom,
        Atop, Atop, Atop, Atop, Atop, Atop, 
        Avertical, Avertical, Avertical, Avertical, Avertical,
        Adiagonal, Adiagonal, Adiagonal, Adiagonal
    };

    // Node ids for force and boundaries are defined in mail file: here define specifically dofs
    // size of blockedNodes in mail file must comply with vector size blockedVectors (similarly forceVectors)
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0},
        {1, 0, 0},
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, -10, 0},
        {0, -10, 0},
        {0, -16, 0},
        {0, -10, 0},
        {0, -10, 0},
    };

    // Set up the simulation with a specific file pat for the meshh and parameters
    SetUpSimulation("../testfiles/bridgetruss3dchp21.mail", E, AVector, blockedVectors, forceVectors);

    // Create and run the simulation
    Simulation simulation(*meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Check assertion about size of deformation and check deformation at midspan point N19 max negative displacement
    std::vector<double> checkedDeformation = {-2.4219383828501};
    ASSERT_EQ(simulation.getDeformations().size(), 12 * 3);
    for (size_t i = 0; i < checkedDeformation.size(); ++i) {
        EXPECT_NEAR(simulation.getDeformations()[i + 19], checkedDeformation[i], 1e-6);
    }
}

TEST_F(SimulationTest, RadialTruss3d) {
    // Define steel materials and section for all 9 beams
    double E = 2.1E8;
    std::vector<double> AVector(9, 0.008770);

    // Node ids for force and boundaries are defined in mail file: here define specifically dofs
    // size of blockedNodes in mail file must comply with vector size blockedVectors (similarly forceVectors)
    std::vector<std::vector<int>> blockedVectors = {
        {0, 0, 0},  // Corresponding to blockedNode 1
        {0, 0, 0},   // Corresponding to blockedNode 2
        {0, 0, 0},  // Corresponding to blockedNode 5
        {0, 0, 0},   // Corresponding to blockedNode 6
    };
    std::vector<std::vector<double>> forceVectors = {
        {0, 0, -100},
    };

    // Set up the simulation with a specific file pat for the meshh and parameters
    SetUpSimulation("../testfiles/radialtruss3d.mail", E, AVector, blockedVectors, forceVectors);

    // Create and run the simulation
    Simulation simulation(*meshReader, blockedVectors, forceVectors, E, AVector);
    simulation.run();

    // Check assertion about size of deformation and check deformation at midspan point N19 max negative displacement
    std::vector<double> checkedDeformation = {6.27156e-20, -3.28405e-20, -0.000156456};
    ASSERT_EQ(simulation.getDeformations().size(), 6 * 3);
    for (size_t i = 0; i < checkedDeformation.size(); ++i) {
        EXPECT_NEAR(simulation.getDeformations()[i + 6], checkedDeformation[i], 1e-6);
    }
}

// Run all tests
int main(int argc, char **argv) {
    // Initialise google test suite
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(print_time) = true;
    ::testing::GTEST_FLAG(throw_on_failure) = true;

    // Run all tests
    int result = RUN_ALL_TESTS();

    // Ensure all output is flushed to the console
    std::cout << std::flush;
    return result;
}
