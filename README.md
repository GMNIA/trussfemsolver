# trussfemsolver
Custom implementation of truss beam solver in cpp 14 environment.

---
## Libraries

- eigen-3.4.0
- googletest
---

## Project scope
This project starts from a structural engineering + python developer background.
The idea of the project is meant as a numerical framework of standard finite element methods to be used in other projects, to compare or to grab results in a estremely fast way.

For example in order to obtain results for training data in a reinforcement learning project, it might be faster usto use this application, rather than install third-party software or proprietary, use API to run the code, submit and solve and then read the results from file.

Example of a mesh analysed and solved by the trussfemsolver (see testfiles\radialtruss3d.mail):

![image](https://github.com/user-attachments/assets/27e80cf7-15ce-48c7-a3c4-c1408e40017e)


# TESTS: How to Build and Run Tests with CMake & GoogleTest

## 1️⃣ Prerequisites
- Ensure you have **CMake (≥3.10)** installed.
- Ensure you have a **C++ compiler** (e.g., `g++` or `clang++`).
- GoogleTest should be cloned into trussfemsolver/googletest (clone official repo) and included in your project via `add_subdirectory(googletest)`.

---

## 2️⃣ Navigate to the Project Directory
```sh
cd /path/to/trussfemsolver
```

---

## 3️⃣ Create and Enter the `build/` Directory
```sh
mkdir -p build && cd build
```
💡 This ensures a clean build by keeping compilation files separate from source code.

---

## 4️⃣ Configure the Project with CMake
```sh
cmake ..
```
🔹 This detects your compiler and prepares the build system.

---

## 5️⃣ Compile the Project
```sh
cmake --build .
```
or if using `make`:
```sh
make
```

📌 **This step compiles**:
- The `simulation_test` executable
- Any linked dependencies (e.g., GoogleTest)

---

## 6️⃣ Run the Test Executable
```sh
./simulation_test
```
For **Windows (PowerShell or CMD)**:
```sh
.\simulation_test.exe
```
✅ This executes all GoogleTest test cases.

---

## 7️⃣ Run a Specific Test (Optional)
To run only a **specific test case** (e.g., `BridgeTruss3d`):
```sh
./simulation_test --gtest_filter=SimulationTest.BridgeTruss3d
```

---

## 8️⃣ Enable Debug Mode (Optional)
If debugging is needed, enable debugging when configuring:
```sh
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build .
```
Run with **GDB** (Linux/macOS):
```sh
gdb ./simulation_test
run
```
For a backtrace after a crash:
```sh
bt
```

---

## 9️⃣ Clean and Rebuild (If Needed)
If you face build issues, **remove the cache** and rebuild:
```sh
rm -rf CMakeCache.txt CMakeFiles/
cmake ..
cmake --build .
```
To completely reset and recompile:
```sh
rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
```

---

## 🔟 Quick Reference Table

| **Step** | **Command** |
|----------|------------|
| **1. Create & enter build directory** | `mkdir -p build && cd build` |
| **2. Run CMake configuration** | `cmake ..` |
| **3. Build the project** | `cmake --build .` |
| **4. Run the tests** | `./simulation_test` |
| **5. Run a specific test** | `./simulation_test --gtest_filter=SimulationTest.BridgeTruss3d` |
| **6. Enable debug mode** | `cmake .. -DCMAKE_BUILD_TYPE=Debug && cmake --build .` |
| **7. Clean and rebuild** | `rm -rf build && mkdir build && cd build && cmake .. && cmake --build .` |

---

## 📌 Notes
- Ensure that **CMakeLists.txt** includes GoogleTest properly.
- Verify that test files (e.g., `bridgetruss3dchp21.mail`) exist in the correct location.
- Verify that GoogleTest repo is cloned and correctly included and built.



# CLI: build and run an example file returning results


Run this command in the terminal, the cli will run the example file under:

> trussfemsolver/testfiles/example.mail

**g++ -o cli src/cli.cpp src/Simulation.cpp src/MeshReader.cpp -Ilib/eigen-3.4.0 -std=c++14 && ./cli**

The output of the deformation and internal forces vector should be:

> Deformations: <br>
> 0 0 0 0.000108595 <br>
> Internal Forces: <br>
> 0 -10 0 10 <br>
