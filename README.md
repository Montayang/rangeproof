# Rarus: Zero-Knowledge Polynomial Range Proofs

This repository contains the official implementation of Rarus, a zero-knowledge polynomial range proof protocol, along with baseline comparisons.

This artifact is prepared for the USENIX Security 2026 Artifact Evaluation. It provides the necessary source code, dependencies, and scripts to reproduce the primary claims made in the paper regarding prover time, verifier time, and proof size.

1. Overview & Features
This repository implements the final protocol for polynomial range proofs based on the mcl library (leveraging its highly optimized MSM and pairing functions). It also includes implementations of baseline protocols for performance comparison.

Core Cryptographic Tools Implemented:

NTT / INTT: Univariate and Bivariate (non-recursive) Number Theoretic Transforms.

Polynomial Operations: Univariate and Bivariate polynomial product and division.

KZG Commitments: Univariate, Bivariate, and their respective batching versions.

Sumset Representation.

Included Baselines:

Plookup

Bulletproofs

2. Prerequisites
The system is developed and tested on Linux (e.g., Ubuntu 22.04 LTS). It requires standard C++ build tools and cryptographic libraries.

Please install the following system dependencies before building:

Bash
sudo apt-get update
sudo apt-get install -y libgmp-dev libssl-dev cmake make g++ git
3. Environment Setup & Build Instructions
The protocol strictly depends on the mcl library (V3+). Please follow these exact steps to clone the dependencies, compile the mcl static library, and build the main Rarus protocol.

Step 3.1: Clone the Repository
Bash
git clone https://github.com/YOUR_USERNAME/rangeproof-main.git
cd rangeproof-main
(Note: Replace YOUR_USERNAME with your actual GitHub username or remove this line if they are already in the directory).

Step 3.2: Fetch and Build mcl
The project expects the mcl library to be located in the root directory.

Bash
# Clone the mcl library directly into the project root
git clone https://github.com/herumi/mcl.git

# Compile the mcl core library
cd mcl
make clean
make -j4
cd ..
Verification: Ensure that mcl/lib/libmcl.a is successfully generated before proceeding.

Step 3.3: Build the Main Protocol
We use CMake to configure and build the project.

Bash
mkdir -p build
cd build
rm -rf * # Ensure a clean build environment
cmake ..
make -j4
4. Kick-the-Tires (Quick Start)
To verify that the system is built correctly and your environment is fully functional, run the main protocol executable from inside the build directory:

Bash
./main
Expected Output:
The program will execute the zero-knowledge PIOP setup and generation phases. Within a few seconds, you should see logs indicating the completion of each step, concluding with the final benchmark metrics:

Plaintext
Setup Completed.
Step1 Completed.
Step2 Completed.
F1 Completed.
F2 Constructing Completed.
F2 Completed.
F3 Constructing Completed.
Step3 Preparing Completed.
Step3 Completed.
Step4 Completed.
1
Prover Time: <time_in_seconds>
Verifier Time: <time_in_seconds>
Proof Size: <size_in_KB> KB
5. Detailed Benchmarks (Results Reproduced)
(Author's Note to AEC: The following instructions map the executable outputs to the claims made in our paper.)

Main Protocol (Rarus)
Running ./main evaluates the core protocol.

Prover Time & Verifier Time: The outputted time (in seconds) corresponds to the performance evaluations shown in Table X / Figure Y of the paper.

Proof Size: The outputted proof size (in KB) corresponds to the communication complexity claims in Section Z.

Baselines (Plookup & Bulletproofs)
(Note: 如果你已经把 plookup 和 bulletproof 写进了 CMakeLists.txt 并且生成了对应的可执行文件，比如 run_plookup 和 run_bulletproof，请在这里写上运行指令。如果没有，可以说还在整理中，或者如果有其他脚本可以运行它们，请补充在这里。例如：)

To run the baseline comparisons evaluated in the paper:

Bash
# Run Bulletproofs baseline
./run_bulletproof

# Run Plookup baseline
./run_plookup
(Please compare these outputs with the baseline data in Table X of our paper.)

6. Directory Structure
zk_PIOP.cpp: The core protocol implementation (Main entry point).

kzg.cpp / kzg.h: Implementation of KZG polynomial commitments.

subgroups.cpp / subgroups.h: Subgroup definitions and operations.

bulletproof.cpp: Baseline implementation of Bulletproofs.

plookup.cpp: Baseline implementation of Plookup.

mcl/: The underlying high-performance cryptographic pairing library.
