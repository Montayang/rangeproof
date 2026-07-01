# Rarus: Zero-Knowledge Polynomial Range Proofs

This repository contains the official implementation of **Rarus**, a zero-knowledge polynomial range proof protocol, along with baseline comparisons. 

## 1. Overview & Features

This repository implements the final protocol for polynomial range proofs based on the `mcl` library (leveraging its highly optimized MSM and pairing functions). It also includes implementations of baseline protocols for performance comparison.

**Core Cryptographic Tools Implemented:**

* **NTT / INTT:** Univariate and Bivariate (non-recursive) Number Theoretic Transforms.
* **Polynomial Operations:** Univariate and Bivariate polynomial product and division.
* **KZG Commitments:** Univariate, Bivariate, and their respective batching versions.
* **Sumset Representation.**

**Included Baselines:**

* Plookup
* Bulletproofs

---

## 2. Prerequisites

The system is developed and tested on Linux (e.g., Ubuntu 22.04 LTS). It requires standard C++ build tools and cryptographic libraries.

Please install the following system dependencies before building:

```bash
sudo apt-get update
sudo apt-get install -y libgmp-dev libssl-dev cmake make g++ git
```

---

## 3. Environment Setup & Build Instructions

The protocol strictly depends on the `mcl` library (V3+). Please follow these exact steps to clone the dependencies, compile the `mcl` static library, and build the main Rarus protocol.

### Step 3.1: Clone the Repository

```bash
git clone https://github.com/Montayang/rangeproof.git
cd rangeproof
```

### Step 3.2: Fetch and Build `mcl`

The project expects the `mcl` library to be located in the root directory. 

```bash
# Clone the mcl library directly into the project root
git clone https://github.com/herumi/mcl.git

# Compile the mcl core library
cd mcl
make clean
make -j4
cd ..
```

*Verification:* Ensure that `mcl/lib/libmcl.a` is successfully generated before proceeding.

### Step 3.3: Build the Main Protocol

We use CMake to configure and build the project.

```bash
mkdir -p build
cd build
rm -rf *
cmake ..
make -j4
```

---

## 4. Kick-the-Tires (Quick Start)

To verify that the system is built correctly and your environment is fully functional, run the main protocol executable from inside the `build` directory:

```bash
./main
```

**Expected Output:**
The program will execute the zero-knowledge PIOP setup and generation phases. Within a few seconds, you should see logs indicating the completion of each step, concluding with the final benchmark metrics:

```text
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
```

---

## 5. Detailed Benchmarks (Results Reproduced)

To maintain absolute transparency and avoid obfuscating the core logic with complex wrapper scripts, we expose the benchmarking parameters directly in the source code. Reviewers can verify and modify the parameters (polynomial size, range, # of elements) directly.

### Parameter Mapping to the Paper

* **Table 3 parameters** map to the results in **Table 8**.
* **Table 4 parameters** map to the results in **Table 9**.
* The outputs correspond to the data points in **Figure 3 and Figure 4**.

### How to Modify Parameters

Before running `make`, you can adjust the testing parameters by modifying the following lines in the source code:

* **Rarus (Main Protocol):** Edit line `211` in `zk_PIOP.cpp`.
* **Bulletproofs Baseline:** Edit line `231` in `bulletproof.cpp`.
* **Plookup Baseline:** Edit line `196` in `plookup.cpp`.

After modifying the parameters, simply rebuild the executables:
```bash
cd build
make -j4
```

### Running the Evaluations
```bash
# Evaluate Rarus
./main

# Evaluate Bulletproofs
./run_bulletproof

# Evaluate Plookup
./run_plookup
```

---

## 6. Directory Structure

* `zk_PIOP.cpp`: The core protocol implementation (Main entry point).
* `kzg.cpp` / `kzg.h`: Implementation of KZG polynomial commitments.
* `subgroups.cpp` / `subgroups.h`: Subgroup definitions and operations.
* `bulletproof.cpp`: Baseline implementation of Bulletproofs.
* `plookup.cpp`: Baseline implementation of Plookup.
* `mcl/`: The underlying high-performance cryptographic pairing library.
