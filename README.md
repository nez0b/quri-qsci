# TE-QSCI: Time-Evolved Quantum Selected Configuration Interaction

This repository implements the **Time-Evolved Quantum Selected Configuration Interaction (TE-QSCI)** algorithm using the QURI Parts framework, based on the research paper "Quantum-selected configuration interaction with time-evolved state" (arXiv:2412.13839v2).

## Overview

TE-QSCI is a quantum-classical hybrid algorithm that uses time-evolved quantum states as input for Quantum Selected Configuration Interaction (QSCI). Unlike traditional approaches that require optimization-based state preparation, TE-QSCI leverages the natural time evolution of quantum states to systematically generate important electronic configurations.

### Key Features

- **Optimization-free state preparation**: Uses time evolution instead of variational optimization
- **Multiple algorithm variants**: Single-time, time-average, and state vector implementations
- **QURI-algo compatible**: Implements the standard algorithm interface from quri-algo
- **Architecture analysis**: Supports LogicalCircuit and ArchLogicalCircuit analysis with QURI VM
- **STAR architecture demonstration**: Optimized for STAR quantum computing architecture

## Algorithm Variants

### 1. Vanilla QSCI
Standard QSCI algorithm for comparison with TE-QSCI methods.

### 2. Single-time TE-QSCI
Uses a single time-evolved state `|ψ(t)⟩ = e^(-iĤt)|ψ_I⟩` as input for QSCI.

### 3. Time-average TE-QSCI
Combines measurements from multiple time-evolved states at different evolution times.

### 4. State vector TE-QSCI
Direct state vector calculation without quantum sampling (for exact simulation).

## Installation

1. Clone this repository:
```bash
git clone https://github.com/your-repo/te-qsci.git
cd te-qsci
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Install QURI Parts with quantum backends:
```bash
pip install "quri-parts[qulacs]"
```

## Usage

### Basic Example

```python
from qsci_algo_interface import create_qsci_algorithm, QSCIVariant
from quri_parts.core.operator import pauli_label
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Create Hamiltonian
hamiltonian = pauli_label("Z0 Z1") + 0.5 * pauli_label("X0 X1")

# Create sampler
sampler = create_qulacs_vector_concurrent_sampler()

# Create TE-QSCI algorithm
algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian=hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_states_pick_out=50
)

# Create initial state (Hartree-Fock)
from quri_parts.circuit import QuantumCircuit
from quri_parts.core.state import CircuitQuantumState

circuit = QuantumCircuit(2)
circuit.add_X_gate(0)  # Fill first orbital
initial_state = CircuitQuantumState(2, circuit)

# Run algorithm
result = algorithm.run(initial_state, total_shots=1000)
print(f"Ground state energy: {result.ground_state_energy}")
```

### H6 Molecule Example

Run the complete H6 molecule study reproducing paper results:

```python
from h6_molecule_example import H6MoleculeStudy

study = H6MoleculeStudy()
study.run_full_study()
```

This will reproduce:
- **Figure 1**: TE-QSCI energy vs evolution time
- **Figure 2**: Architecture analysis (LogicalCircuit vs ArchLogicalCircuit)
- **Table II**: Subspace dimension study
- **Table III**: Method comparison

### Architecture Analysis

Analyze quantum resource requirements on different architectures:

```python
from qsci_vm_analysis import create_vm_enabled_algorithm
from qsci_algo_interface import LoweringLevel

# Create VM-enabled algorithm
vm_algorithm = create_vm_enabled_algorithm(base_algorithm, "STAR")

# Analyze at different levels
logical_analysis = vm_algorithm.analyze(
    [input_state], 1000, LoweringLevel.LogicalCircuit
)

arch_analysis = vm_algorithm.analyze(
    [input_state], 1000, LoweringLevel.ArchLogicalCircuit
)

print(f"Logical circuit latency: {logical_analysis.total_latency}")
print(f"ArchLogical circuit latency: {arch_analysis.total_latency}")
print(f"SWAP overhead: {arch_analysis.total_swap_overhead}")
```

## File Structure

```
te-qsci/
├── qsci_algorithms.py          # Core QSCI and TE-QSCI implementations
├── qsci_algo_interface.py      # QURI-algo compatible interfaces
├── qsci_vm_analysis.py         # QURI VM integration for circuit analysis
├── h6_molecule_example.py      # H6 molecule study reproducing paper results
├── requirements.txt            # Python dependencies
├── README.md                   # This file
├── 2302.11320v1.pdf           # Original QSCI paper
└── arXiv-2412.13839v2/        # TE-QSCI paper and figures
    ├── 2412.13839v2.pdf
    ├── main.tex
    └── [figure files]
```

## Algorithm Classes

### Core Algorithm Classes

- `QSCIBase`: Base class for all QSCI algorithms
- `VanillaQSCI`: Standard QSCI implementation
- `TimeEvolvedQSCI`: TE-QSCI implementation with multiple variants

### QURI-algo Compatible Classes

- `QSCIAlgorithmBase`: Base class implementing quri-algo interface
- `VanillaQSCIAlgorithm`: Vanilla QSCI with quri-algo interface
- `SingleTimeTeQSCIAlgorithm`: Single-time TE-QSCI
- `TimeAverageTeQSCIAlgorithm`: Time-average TE-QSCI
- `StateVectorTeQSCIAlgorithm`: State vector TE-QSCI

### VM Analysis Classes

- `VMCircuitAnalyzer`: Circuit analysis with architecture awareness
- `QSCIVMAnalysis`: Enhanced QSCI analysis using QURI VM
- `VMEnabledQSCIAlgorithm`: QSCI with VM analysis capabilities

## Supported Architectures

- **LogicalCircuit**: Architecture-independent analysis
- **ArchLogicalCircuit**: Architecture-specific analysis with STAR topology
- **STAR Architecture**: Optimized for star-connected quantum processors

## Key Results from H6 Study

When running the H6 molecule example, you should observe:

1. **Optimal evolution time**: Around t = 1.0 for single-time TE-QSCI
2. **Subspace scaling**: Energy error decreases with larger subspace dimension R
3. **Method comparison**: TE-QSCI outperforms vanilla QSCI with HF input
4. **Architecture overhead**: STAR topology reduces SWAP gate requirements

## Theoretical Background

TE-QSCI leverages the intuition that time evolution by the Hamiltonian creates electron excitations of various orders:

```
|ψ(t)⟩ = e^(-iĤt)|ψ_I⟩ = |ψ_I⟩ - iĤt|ψ_I⟩ + (−iĤt)²/2!|ψ_I⟩ + ...
```

The k-th order term includes excitations up to 2k-th order, naturally generating important configurations for selected CI calculations.

## References

1. **TE-QSCI Paper**: "Quantum-selected configuration interaction with time-evolved state" (arXiv:2412.13839v2)
2. **Original QSCI**: "Subspace diagonalization by quantum-selected configuration interaction" (arXiv:2302.11320v1)
3. **QURI Parts**: https://quri-parts.qunasys.com/
4. **QURI SDK**: https://quri-sdk.qunasys.com/

## License

This implementation is provided for research and educational purposes. Please cite the original papers when using this code.

## Contributing

Contributions are welcome! Please ensure your code follows the existing style and includes appropriate tests.

## Contact

For questions about this implementation, please open an issue in the repository.