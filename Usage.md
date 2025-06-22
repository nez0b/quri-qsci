# TE-QSCI Usage Guide

This guide provides comprehensive documentation for using the Time-Evolved Quantum Selected Configuration Interaction (TE-QSCI) implementation with the QURI ecosystem.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Basic Usage](#basic-usage)
4. [Algorithm Variants](#algorithm-variants)
5. [Quantum Chemistry Integration](#quantum-chemistry-integration)
6. [Advanced Features](#advanced-features)
7. [Examples](#examples)
8. [Performance Tuning](#performance-tuning)
9. [Troubleshooting](#troubleshooting)
10. [API Reference](#api-reference)

## Quick Start

Here's a minimal example to get you started with TE-QSCI:

```python
from src import create_qsci_algorithm, QSCIVariant
from quri_parts.core.operator import Operator, pauli_label
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Create a simple H2 Hamiltonian
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})
hamiltonian += Operator({pauli_label("Z1"): -1.25}) 
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})
hamiltonian.constant = -0.5

# Create Hartree-Fock initial state
circuit = QuantumCircuit(2)
circuit.add_X_gate(0)  # Fill first orbital
circuit.add_X_gate(1)  # Fill second orbital
initial_state = GeneralCircuitQuantumState(2, circuit)

# Create sampler
sampler = create_qulacs_vector_concurrent_sampler()

# Run single-time TE-QSCI
algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_states_pick_out=50
)

result = algorithm.run(initial_state, total_shots=1000)
print(f"Ground state energy: {result.ground_state_energy:.6f}")
```

## Installation

### Prerequisites

- Python ≥ 3.12
- QURI Parts with quantum chemistry support
- PySCF for molecular calculations

### Setup

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd te-qsci
   ```

2. **Create virtual environment:**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   uv add quri-parts[qulacs,pyscf,openfermion]
   uv add pyscf
   uv add matplotlib scipy numpy
   ```

### Verify Installation

```python
from src import create_qsci_algorithm, QSCIVariant
print("TE-QSCI successfully installed!")
```

## Basic Usage

### 1. Creating Hamiltonians

#### Manual Construction
```python
from quri_parts.core.operator import Operator, pauli_label

# Build Hamiltonian term by term
hamiltonian = Operator()

# Single-qubit terms
for i in range(4):
    hamiltonian += Operator({pauli_label(f"Z{i}"): -1.0})

# Two-qubit interactions
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})
hamiltonian += Operator({pauli_label("X0 X1"): 0.25})

# Constant energy offset
hamiltonian.constant = -2.0
```

#### From Quantum Chemistry
```python
from pyscf import gto, scf
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner

# Define molecule
mol = gto.M(
    atom="H 0 0 0; H 0 0 0.74",  # H2 molecule
    basis="sto-3g",
    charge=0,
    spin=0
)

# Hartree-Fock calculation
mf = scf.RHF(mol)
mf.run()

# Generate molecular orbitals
active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mol, mf.mo_coeff)

# Create qubit Hamiltonian
hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
    active_space, mo_eint_set, fermion_qubit_mapping=jordan_wigner
)
```

### 2. Creating Initial States

#### Hartree-Fock State
```python
def create_hf_state(n_qubits: int, n_electrons: int) -> GeneralCircuitQuantumState:
    """Create Hartree-Fock state with n_electrons in lowest orbitals."""
    circuit = QuantumCircuit(n_qubits)
    for i in range(n_electrons):
        circuit.add_X_gate(i)
    return GeneralCircuitQuantumState(n_qubits, circuit)

# For H2: 2 electrons, 4 spin orbitals
hf_state = create_hf_state(4, 2)
```

#### Custom Initial States
```python
# Superposition state
circuit = QuantumCircuit(4)
circuit.add_H_gate(0)  # Superposition on first qubit
circuit.add_X_gate(1)  # Occupy second orbital
custom_state = GeneralCircuitQuantumState(4, circuit)

# Entangled state
circuit = QuantumCircuit(4)
circuit.add_H_gate(0)
circuit.add_CNOT_gate(0, 1)
entangled_state = GeneralCircuitQuantumState(4, circuit)
```

### 3. Setting Up Samplers

```python
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
from quri_parts.core.sampling import create_sampler_from_concurrent_sampler

# Concurrent sampler for multiple circuits
concurrent_sampler = create_qulacs_vector_concurrent_sampler()

# Regular sampler for single circuits  
sampler = create_sampler_from_concurrent_sampler(concurrent_sampler)
```

## Algorithm Variants

### 1. Vanilla QSCI

Standard QSCI without time evolution:

```python
algorithm = create_qsci_algorithm(
    QSCIVariant.VANILLA,
    hamiltonian,
    sampler=concurrent_sampler,
    num_states_pick_out=100
)

# Run with multiple input states
result = algorithm.run(
    input_states=[hf_state, custom_state],
    total_shots=5000
)
```

### 2. Single-Time TE-QSCI

Time evolution at a single time point:

```python
algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    sampler=concurrent_sampler,
    evolution_time=1.5,          # Evolution time
    trotter_steps=10,            # Trotter decomposition steps
    num_states_pick_out=85
)

result = algorithm.run(initial_state, total_shots=10000)
```

### 3. Time-Average TE-QSCI

Average over multiple time points:

```python
algorithm = create_qsci_algorithm(
    QSCIVariant.TIME_AVERAGE_TE,
    hamiltonian,
    sampler=concurrent_sampler,
    evolution_times=[0.5, 1.0, 1.5, 2.0],  # Multiple time points
    trotter_steps=8,
    num_states_pick_out=100
)

result = algorithm.run(
    initial_state, 
    shots_per_time=2500  # Shots per time point
)
```

### 4. State Vector TE-QSCI

Exact state vector calculation (requires exact simulator):

```python
algorithm = create_qsci_algorithm(
    QSCIVariant.STATE_VECTOR,
    hamiltonian,
    evolution_time=1.0,
    num_states_pick_out=50
)

# Note: Currently raises NotImplementedError
# Requires exact simulation backend
```

## Quantum Chemistry Integration

### Complete Molecular Study

Here's a complete example for studying a molecular system:

```python
from src import create_qsci_algorithm, QSCIVariant
from pyscf import gto, scf
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner

class MolecularStudy:
    def __init__(self, molecule_geometry: str, basis: str = "sto-3g"):
        self.molecule_geometry = molecule_geometry
        self.basis = basis
        self.setup_molecule()
    
    def setup_molecule(self):
        """Setup molecular calculation."""
        # Parse geometry
        atoms = []
        for line in self.molecule_geometry.split(";"):
            parts = line.strip().split()
            if len(parts) == 4:
                atom, x, y, z = parts
                atoms.append([atom, (float(x), float(y), float(z))])
        
        # Create molecule
        self.mol = gto.M(atom=atoms, basis=self.basis, charge=0, spin=0)
        
        # Hartree-Fock calculation
        self.mf = scf.RHF(self.mol)
        self.mf.run()
        
        # Generate Hamiltonian
        active_space, mo_eint_set = get_spin_mo_integrals_from_mole(
            self.mol, self.mf.mo_coeff
        )
        
        self.hamiltonian, self.jw_mapping = get_qubit_mapped_hamiltonian(
            active_space, mo_eint_set, fermion_qubit_mapping=jordan_wigner
        )
        
        self.active_space = active_space
        self.hf_energy = self.mf.e_tot
        
        print(f"Molecule: {self.molecule_geometry}")
        print(f"HF Energy: {self.hf_energy:.6f} Ha")
        print(f"Hamiltonian terms: {self.hamiltonian.n_terms}")
        print(f"Qubits needed: {2 * active_space.n_active_orb}")
    
    def create_hf_state(self):
        """Create Hartree-Fock initial state."""
        n_qubits = 2 * self.active_space.n_active_orb
        n_electrons = self.active_space.n_active_ele
        
        circuit = QuantumCircuit(n_qubits)
        for i in range(n_electrons):
            circuit.add_X_gate(i)
        
        return GeneralCircuitQuantumState(n_qubits, circuit)
    
    def run_te_qsci_study(self):
        """Run comprehensive TE-QSCI study."""
        initial_state = self.create_hf_state()
        sampler = create_qulacs_vector_concurrent_sampler()
        
        # Single-time study
        results = {}
        evolution_times = [0.5, 1.0, 1.5, 2.0]
        
        for t in evolution_times:
            algorithm = create_qsci_algorithm(
                QSCIVariant.SINGLE_TIME_TE,
                self.hamiltonian,
                sampler=sampler,
                evolution_time=t,
                num_states_pick_out=85
            )
            
            result = algorithm.run(initial_state, total_shots=5000)
            results[t] = result.ground_state_energy
            
            print(f"t={t:.1f}: E={result.ground_state_energy:.6f} Ha")
        
        return results

# Example usage
study = MolecularStudy("H 0 0 0; H 0 0 0.74")  # H2 molecule
energies = study.run_te_qsci_study()
```

### Supported Molecules

The implementation works with any molecule that PySCF can handle:

```python
# Common molecules
molecules = {
    "H2": "H 0 0 0; H 0 0 0.74",
    "H4": "H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 
    "H6": "H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3; H 0 0 4; H 0 0 5",
    "LiH": "Li 0 0 0; H 0 0 1.6",
    "BeH2": "Be 0 0 0; H 0 0 -1.3; H 0 0 1.3",
    "H2O": "O 0 0 0; H 0.757 0.587 0; H -0.757 0.587 0"
}
```

## Advanced Features

### 1. Circuit Analysis with QURI VM

```python
from src import create_vm_enabled_algorithm, LoweringLevel

# Create VM-enabled algorithm
base_algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    evolution_time=1.0
)

vm_algorithm = create_vm_enabled_algorithm(base_algorithm, "STAR")

# Analyze at different levels
logical_analysis = vm_algorithm.analyze(
    [initial_state], 
    total_shots=1000,
    lowering_level=LoweringLevel.LogicalCircuit
)

arch_analysis = vm_algorithm.analyze(
    [initial_state],
    total_shots=1000, 
    lowering_level=LoweringLevel.ArchLogicalCircuit
)

print(f"Logical qubits: {logical_analysis.max_physical_qubit_count}")
print(f"Physical qubits: {arch_analysis.max_physical_qubit_count}")
print(f"Execution time: {arch_analysis.total_latency.value} μs")
```

### 2. Custom Configuration

```python
from src.qsci_algorithms import TimeEvolvedQSCI

# Direct algorithm instantiation with custom parameters
te_qsci = TimeEvolvedQSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    num_states_pick_out=150,
    time_evolution_method="trotter"
)

# Custom Trotter parameters
result = te_qsci.run_single_time(
    initial_state=initial_state,
    evolution_time=2.5,
    total_shots=8000,
    trotter_steps=20  # Higher precision
)
```

### 3. Parameter Sweeps

```python
def parameter_sweep_study(hamiltonian, initial_state, sampler):
    """Study parameter dependencies."""
    results = {}
    
    # Evolution time sweep
    for t in np.linspace(0.1, 3.0, 30):
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=sampler,
            evolution_time=t,
            num_states_pick_out=85
        )
        
        result = algorithm.run(initial_state, total_shots=2000)
        results[f"t_{t:.2f}"] = result.ground_state_energy
    
    # Subspace dimension sweep
    for R in [25, 50, 85, 100, 150, 200]:
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=sampler,
            evolution_time=1.0,
            num_states_pick_out=R
        )
        
        result = algorithm.run(initial_state, total_shots=5000)
        results[f"R_{R}"] = result.ground_state_energy
    
    return results
```

## Examples

### Example 1: H2 Molecule Basic Study

```python
"""Basic H2 molecule study with TE-QSCI."""
import numpy as np
import matplotlib.pyplot as plt
from src import create_qsci_algorithm, QSCIVariant
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Load example (assuming proper quantum chemistry setup)
from examples.h6_molecule_example import H6MoleculeStudy

class H2Study(H6MoleculeStudy):
    def __init__(self):
        self.molecule_geometry = "H 0 0 0; H 0 0 0.74"
        self.jw_mapping = None
        self.active_space = None
        self.exact_ground_state_energy = None
        
        print("Creating H2 Hamiltonian...")
        self.hamiltonian = self._create_h6_hamiltonian()  # Reuses method
        self.hartree_fock_state = self._create_hartree_fock_state()
        self.sampler = create_qulacs_vector_concurrent_sampler()

# Run study
h2_study = H2Study()
times, energies = h2_study.run_single_time_te_qsci_vs_time()

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(times, energies, 'bo-', label='TE-QSCI')
plt.axhline(y=h2_study.exact_ground_state_energy, color='r', 
           linestyle='--', label='HF Reference')
plt.xlabel('Evolution Time')
plt.ylabel('Energy (Ha)')
plt.title('H2 TE-QSCI Energy vs Evolution Time')
plt.legend()
plt.grid(True)
plt.show()
```

### Example 2: Convergence Analysis

```python
"""Analyze TE-QSCI convergence properties."""

def convergence_analysis(hamiltonian, initial_state, sampler):
    """Analyze convergence with respect to shots and subspace size."""
    
    # Shot count convergence
    shot_counts = [500, 1000, 2000, 5000, 10000]
    shot_energies = []
    
    algorithm = create_qsci_algorithm(
        QSCIVariant.SINGLE_TIME_TE,
        hamiltonian,
        sampler=sampler,
        evolution_time=1.0,
        num_states_pick_out=85
    )
    
    for shots in shot_counts:
        result = algorithm.run(initial_state, total_shots=shots)
        shot_energies.append(result.ground_state_energy)
        print(f"Shots: {shots:5d}, Energy: {result.ground_state_energy:.6f}")
    
    # Subspace size convergence
    subspace_sizes = [25, 50, 75, 100, 125, 150]
    subspace_energies = []
    
    for R in subspace_sizes:
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=sampler,
            evolution_time=1.0,
            num_states_pick_out=R
        )
        
        result = algorithm.run(initial_state, total_shots=5000)
        subspace_energies.append(result.ground_state_energy)
        print(f"Subspace: {R:3d}, Energy: {result.ground_state_energy:.6f}")
    
    return (shot_counts, shot_energies), (subspace_sizes, subspace_energies)

# Example usage
(shots_data, subspace_data) = convergence_analysis(
    hamiltonian, initial_state, sampler
)
```

### Example 3: Method Comparison

```python
"""Compare different QSCI variants."""

def method_comparison(hamiltonian, initial_state, sampler):
    """Compare Vanilla QSCI vs TE-QSCI variants."""
    results = {}
    
    # Vanilla QSCI
    vanilla = create_qsci_algorithm(
        QSCIVariant.VANILLA,
        hamiltonian,
        sampler=sampler,
        num_states_pick_out=85
    )
    result = vanilla.run([initial_state], total_shots=5000)
    results["Vanilla QSCI"] = result.ground_state_energy
    
    # Single-time TE-QSCI
    single_time = create_qsci_algorithm(
        QSCIVariant.SINGLE_TIME_TE,
        hamiltonian,
        sampler=sampler,
        evolution_time=1.0,
        num_states_pick_out=85
    )
    result = single_time.run(initial_state, total_shots=5000)
    results["Single-time TE-QSCI"] = result.ground_state_energy
    
    # Time-average TE-QSCI
    time_average = create_qsci_algorithm(
        QSCIVariant.TIME_AVERAGE_TE,
        hamiltonian,
        sampler=sampler,
        evolution_times=[0.5, 1.0, 1.5, 2.0],
        num_states_pick_out=85
    )
    result = time_average.run(initial_state, shots_per_time=1250)
    results["Time-average TE-QSCI"] = result.ground_state_energy
    
    # Print comparison
    print("Method Comparison:")
    print("-" * 40)
    for method, energy in results.items():
        print(f"{method:<25}: {energy:.6f} Ha")
    
    return results
```

## Performance Tuning

### 1. Optimal Parameters

```python
# General guidelines for parameter selection

def get_optimal_parameters(molecule_size, available_shots):
    """Get recommended parameters based on system size."""
    
    if molecule_size <= 4:  # Small molecules (H2, H4)
        return {
            "num_states_pick_out": min(50, available_shots // 50),
            "trotter_steps": 5,
            "evolution_time": 1.0
        }
    elif molecule_size <= 8:  # Medium molecules (H6, LiH)
        return {
            "num_states_pick_out": min(85, available_shots // 100),
            "trotter_steps": 8,
            "evolution_time": 1.0
        }
    else:  # Large molecules
        return {
            "num_states_pick_out": min(150, available_shots // 150),
            "trotter_steps": 10,
            "evolution_time": 1.2
        }

# Example usage
n_qubits = 6  # H6 molecule
total_shots = 10000
params = get_optimal_parameters(n_qubits, total_shots)

algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    sampler=sampler,
    **params
)
```

### 2. Memory Optimization

```python
# For large systems, use streaming approach
def memory_efficient_te_qsci(hamiltonian, initial_state, sampler, 
                            total_shots, batch_size=1000):
    """Memory-efficient TE-QSCI for large systems."""
    
    all_counts = {}
    remaining_shots = total_shots
    
    algorithm = create_qsci_algorithm(
        QSCIVariant.SINGLE_TIME_TE,
        hamiltonian,
        sampler=sampler,
        evolution_time=1.0,
        num_states_pick_out=100
    )
    
    while remaining_shots > 0:
        current_shots = min(batch_size, remaining_shots)
        
        # Run smaller batch
        result = algorithm.run(initial_state, total_shots=current_shots)
        
        # This is a simplified example - actual implementation would
        # need access to measurement counts for accumulation
        print(f"Processed {current_shots} shots, "
              f"Energy: {result.ground_state_energy:.6f}")
        
        remaining_shots -= current_shots
    
    return result
```

### 3. Parallel Execution

```python
import concurrent.futures
from typing import List

def parallel_time_sweep(hamiltonian, initial_state, sampler, 
                       evolution_times: List[float], 
                       shots_per_time: int = 2000):
    """Run TE-QSCI for multiple times in parallel."""
    
    def run_single_time(t):
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=sampler,
            evolution_time=t,
            num_states_pick_out=85
        )
        result = algorithm.run(initial_state, total_shots=shots_per_time)
        return t, result.ground_state_energy
    
    # Parallel execution
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(run_single_time, t) for t in evolution_times]
        results = [future.result() for future in futures]
    
    # Sort by time and return
    results.sort(key=lambda x: x[0])
    times, energies = zip(*results)
    
    return np.array(times), np.array(energies)
```

## Troubleshooting

### Common Issues

#### 1. Import Errors
```bash
# Error: ModuleNotFoundError: No module named 'quri_parts.pyscf'
# Solution: Install quantum chemistry dependencies
uv add quri-parts[pyscf,openfermion]
uv add pyscf
```

#### 2. Complex Number Warnings
```python
# Warning: ComplexWarning: Casting complex values to real discards imaginary part
# This is expected for Hamiltonian coefficients and is handled automatically
# The implementation extracts the real part of angles for rotation gates
```

#### 3. Memory Issues with Large Hamiltonians
```python
# For molecules with >1000 Hamiltonian terms:
# 1. Reduce num_states_pick_out
# 2. Use smaller batch sizes
# 3. Consider state vector methods for exact small systems

# Example for large systems:
algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_states_pick_out=50  # Reduced from default 100
)
```

#### 4. Convergence Issues
```python
# If energies don't converge:
# 1. Increase total_shots
# 2. Increase num_states_pick_out  
# 3. Try different evolution_time values
# 4. Increase trotter_steps for better time evolution

# Convergence check example:
def check_convergence(energies, threshold=1e-4):
    if len(energies) < 3:
        return False
    recent_energies = energies[-3:]
    return np.std(recent_energies) < threshold
```

### Performance Issues

#### 1. Slow Execution
```python
# Optimization strategies:
# 1. Use concurrent samplers
sampler = create_qulacs_vector_concurrent_sampler()

# 2. Reduce precision for preliminary studies
preliminary_params = {
    "num_states_pick_out": 25,
    "total_shots": 1000,
    "trotter_steps": 3
}

# 3. Profile bottlenecks
import time

start_time = time.time()
result = algorithm.run(initial_state, total_shots=5000)
execution_time = time.time() - start_time
print(f"Execution time: {execution_time:.2f} seconds")
```

#### 2. Memory Usage
```python
# Monitor memory usage
import psutil
import os

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # MB

print(f"Memory usage: {get_memory_usage():.1f} MB")
```

### Debugging Tips

```python
# Enable debug output
import logging
logging.basicConfig(level=logging.DEBUG)

# Check Hamiltonian properties
print(f"Hamiltonian terms: {hamiltonian.n_terms}")
print(f"Hermitian: {is_hermitian(hamiltonian)}")
print(f"Largest coefficient: {max(abs(c) for c in hamiltonian.values())}")

# Verify initial state
print(f"Initial state qubits: {initial_state.qubit_count}")
print(f"Circuit depth: {initial_state.circuit.depth}")
print(f"Gate count: {len(initial_state.circuit.gates)}")

# Monitor result properties
print(f"Selected states: {len(result.qsci_result.selected_states)}")
print(f"Subspace dimension: {result.qsci_result.subspace_dimension}")
print(f"Eigenvalue spread: {np.ptp(result.qsci_result.eigenvalues):.6f}")
```

## API Reference

### Core Classes

#### `QSCIVariant`
```python
class QSCIVariant(Enum):
    VANILLA = "vanilla"
    SINGLE_TIME_TE = "single_time_te"
    TIME_AVERAGE_TE = "time_average_te"
    STATE_VECTOR = "state_vector"
```

#### `QSCIResult`
```python
@dataclass
class QSCIResult:
    eigenvalues: Sequence[float]
    eigenstates: Sequence[ComputationalBasisSuperposition]
    selected_states: Sequence[ComputationalBasisState]
    subspace_dimension: int
    total_shots: int
    algorithm_variant: QSCIVariant
    execution_time: Optional[float] = None
    
    @property
    def ground_state_energy(self) -> float:
        """Ground state energy (lowest eigenvalue)."""
        
    @property  
    def ground_state(self) -> ComputationalBasisSuperposition:
        """Ground state (eigenstate for lowest eigenvalue)."""
```

### Factory Functions

#### `create_qsci_algorithm`
```python
def create_qsci_algorithm(
    variant: QSCIVariant,
    hamiltonian: Operator,
    sampler: Optional[ConcurrentSampler] = None,
    num_states_pick_out: Optional[int] = None,
    **variant_kwargs
) -> QSCIAlgorithmBase:
    """
    Factory function to create QSCI algorithm variants.
    
    Args:
        variant: QSCI algorithm variant
        hamiltonian: Target Hamiltonian
        sampler: Quantum sampler (required for sampling variants)
        num_states_pick_out: Number of states to select for subspace
        **variant_kwargs: Variant-specific parameters:
            - evolution_time: Time for evolution (TE variants)
            - evolution_times: List of times (time-average variant)
            - trotter_steps: Trotter decomposition steps
    
    Returns:
        QSCI algorithm instance
    """
```

### Algorithm Classes

#### `VanillaQSCIAlgorithm`
```python
class VanillaQSCIAlgorithm(QSCIAlgorithmBase):
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run vanilla QSCI algorithm."""
```

#### `SingleTimeTeQSCIAlgorithm`
```python
class SingleTimeTeQSCIAlgorithm(QSCIAlgorithmBase):
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_time: float,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        trotter_steps: Optional[int] = None
    ):
        
    def run(
        self,
        initial_state: CircuitQuantumState,
        total_shots: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run single-time TE-QSCI algorithm."""
```

#### `TimeAverageTeQSCIAlgorithm`
```python
class TimeAverageTeQSCIAlgorithm(QSCIAlgorithmBase):
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_times: Sequence[float],
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        trotter_steps: Optional[int] = None
    ):
        
    def run(
        self,
        initial_state: CircuitQuantumState,
        shots_per_time: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run time-average TE-QSCI algorithm."""
```

### Analysis Functions

#### `create_vm_enabled_algorithm`
```python
def create_vm_enabled_algorithm(
    algorithm: QSCIAlgorithmBase,
    architecture: str = "STAR"
) -> QSCIAlgorithmBase:
    """
    Enable VM analysis for QSCI algorithm.
    
    Args:
        algorithm: Base QSCI algorithm
        architecture: Quantum architecture ("STAR", "LINEAR", etc.)
    
    Returns:
        VM-enabled algorithm with analysis capabilities
    """
```

### Parameter Guidelines

| System Size | num_states_pick_out | trotter_steps | evolution_time |
|-------------|-------------------|---------------|----------------|
| 2-4 qubits  | 25-50            | 3-5           | 0.5-1.0        |
| 4-8 qubits  | 50-85            | 5-8           | 1.0-1.5        |
| 8-12 qubits | 85-150           | 8-12          | 1.0-2.0        |
| >12 qubits  | 100-200          | 10-15         | 1.5-2.5        |

### Default Values

- `num_states_pick_out`: 100
- `trotter_steps`: max(1, int(evolution_time * 10))
- `evolution_time`: 1.0
- `total_shots`: 5000 (recommended minimum)

---

This completes the comprehensive usage guide for the TE-QSCI implementation. The framework provides a robust foundation for quantum chemistry research with the QURI ecosystem, offering both ease of use for beginners and advanced features for research applications.