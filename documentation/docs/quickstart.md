# Quick Start

Get up and running with QSCI in minutes! This guide shows you how to run your first QSCI calculation.

## Basic Example

Here's a complete example that runs VanillaQSCI on a simple 2-qubit Hamiltonian:

```python
import numpy as np
from src.qsci_algorithms import VanillaQSCI
from quri_parts.core.operator import Operator, pauli_label
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Create a simple 2-qubit Hamiltonian
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})
hamiltonian += Operator({pauli_label("Z1"): -1.25}) 
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})

# Create uniform superposition initial state |++⟩
circuit = QuantumCircuit(2)
circuit.add_H_gate(0)  # Apply Hadamard to qubit 0
circuit.add_H_gate(1)  # Apply Hadamard to qubit 1
initial_state = GeneralCircuitQuantumState(2, circuit)

# Set up QSCI algorithm
sampler = create_qulacs_vector_concurrent_sampler()
qsci = VanillaQSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    num_states_pick_out=4  # Complete subspace for 2 qubits
)

# Run the algorithm
result = qsci.run([initial_state], total_shots=2000)

# Print results
print(f"Ground state energy: {result.ground_state_energy:.6f}")
print(f"All eigenvalues: {result.eigenvalues}")
print(f"Selected {len(result.selected_states)} configurations")
```

Expected output:
```
Ground state energy: -2.750000
All eigenvalues: [-2.75, -0.75, -0.75, 1.25]
Selected 4 configurations
```

## TE-QSCI Example

Now let's try a time-evolved QSCI calculation:

```python
from src.qsci_algorithms import SingleTimeTE_QSCI

# Same Hamiltonian and initial state as above
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})
hamiltonian += Operator({pauli_label("Z1"): -1.25})
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})

# Create initial state
circuit = QuantumCircuit(2)
circuit.add_H_gate(0)
circuit.add_H_gate(1)
initial_state = GeneralCircuitQuantumState(2, circuit)

# Set up TE-QSCI algorithm
sampler = create_qulacs_vector_concurrent_sampler()
te_qsci = SingleTimeTE_QSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    evolution_time=1.0,  # Evolution time parameter
    num_states_pick_out=4
)

# Run the algorithm
result = te_qsci.run([initial_state], total_shots=2000)

print(f"TE-QSCI Ground state energy: {result.ground_state_energy:.6f}")
```

## Step-by-Step Explanation

### 1. Create the Hamiltonian

```python
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})  # σ_z^0 term
hamiltonian += Operator({pauli_label("Z1"): -1.25})  # σ_z^1 term  
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5}) # σ_z^0 ⊗ σ_z^1 term
```

This creates the Hamiltonian:
```
H = -1.25 σ_z^0 - 1.25 σ_z^1 + 0.5 σ_z^0 σ_z^1
```

### 2. Create the Initial State

```python
circuit = QuantumCircuit(2)
circuit.add_H_gate(0)  # |+⟩ = (|0⟩ + |1⟩)/√2
circuit.add_H_gate(1)  # |+⟩ = (|0⟩ + |1⟩)/√2
initial_state = GeneralCircuitQuantumState(2, circuit)
```

This creates the uniform superposition state:

$$|\psi\rangle = |++\rangle = \frac{|00\rangle + |01\rangle + |10\rangle + |11\rangle}{2}$$

### 3. Configure the Algorithm

```python
qsci = VanillaQSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    num_states_pick_out=4  # Use complete subspace
)
```

Key parameters:
- `hamiltonian`: The quantum Hamiltonian to diagonalize
- `sampler`: Quantum measurement simulator
- `num_states_pick_out`: Number of configurations to select (4 = 2² for complete subspace)

### 4. Run and Analyze

```python
result = qsci.run([initial_state], total_shots=2000)
```

The result contains:
- `ground_state_energy`: Lowest eigenvalue
- `eigenvalues`: All computed eigenvalues
- `selected_states`: Configurations used in the calculation

## More Examples

### Time-Average TE-QSCI

```python
from src.qsci_algorithms import TimeAverageTE_QSCI

te_qsci_avg = TimeAverageTE_QSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    evolution_times=[0.5, 1.0, 1.5],  # Multiple time points
    num_states_pick_out=4
)

result = te_qsci_avg.run([initial_state], total_shots=3000)
print(f"Time-average energy: {result.ground_state_energy:.6f}")
```

### State Vector TE-QSCI

```python
from src.qsci_algorithms import StateVectorTE_QSCI

te_qsci_exact = StateVectorTE_QSCI(
    hamiltonian=hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_states_pick_out=4
)

result = te_qsci_exact.run([initial_state], total_shots=0)  # No shots needed
print(f"Exact energy: {result.ground_state_energy:.10f}")
```

## Best Practices

### 1. Use Uniform Superposition

For optimal performance, always use uniform superposition initial states:

```python
def create_uniform_superposition(n_qubits):
    circuit = QuantumCircuit(n_qubits)
    for i in range(n_qubits):
        circuit.add_H_gate(i)
    return GeneralCircuitQuantumState(n_qubits, circuit)
```

### 2. Complete Subspace Coverage

For small systems, use complete subspace coverage:

```python
num_states_pick_out = 2**n_qubits  # Complete Hilbert space
```

### 3. Sufficient Shot Counts

Use adequate shot counts for stable results:

```python
total_shots = 2000  # Minimum for 2-qubit systems
total_shots = 3000  # Better for TE-QSCI variants
```

### 4. Validate Results

Always validate against exact solutions for small systems:

```python
from scipy.sparse.linalg import eigsh

# Convert to matrix and solve exactly
H_matrix = hamiltonian_to_matrix(hamiltonian)  # Helper function
exact_energy, _ = eigsh(H_matrix, k=1, which='SA')
print(f"QSCI energy: {result.ground_state_energy:.6f}")
print(f"Exact energy: {exact_energy[0]:.6f}")
print(f"Error: {abs(result.ground_state_energy - exact_energy[0]):.2e}")
```

## Common Issues

### 1. Import Errors
Make sure QURI Parts is installed:
```bash
pip install quri-parts qulacs
```

### 2. Numerical Precision
Use appropriate tolerances for comparisons:
```python
assert abs(energy_error) < 1e-8  # Machine precision for uniform superposition
```

### 3. Memory Issues
For large systems, reduce the number of selected states:
```python
num_states_pick_out = min(100, 2**n_qubits)  # Limit subspace size
```

## Next Steps

Now that you've run your first QSCI calculation:

1. **[Algorithm Overview](algorithms/overview.md)** - Learn about different QSCI variants
2. **[Examples](examples/basic_usage.md)** - Explore more complex examples
3. **[Testing](testing/strategy.md)** - Understand the validation framework
4. **[API Reference](api/)** - Dive into the technical details

## Interactive Examples

Try these examples in a Jupyter notebook:

- **[H2 Molecule](examples/h2_molecule.md)** - Quantum chemistry application
- **[TFIM Models](examples/tfim_models.md)** - Condensed matter physics
- **[Performance Analysis](examples/performance.md)** - Scaling and optimization