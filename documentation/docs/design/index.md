# Design & Implementation

This section covers the technical design, architecture decisions, and implementation details of the QSCI project.

## Architecture Overview

The QSCI project is built with a modular, extensible architecture that integrates seamlessly with the QURI ecosystem while providing clear separation of concerns.

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    QSCI Application Layer                   │
├─────────────────────────────────────────────────────────────┤
│  VanillaQSCI  │  SingleTimeTE_QSCI  │  TimeAverageTE_QSCI  │
├─────────────────────────────────────────────────────────────┤
│                 QSCIBase (Abstract Base)                   │
├─────────────────────────────────────────────────────────────┤
│             QURI Algo Interface Layer                      │
├─────────────────────────────────────────────────────────────┤
│        QURI Parts (Operators, States, Circuits)           │
├─────────────────────────────────────────────────────────────┤
│              QURI VM (Circuit Analysis)                    │
├─────────────────────────────────────────────────────────────┤
│         Backend Layer (Qulacs, IBM, etc.)                 │
└─────────────────────────────────────────────────────────────┘
```

## Design Principles

### 1. **Modularity and Extensibility**
- Clear separation between algorithm variants
- Abstract base classes for common functionality
- Plugin-style architecture for new variants

### 2. **QURI Ecosystem Integration**
- Native QURI Parts operator and state handling
- Compatible with quri-algo algorithm interfaces
- QURI VM integration for circuit analysis

### 3. **Performance Optimization**
- Efficient sparse matrix operations
- Concurrent sampling strategies
- Memory-efficient state vector handling

### 4. **Scientific Rigor**
- Comprehensive validation against exact solutions
- Machine precision accuracy verification
- Robust numerical handling

## Core Components

### Algorithm Base Classes

#### QSCIBase
The foundation class providing common functionality:

```python
class QSCIBase(ABC):
    """Abstract base class for all QSCI variants."""
    
    def __init__(self, hamiltonian, sampler, num_states_pick_out=None):
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        self.num_states_pick_out = num_states_pick_out
    
    @abstractmethod
    def run(self, input_states, total_shots, **kwargs):
        """Execute the QSCI algorithm."""
        pass
```

#### TimeEvolvedQSCI
Base class for time evolution variants:

```python
class TimeEvolvedQSCI(QSCIBase):
    """Base class for time-evolved QSCI algorithms."""
    
    def __init__(self, hamiltonian, sampler, evolution_time, 
                 time_evolution_method="exact", **kwargs):
        super().__init__(hamiltonian, sampler, **kwargs)
        self.evolution_time = evolution_time
        self.time_evolution_method = time_evolution_method
```

### State Preparation Strategy

#### Uniform Superposition Breakthrough
Our key innovation: uniform superposition states achieve optimal performance:

```python
def create_uniform_superposition(n_qubits: int) -> GeneralCircuitQuantumState:
    """Create uniform superposition |++...+⟩ state."""
    circuit = QuantumCircuit(n_qubits)
    for i in range(n_qubits):
        circuit.add_H_gate(i)  # Apply Hadamard to all qubits
    return GeneralCircuitQuantumState(n_qubits, circuit)
```

**Benefits:**
- Machine precision accuracy (< 1e-15)
- Complete subspace coverage when `num_states_pick_out = 2**n_qubits`
- Simple state preparation using only Hadamard gates
- Robust across all QSCI variants

### Time Evolution Implementation

#### Trotter Decomposition Integration
```python
from quri_algo.core.time_evolution import TrotterTimeEvolutionCircuitFactory

def _create_trotter_time_evolution(self, initial_state, evolution_time):
    """Create time evolution using Trotter decomposition."""
    time_evolution_factory = TrotterTimeEvolutionCircuitFactory(
        hamiltonian=self.hamiltonian,
        n_trotter_steps=self.num_trotter_steps
    )
    
    evolved_circuit = time_evolution_factory.create(evolution_time)
    return evolved_circuit(initial_state)
```

#### Exact Unitary Evolution
```python
def _create_exact_time_evolution(self, initial_state, evolution_time):
    """Create exact unitary time evolution."""
    # Convert to state vector, apply exact evolution, convert back
    state_vector = initial_state.to_state_vector()
    evolution_operator = expm(-1j * evolution_time * hamiltonian_matrix)
    evolved_vector = evolution_operator @ state_vector
    return create_circuit_state_from_vector(evolved_vector)
```

## QURI Integration Details

### Operator Handling
```python
# Native QURI Parts operator support
from quri_parts.core.operator import Operator, pauli_label

# Create operators using QURI Parts syntax
hamiltonian = Operator({
    pauli_label("Z0 Z1"): -1.0,
    pauli_label("X0"): -0.5,
    pauli_label("X1"): -0.5
})
```

### Sampler Integration
```python
# Use QURI Parts samplers
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

sampler = create_qulacs_vector_concurrent_sampler()
```

### Circuit Analysis with QURI VM
```python
from quri_parts.quri_vm import QVMCircuitResourceAnalyzer

def analyze_circuit_resources(circuit, architecture="star"):
    """Analyze circuit with QURI VM."""
    analyzer = QVMCircuitResourceAnalyzer(architecture)
    analysis = analyzer.analyze(circuit)
    
    return {
        'gate_count': analysis.gate_count,
        'circuit_depth': analysis.circuit_depth,
        'estimated_fidelity': analysis.circuit_fidelity,
        'execution_time': analysis.latency
    }
```

## Algorithm Variant Design

### VanillaQSCI
Standard QSCI with computational basis measurement:

```python
class VanillaQSCI(QSCIBase):
    def run(self, input_states, total_shots, **kwargs):
        # 1. Prepare measurement circuit
        # 2. Sample computational basis
        # 3. Apply QSCI selection
        # 4. Compute eigenvalues
        return QSCIResult(...)
```

### SingleTimeTE_QSCI
Time evolution at a single time point:

```python
class SingleTimeTE_QSCI(TimeEvolvedQSCI):
    def run(self, input_states, total_shots, **kwargs):
        # 1. Time evolve initial states
        # 2. Apply VanillaQSCI to evolved states
        # 3. Return ground state energy
        return QSCIResult(...)
```

### TimeAverageTE_QSCI
Time evolution averaged over multiple time points:

```python
class TimeAverageTE_QSCI(TimeEvolvedQSCI):
    def run(self, input_states, total_shots, **kwargs):
        # 1. Generate time points
        # 2. Evolve and measure at each time
        # 3. Average results
        # 4. Apply QSCI selection
        return QSCIResult(...)
```

### StateVectorTE_QSCI
Direct state vector processing for validation:

```python
class StateVectorTE_QSCI(TimeEvolvedQSCI):
    def run(self, input_states, total_shots, **kwargs):
        # 1. Direct state vector time evolution
        # 2. Exact probability computation
        # 3. QSCI selection and eigenvalue computation
        return QSCIResult(...)
```

## Performance Optimizations

### Sparse Matrix Operations
```python
import scipy.sparse as sp

def optimize_hamiltonian_operations(hamiltonian_matrix):
    """Convert to optimal sparse format."""
    if not sp.issparse(hamiltonian_matrix):
        hamiltonian_matrix = sp.csr_matrix(hamiltonian_matrix)
    return hamiltonian_matrix.tocsr()  # CSR for fast matrix-vector products
```

### Concurrent Sampling
```python
from concurrent.futures import ThreadPoolExecutor

def parallel_measurement_sampling(circuits, sampler, shots_per_circuit):
    """Parallelize measurement sampling."""
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(sampler.sample, circuit, shots_per_circuit)
            for circuit in circuits
        ]
        return [future.result() for future in futures]
```

### Memory Management
```python
def memory_efficient_eigenvalue_computation(selected_states, hamiltonian):
    """Compute eigenvalues with minimal memory footprint."""
    # Use iterative methods for large matrices
    if len(selected_states) > 1000:
        from scipy.sparse.linalg import eigsh
        return eigsh(hamiltonian_subspace, k=1, which='SA')[0]
    else:
        # Direct diagonalization for small matrices
        return np.linalg.eigvals(hamiltonian_subspace).min()
```

## Testing Architecture

### Uniform Superposition Testing Framework
```python
def create_test_framework(algorithm_class, test_systems):
    """Comprehensive testing with uniform superposition."""
    for system_name, system_data in test_systems.items():
        # Create uniform superposition
        uniform_state = create_uniform_superposition(system_data['n_qubits'])
        
        # Run with complete subspace coverage
        algorithm = algorithm_class(
            hamiltonian=system_data['hamiltonian'],
            sampler=create_qulacs_vector_concurrent_sampler(),
            num_states_pick_out=2**system_data['n_qubits']
        )
        
        result = algorithm.run([uniform_state], total_shots=2500)
        
        # Validate machine precision
        energy_error = abs(result.ground_state_energy - system_data['exact_energy'])
        assert energy_error < 1e-8, f"Failed precision test for {system_name}"
```

## Error Handling Strategy

### Comprehensive Error Checking
```python
class QSCIValidationError(Exception):
    """Custom exception for QSCI validation errors."""
    pass

def validate_algorithm_inputs(hamiltonian, sampler, num_states_pick_out):
    """Validate all algorithm inputs."""
    if not isinstance(hamiltonian, Operator):
        raise TypeError("Hamiltonian must be a QURI Parts Operator")
    
    if num_states_pick_out is not None and num_states_pick_out <= 0:
        raise ValueError("num_states_pick_out must be positive")
    
    # Additional validation logic...
```

## Future Design Considerations

### Extensibility Points
1. **New Algorithm Variants**: Inherit from `QSCIBase` or `TimeEvolvedQSCI`
2. **Custom Time Evolution**: Implement new time evolution methods
3. **Alternative Samplers**: Support for different quantum backends
4. **Advanced Analysis**: Extended QURI VM integration

### Scalability Considerations
1. **Large Molecule Support**: Memory-efficient handling of large Hamiltonians
2. **Distributed Computing**: Parallel execution across multiple nodes
3. **Hardware Integration**: Real quantum device support
4. **Error Mitigation**: Built-in error correction strategies

## Related Documentation

- **[Architecture](architecture.md)** - Detailed architectural diagrams and component relationships
- **[QURI Integration](quri_integration.md)** - Deep dive into QURI ecosystem integration
- **[Time Evolution](time_evolution.md)** - Comprehensive time evolution implementation details
- **[Performance](performance.md)** - Performance optimization strategies and benchmarks
- **[Future Plans](future.md)** - Roadmap and future development directions