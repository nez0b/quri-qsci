# QURI Ecosystem Integration

The TE-QSCI implementation is deeply integrated with the QURI ecosystem, leveraging QURI Parts for quantum operations, QURI Algo for algorithm interfaces, and QURI VM for resource analysis.

## QURI Parts Integration

### Operator Representation

We use QURI Parts' `Operator` class for comprehensive Hamiltonian representation, providing native Pauli string support and efficient operator arithmetic.

```python
from quri_parts.core.operator import Operator, pauli_label

# Proper operator construction
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})
hamiltonian += Operator({pauli_label("Z1"): -1.25})
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})
hamiltonian.constant = -2.5  # Identity term
```

**Key Benefits:**
- Native Pauli string representation for quantum chemistry
- Efficient operator arithmetic for Hamiltonian construction
- Built-in Hermiticity checking for physical validity
- Sparse representation optimized for large molecular systems

### Quantum State Management

The framework uses `GeneralCircuitQuantumState` for quantum state representation, enabling circuit-based state preparation and automatic state vector computation.

```python
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit

# Hartree-Fock state preparation
def create_hf_state(n_qubits: int, n_electrons: int) -> GeneralCircuitQuantumState:
    circuit = QuantumCircuit(n_qubits)
    for i in range(n_electrons):
        circuit.add_X_gate(i)
    return GeneralCircuitQuantumState(n_qubits, circuit)
```

**Key Benefits:**
- Circuit-based state representation for quantum hardware compatibility
- Automatic state vector computation for classical simulation
- Unified sampling interface across different backends
- Full compatibility with QURI Parts backend ecosystem

### Sampling and Measurement

Integration with QURI Parts sampling framework provides efficient quantum measurement simulation:

```python
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
from quri_parts.core.sampling import create_sampler_from_concurrent_sampler

# Sampler setup for quantum measurements
concurrent_sampler = create_qulacs_vector_concurrent_sampler()
sampler = create_sampler_from_concurrent_sampler(concurrent_sampler)

# Use in TE-QSCI measurements
measurement_counts = sampler(evolved_state, n_shots)
```

**Features:**
- High-performance concurrent sampling using Qulacs backend
- Support for multiple quantum backends (Qulacs, Cirq, etc.)
- Optimized measurement collection for large shot counts
- Statistical analysis tools for measurement data

## QURI Algo Integration

### Algorithm Interface Compliance

Our implementation fully complies with quri-algo's standardized `Algorithm` interface, ensuring compatibility with the broader algorithm ecosystem.

```python
from quri_parts.algo.interface import Algorithm, AlgorithmResult

class QSCIAlgorithmBase(Algorithm):
    @property
    def name(self) -> str:
        """Algorithm name for identification."""
        
    def run(self, *args, **kwargs) -> QSCIAlgorithmResult:
        """Execute the algorithm and return results."""
```

**Key Benefits:**
- Standardized algorithm interface for tool compatibility
- Consistent result type system across algorithms
- Integration with quri-algo analysis framework
- Benchmarking and comparison capabilities

### Factory Pattern Implementation

Factory functions provide easy algorithm creation with sensible defaults:

```python
def create_qsci_algorithm(
    variant: QSCIVariant,
    hamiltonian: Operator,
    **kwargs
) -> QSCIAlgorithmBase:
    """Factory function for creating QSCI algorithms."""
    
    if variant == QSCIVariant.VANILLA:
        return VanillaQSCIAlgorithm(hamiltonian, **kwargs)
    elif variant == QSCIVariant.SINGLE_TIME_TE:
        return SingleTimeTeQSCIAlgorithm(hamiltonian, **kwargs)
    elif variant == QSCIVariant.TIME_AVERAGE_TE:
        return TimeAverageTeQSCIAlgorithm(hamiltonian, **kwargs)
    elif variant == QSCIVariant.STATE_VECTOR_TE:
        return StateVectorTeQSCIAlgorithm(hamiltonian, **kwargs)
    else:
        raise ValueError(f"Unknown QSCI variant: {variant}")
```

### Analysis Framework Integration

Custom analysis classes provide comprehensive algorithm performance evaluation:

```python
class QSCIAnalysis(Analysis):
    """Analysis specific to QSCI algorithms."""
    
    def analyze_convergence(self) -> Dict[str, Any]:
        """Analyze algorithm convergence properties."""
        return {
            "convergence_rate": self._calculate_convergence_rate(),
            "energy_variance": self._calculate_energy_variance(),
            "configuration_stability": self._analyze_configuration_stability()
        }
        
    def analyze_resource_requirements(self) -> Dict[str, Any]:
        """Analyze computational resource requirements."""
        return {
            "quantum_resources": self._estimate_quantum_resources(),
            "classical_overhead": self._estimate_classical_overhead(),
            "scaling_behavior": self._analyze_scaling_behavior()
        }
```

## QURI VM Integration

### Circuit Analysis at Multiple Levels

QURI VM integration provides comprehensive circuit analysis at different abstraction levels:

```python
class VMCircuitAnalyzer:
    def analyze_logical_circuit(
        self, 
        circuit: QuantumCircuit
    ) -> CircuitResourceEstimate:
        """Analyze circuit at logical level."""
        gate_count = self._count_gates_by_type(circuit)
        circuit_depth = circuit.depth()
        
        return CircuitResourceEstimate(
            gate_count=gate_count,
            circuit_depth=circuit_depth,
            logical_qubit_count=circuit.qubit_count
        )
        
    def analyze_arch_logical_circuit(
        self, 
        circuit: QuantumCircuit,
        architecture: str = "STAR"
    ) -> CircuitResourceEstimate:
        """Analyze circuit with architecture constraints."""
        logical_estimate = self.analyze_logical_circuit(circuit)
        
        # Apply architecture-specific overheads
        swap_overhead = self._estimate_swap_overhead(circuit, architecture)
        physical_qubits = self._estimate_physical_qubits(circuit, architecture)
        
        return CircuitResourceEstimate(
            gate_count=logical_estimate.gate_count,
            circuit_depth=logical_estimate.circuit_depth + swap_overhead,
            logical_qubit_count=logical_estimate.logical_qubit_count,
            physical_qubit_count=physical_qubits,
            execution_time=self._estimate_execution_time(circuit, architecture),
            fidelity_estimate=self._estimate_fidelity(circuit, architecture)
        )
```

**Analysis Metrics:**
- Gate count by type (single-qubit, two-qubit, multi-qubit)
- Circuit depth and critical path analysis
- Execution time estimates for different hardware
- Fidelity estimates including noise models
- Physical qubit requirements for specific architectures

### Architecture-Specific Optimizations

The system includes optimizations for common quantum architectures:

```python
def _estimate_star_swap_overhead(
    self, 
    circuit: QuantumCircuit
) -> int:
    """Estimate SWAP gate overhead for STAR topology."""
    two_qubit_gates = sum(
        1 for gate in circuit.gates 
        if len(gate.target_indices) == 2
    )
    
    # STAR topology reduces SWAP overhead due to central connectivity
    swap_fraction = 0.1  # Optimized for STAR connectivity
    return int(two_qubit_gates * swap_fraction)

def _estimate_grid_swap_overhead(
    self, 
    circuit: QuantumCircuit
) -> int:
    """Estimate SWAP gate overhead for 2D grid topology."""
    two_qubit_gates = sum(
        1 for gate in circuit.gates 
        if len(gate.target_indices) == 2
    )
    
    # Grid topology requires more SWAPs for long-range connectivity
    swap_fraction = 0.3
    return int(two_qubit_gates * swap_fraction)
```

### Resource Estimation Framework

Comprehensive resource estimation provides detailed insights into algorithm requirements:

```python
@dataclass
class CircuitResourceEstimate:
    gate_count: Dict[str, int]           # Gates by type
    circuit_depth: int                   # Circuit depth
    physical_qubit_count: int            # Physical qubits needed
    logical_qubit_count: int             # Logical qubits
    execution_time: TimeValue            # Estimated execution time
    fidelity_estimate: float             # Expected fidelity
    
    def total_gate_count(self) -> int:
        """Total number of gates in the circuit."""
        return sum(self.gate_count.values())
    
    def two_qubit_gate_count(self) -> int:
        """Number of two-qubit gates (most resource-intensive)."""
        return sum(
            count for gate_type, count in self.gate_count.items()
            if gate_type in ["CNOT", "CZ", "RZZ", "RYY", "RXX"]
        )
```

## Integration Examples

### Simple Algorithm Creation

```python
from src import create_qsci_algorithm, QSCIVariant
from quri_parts.core.operator import Operator, pauli_label

# Create Hamiltonian
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.25})
hamiltonian += Operator({pauli_label("Z1"): -1.25})
hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})

# Create algorithm
algorithm = create_qsci_algorithm(
    QSCIVariant.SINGLE_TIME_TE,
    hamiltonian,
    evolution_time=1.0,
    num_states_pick_out=100
)

# Run algorithm
result = algorithm.run(initial_state, total_shots=1000)
```

### Advanced Analysis Workflow

```python
from src import create_vm_enabled_algorithm, LoweringLevel

# Create VM-enabled algorithm
vm_algorithm = create_vm_enabled_algorithm(algorithm, "STAR")

# Analyze at different levels
logical_analysis = vm_algorithm.analyze(
    [initial_state], 1000, LoweringLevel.LogicalCircuit
)

arch_analysis = vm_algorithm.analyze(
    [initial_state], 1000, LoweringLevel.ArchLogicalCircuit
)

# Compare resource requirements
print(f"Logical qubits: {logical_analysis.logical_qubit_count}")
print(f"Physical qubits: {arch_analysis.physical_qubit_count}")
print(f"Circuit depth increase: {arch_analysis.circuit_depth - logical_analysis.circuit_depth}")
print(f"Execution time: {arch_analysis.execution_time}")
```

## Benefits of QURI Integration

### 1. Ecosystem Compatibility
- Seamless integration with existing QURI tools and workflows
- Standard interfaces enable tool interoperability
- Consistent data formats across the ecosystem

### 2. Performance Optimization
- Optimized backend implementations (Qulacs, etc.)
- Efficient operator representations
- Hardware-aware circuit compilation

### 3. Research Enablement
- Access to cutting-edge quantum algorithms
- Comprehensive analysis and benchmarking tools
- Support for emerging quantum hardware architectures

### 4. Maintenance and Support
- Benefit from ecosystem-wide updates and improvements
- Community-driven development and testing
- Long-term sustainability through QURI foundation

## Navigation

- **Previous**: [Architecture](architecture.md) - System architecture overview
- **Next**: [Time Evolution](time_evolution.md) - Time evolution implementation
- **See Also**: [Performance](performance.md) - Performance considerations