# Time Evolution Implementation

The Time-Evolved Quantum Selected Configuration Interaction (TE-QSCI) algorithm leverages time evolution to systematically explore electronic configurations without requiring variational optimization. This page details the implementation of time evolution in our framework.

## Scientific Background

TE-QSCI generates important electronic configurations through time evolution of an initial state:

```
|ψ(t)⟩ = e^(-iĤt)|ψ_I⟩ = |ψ_I⟩ - iĤt|ψ_I⟩ + (−iĤt)²/2!|ψ_I⟩ + ...
```

Where the k-th order term includes up to 2k-th order excitations, systematically exploring the configuration space needed for accurate quantum chemistry calculations.

**Key Insight**: Time evolution naturally generates the most important electronic configurations through the Taylor expansion, providing a non-variational path to ground state preparation.

## TimeEvolvedQSCI Class

The `TimeEvolvedQSCI` class serves as the base implementation for all TE-QSCI variants, providing three distinct execution modes:

### Class Definition

```python
class TimeEvolvedQSCI(QSCIBase):
    """Time-Evolved QSCI (TE-QSCI) algorithm implementation."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        time_evolution_method: str = "trotter"
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.time_evolution_method = time_evolution_method
        self._validate_time_evolution_method()
```

### Three Execution Modes

#### 1. Single-Time Evolution

Evolution at a specific time point for systematic configuration exploration:

```python
def run_single_time(
    self,
    initial_state: CircuitQuantumState,
    evolution_time: float,
    total_shots: int,
    trotter_steps: Optional[int] = None,
    **kwargs
) -> QSCIResult:
    """Run single-time TE-QSCI algorithm.
    
    Args:
        initial_state: Initial quantum state (typically Hartree-Fock)
        evolution_time: Time for evolution (t in e^(-iĤt))
        total_shots: Number of measurement shots
        trotter_steps: Steps for Trotter decomposition
        
    Returns:
        QSCIResult with ground state energy and configurations
    """
    # Create time-evolved quantum state
    evolved_state = self._create_time_evolved_state(
        initial_state, evolution_time, trotter_steps
    )
    
    # Sample evolved state to get configurations
    measurement_counts = self._sample_state(evolved_state, total_shots)
    
    # Solve QSCI in the sampled subspace
    return self._solve_qsci_problem(measurement_counts)
```

#### 2. Time-Average Evolution

Evolution averaged over multiple time points to improve sampling diversity:

```python
def run_time_average(
    self,
    initial_state: CircuitQuantumState,
    evolution_times: Sequence[float],
    shots_per_time: int,
    trotter_steps: Optional[int] = None,
    **kwargs
) -> QSCIResult:
    """Run time-average TE-QSCI algorithm.
    
    Args:
        initial_state: Initial quantum state
        evolution_times: List of evolution times
        shots_per_time: Shots per time point
        trotter_steps: Steps for Trotter decomposition
        
    Returns:
        QSCIResult with averaged configurations
    """
    all_counts = {}
    
    # Collect measurements from all time points
    for t in evolution_times:
        evolved_state = self._create_time_evolved_state(
            initial_state, t, trotter_steps
        )
        counts = self._sample_state(evolved_state, shots_per_time)
        
        # Merge counts from this time point
        for config, count in counts.items():
            all_counts[config] = all_counts.get(config, 0) + count
    
    return self._solve_qsci_problem(all_counts)
```

#### 3. State Vector Evolution

Direct state vector calculation without sampling (for small systems):

```python
def run_state_vector(
    self,
    initial_state: QuantumState,
    evolution_time: float,
    num_eigenstates: int = 1,
    **kwargs
) -> QSCIResult:
    """Run TE-QSCI with direct state vector calculation.
    
    Args:
        initial_state: Initial quantum state
        evolution_time: Evolution time
        num_eigenstates: Number of eigenstates to compute
        
    Returns:
        QSCIResult with exact state vector results
    """
    # Create time evolution circuit
    evolution_circuit = self._create_time_evolved_circuit(
        initial_state.circuit, evolution_time
    )
    
    # Calculate state vector probabilities directly
    probabilities = self._calculate_state_vector_probabilities(
        evolution_circuit
    )
    
    # Use all significant configurations
    configurations = self._extract_significant_configurations(probabilities)
    
    return self._solve_qsci_problem(configurations)
```

## Time Evolution Implementation

### Trotter Decomposition

The core time evolution uses first-order Trotter decomposition to approximate the matrix exponential:

```python
def _create_trotter_circuit(
    self,
    initial_circuit: QuantumCircuit,
    evolution_time: float,
    trotter_steps: int
) -> QuantumCircuit:
    """Create Trotter-decomposed time evolution circuit."""
    
    evolved_circuit = initial_circuit.copy()
    dt = evolution_time / trotter_steps
    
    for step in range(trotter_steps):
        # First-order Trotter decomposition: e^(-iĤt) ≈ ∏ e^(-iĤ_j dt)
        for pauli_term, coeff in self.hamiltonian.items():
            if pauli_term.is_identity():
                continue  # Skip identity terms
                
            # Add exp(-i * coeff * dt * PauliTerm)
            self._add_pauli_evolution_gate(
                evolved_circuit, 
                pauli_term, 
                coeff * dt
            )
    
    return evolved_circuit
```

### Pauli Evolution Gates

Individual Pauli terms are evolved using rotation gates:

```python
def _add_pauli_evolution_gate(
    self,
    circuit: QuantumCircuit,
    pauli_term: PauliLabel,
    angle: float
) -> None:
    """Add Pauli evolution gate to circuit.
    
    Implements exp(-i * angle * PauliTerm) using rotations.
    """
    
    # Convert to Z-basis if needed
    pre_gates, post_gates = self._pauli_to_z_conversion(pauli_term)
    
    # Apply pre-rotation gates (X→Z, Y→Z conversions)
    for gate in pre_gates:
        circuit.add_gate(gate)
    
    # Apply controlled Z-rotations
    qubits = [idx for idx, _ in pauli_term]
    if len(qubits) == 1:
        # Single-qubit Z rotation
        circuit.add_RZ_gate(qubits[0], angle)
    else:
        # Multi-qubit Pauli evolution
        self._add_multi_qubit_pauli_evolution(circuit, qubits, angle)
    
    # Apply post-rotation gates (Z→X, Z→Y conversions)
    for gate in post_gates:
        circuit.add_gate(gate)

def _pauli_to_z_conversion(self, pauli_term: PauliLabel) -> Tuple[List, List]:
    """Convert Pauli operators to Z-basis for evolution."""
    pre_gates = []
    post_gates = []
    
    for qubit_idx, pauli_op in pauli_term:
        if pauli_op == "X":
            # X → Z: apply H before and after
            pre_gates.append(("H", qubit_idx))
            post_gates.append(("H", qubit_idx))
        elif pauli_op == "Y":
            # Y → Z: apply RX(π/2) before and RX(-π/2) after
            pre_gates.append(("RX", qubit_idx, np.pi/2))
            post_gates.append(("RX", qubit_idx, -np.pi/2))
        # Z operators need no conversion
    
    return pre_gates, list(reversed(post_gates))
```

### Multi-Qubit Pauli Evolution

For multi-qubit Pauli terms, we implement controlled evolution:

```python
def _add_multi_qubit_pauli_evolution(
    self,
    circuit: QuantumCircuit,
    qubits: List[int],
    angle: float
) -> None:
    """Add multi-qubit Pauli evolution using CNOT ladders."""
    
    if len(qubits) == 2:
        # Two-qubit case: use CNOT-RZ-CNOT pattern
        control, target = qubits[0], qubits[1]
        circuit.add_CNOT_gate(control, target)
        circuit.add_RZ_gate(target, angle)
        circuit.add_CNOT_gate(control, target)
    
    else:
        # Multi-qubit case: CNOT ladder approach
        # Forward CNOT ladder
        for i in range(len(qubits) - 1):
            circuit.add_CNOT_gate(qubits[i], qubits[i + 1])
        
        # Z rotation on last qubit
        circuit.add_RZ_gate(qubits[-1], angle)
        
        # Backward CNOT ladder
        for i in range(len(qubits) - 2, -1, -1):
            circuit.add_CNOT_gate(qubits[i], qubits[i + 1])
```

## Integration with QURI-Algo

When quri-algo is available, we leverage its optimized time evolution implementations:

```python
def _create_quri_algo_time_evolution(
    self,
    initial_circuit: QuantumCircuit,
    evolution_time: float,
    trotter_steps: Optional[int] = None
) -> QuantumCircuit:
    """Create time evolution using quri-algo factories."""
    
    try:
        from quri_algo.circuit.time_evolution.trotter_time_evo import (
            TrotterTimeEvolutionCircuitFactory
        )
        from quri_algo.problem import QubitHamiltonian
        
        # Convert to quri-algo format
        qubit_hamiltonian = self._convert_to_qubit_hamiltonian()
        
        # Create Trotter factory
        if trotter_steps is None:
            trotter_steps = max(10, int(evolution_time * 10))
        
        factory = TrotterTimeEvolutionCircuitFactory(
            hamiltonian=qubit_hamiltonian,
            delta_t=evolution_time / trotter_steps
        )
        
        # Generate time evolution circuit
        time_evo_circuit = factory.create_circuit(trotter_steps)
        
        # Combine with initial state circuit
        combined_circuit = initial_circuit.copy()
        combined_circuit.extend(time_evo_circuit)
        
        return combined_circuit
        
    except ImportError:
        # Fallback to internal implementation
        return self._create_trotter_circuit(
            initial_circuit, evolution_time, trotter_steps or 20
        )
```

## Algorithm Variants

The framework provides three wrapper classes for different TE-QSCI applications:

### SingleTimeTE_QSCI

For studies at a specific evolution time:

```python
class SingleTimeTE_QSCI(TimeEvolvedQSCI):
    """Single-time TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_time, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(self, input_states, total_shots, **kwargs):
        """Run single-time TE-QSCI at fixed evolution time."""
        if len(input_states) != 1:
            raise ValueError("SingleTimeTE_QSCI expects exactly one initial state")
        return self.run_single_time(
            input_states[0], self.evolution_time, total_shots, **kwargs
        )
```

### TimeAverageTE_QSCI

For averaging over multiple time points:

```python
class TimeAverageTE_QSCI(TimeEvolvedQSCI):
    """Time-average TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_times, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_times = evolution_times
    
    def run(self, input_states, total_shots, **kwargs):
        """Run time-average TE-QSCI over multiple evolution times."""
        if len(input_states) != 1:
            raise ValueError("TimeAverageTE_QSCI expects exactly one initial state")
        shots_per_time = total_shots // len(self.evolution_times)
        return self.run_time_average(
            input_states[0], self.evolution_times, shots_per_time, **kwargs
        )
```

### StateVectorTE_QSCI

For exact simulation without sampling noise:

```python
class StateVectorTE_QSCI(TimeEvolvedQSCI):
    """State vector TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_time, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(self, input_states, total_shots, **kwargs):
        """Run TE-QSCI with direct state vector calculation."""
        if len(input_states) != 1:
            raise ValueError("StateVectorTE_QSCI expects exactly one initial state")
        
        initial_state = input_states[0]
        
        # Validate circuit attribute requirement
        if not hasattr(initial_state, 'circuit'):
            raise TypeError(
                "StateVectorTE_QSCI requires a GeneralCircuitQuantumState with a 'circuit' attribute"
            )
        
        return self.run_state_vector(
            initial_state, self.evolution_time, **kwargs
        )
```

## Performance Optimizations

### Adaptive Trotter Steps

```python
def _adaptive_trotter_steps(
    self,
    evolution_time: float,
    target_error: float = 1e-3
) -> int:
    """Determine optimal Trotter steps for given error tolerance."""
    
    # Estimate Trotter error scaling
    hamiltonian_norm = self._estimate_hamiltonian_norm()
    
    # First-order Trotter error: O((Ĥt)²/n)
    # Solve for n: target_error = C * (Ĥt)² / n
    required_steps = int(np.ceil(
        (hamiltonian_norm * evolution_time) ** 2 / target_error
    ))
    
    # Ensure minimum and maximum bounds
    return max(10, min(required_steps, 100))
```

### Circuit Optimization

```python
def _optimize_evolution_circuit(self, circuit: QuantumCircuit) -> QuantumCircuit:
    """Optimize time evolution circuit for better performance."""
    
    # Combine adjacent rotation gates
    optimized_circuit = self._combine_rotation_gates(circuit)
    
    # Remove redundant CNOT pairs
    optimized_circuit = self._remove_redundant_cnots(optimized_circuit)
    
    # Optimize gate ordering for hardware constraints
    optimized_circuit = self._optimize_gate_ordering(optimized_circuit)
    
    return optimized_circuit
```

## Usage Examples

### Basic Single-Time Evolution

```python
from src import TimeEvolvedQSCI
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Setup
sampler = create_qulacs_vector_concurrent_sampler()
te_qsci = TimeEvolvedQSCI(hamiltonian, sampler)

# Run single-time evolution
result = te_qsci.run_single_time(
    initial_state=hf_state,
    evolution_time=1.0,
    total_shots=2000,
    trotter_steps=20
)

print(f"Ground state energy: {result.ground_state_energy}")
```

### Time-Average Evolution

```python
# Multiple time points for better sampling
evolution_times = [0.5, 1.0, 1.5, 2.0]

result = te_qsci.run_time_average(
    initial_state=hf_state,
    evolution_times=evolution_times,
    shots_per_time=500,  # Total: 2000 shots
    trotter_steps=20
)
```

### State Vector Evolution

```python
# Exact evolution for small systems
result = te_qsci.run_state_vector(
    initial_state=hf_state,
    evolution_time=1.0,
    num_eigenstates=3  # Ground + 2 excited states
)
```

## Navigation

- **Previous**: [QURI Integration](quri_integration.md) - Ecosystem integration details  
- **Next**: [Performance](performance.md) - Performance analysis and optimization
- **Related**: [Algorithm Overview](../algorithms/overview.md) - All QSCI variants
- **See Also**: [Architecture](architecture.md) - System architecture overview