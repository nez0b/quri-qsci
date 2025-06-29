# Algorithm Implementation Overview

This page provides a comprehensive overview of the four QSCI algorithm variants implemented in this framework, each designed for specific use cases and computational requirements.

## QSCI Variants Summary

The implementation provides four distinct QSCI algorithm variants:

| Variant | Use Case | Time Evolution | Key Features |
|---------|----------|----------------|--------------|
| **VanillaQSCI** | Standard QSCI baseline | None | Computational basis sampling, proven convergence |
| **SingleTimeTE_QSCI** | Systematic configuration exploration | Single time point | Fixed evolution time, deterministic results |
| **TimeAverageTE_QSCI** | Enhanced sampling diversity | Multiple time points | Reduced statistical fluctuations |
| **StateVectorTE_QSCI** | Exact validation and small systems | Direct state vector | No sampling noise, exact results |

## Core Algorithm Architecture

### QSCIBase: Abstract Foundation

All QSCI algorithms inherit from the abstract `QSCIBase` class, which defines the common interface and shared functionality:

```python
class QSCIBase(ABC):
    """Abstract base class for all QSCI algorithms."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None
    ):
        """Initialize QSCI algorithm.
        
        Args:
            hamiltonian: Target Hamiltonian to diagonalize
            sampler: Quantum sampler for measurement
            num_states_pick_out: Number of states to select for subspace
        """
        if not is_hermitian(hamiltonian):
            raise ValueError("Hamiltonian must be Hermitian")
        
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        self.num_states_pick_out = num_states_pick_out
    
    @abstractmethod
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """Run the QSCI algorithm."""
        pass
```

**Common Methods:**
- `_pick_out_states()`: Select most frequent states from measurement counts
- `_generate_truncated_hamiltonian()`: Build Hamiltonian matrix in selected subspace
- `_diagonalize_truncated_hamiltonian()`: Solve eigenvalue problem using scipy
- `_calculate_computational_basis_probabilities()`: Calculate exact state probabilities

## 1. VanillaQSCI: Standard Implementation

The baseline QSCI algorithm that uses computational basis state sampling from input quantum states.

### Class Definition

```python
class VanillaQSCI(QSCIBase):
    """Vanilla QSCI algorithm implementation.
    
    Uses computational basis sampling from input states to build
    the configuration subspace for eigenvalue decomposition.
    """
    
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """Run vanilla QSCI algorithm."""
        
        # Sample all input states
        all_counts = {}
        shots_per_state = total_shots // len(input_states)
        
        for state in input_states:
            counts = self._sample_state(state, shots_per_state)
            for config, count in counts.items():
                all_counts[config] = all_counts.get(config, 0) + count
        
        # Select most probable configurations
        selected_configs = self._pick_out_states(all_counts)
        
        # Build and diagonalize truncated Hamiltonian
        return self._solve_qsci_problem(selected_configs)
```

### Implementation Details

**Sampling Strategy:**
- Direct computational basis measurement from input quantum states
- Equal shot distribution across multiple input states
- Configuration selection based on measurement frequency

**Matrix Construction:**
```python
def _build_hamiltonian_matrix(
    self,
    configurations: Mapping[int, int]
) -> np.ndarray:
    """Build Hamiltonian matrix in the selected configuration subspace."""
    
    config_list = list(configurations.keys())
    n_configs = len(config_list)
    hamiltonian_matrix = np.zeros((n_configs, n_configs), dtype=complex)
    
    for i, config_i in enumerate(config_list):
        for j, config_j in enumerate(config_list):
            # Calculate matrix element ⟨config_i|Ĥ|config_j⟩
            matrix_element = self._compute_matrix_element(config_i, config_j)
            hamiltonian_matrix[i, j] = matrix_element
    
    return hamiltonian_matrix
```

**Use Cases:**
- Baseline comparison for TE-QSCI variants
- Systems where good initial state guesses are available
- Validation of QSCI methodology

## 2. TimeEvolvedQSCI: TE-QSCI Base Implementation

The foundational class for all time-evolved variants, providing three execution modes.

### Core Methods

#### Single-Time Evolution
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
    
    Evolves the initial state for a fixed time and samples
    the resulting configuration distribution.
    """
    
    # Create time-evolved quantum state
    evolved_state = self._create_time_evolved_state(
        initial_state, evolution_time, trotter_steps
    )
    
    # Sample evolved state
    measurement_counts = self._sample_state(evolved_state, total_shots)
    
    # Solve QSCI in the sampled subspace
    return self._solve_qsci_problem(measurement_counts)
```

#### Time-Average Evolution
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
    
    Samples multiple evolution times and averages the
    resulting configuration distributions.
    """
    
    all_counts = {}
    
    for t in evolution_times:
        # Evolve at this time point
        evolved_state = self._create_time_evolved_state(
            initial_state, t, trotter_steps
        )
        
        # Sample and accumulate counts
        counts = self._sample_state(evolved_state, shots_per_time)
        for config, count in counts.items():
            all_counts[config] = all_counts.get(config, 0) + count
    
    return self._solve_qsci_problem(all_counts)
```

#### State Vector Evolution
```python
def run_state_vector(
    self,
    initial_state: QuantumState,
    evolution_time: float,
    num_eigenstates: int = 1,
    **kwargs
) -> QSCIResult:
    """Run TE-QSCI with direct state vector calculation.
    
    Computes exact probabilities without sampling noise.
    """
    
    # Create time evolution circuit
    evolution_circuit = self._create_time_evolved_circuit(
        initial_state.circuit, evolution_time
    )
    
    # Calculate exact state vector probabilities
    probabilities = self._calculate_state_vector_probabilities(evolution_circuit)
    
    # Convert probabilities to configuration counts
    configurations = self._probabilities_to_configurations(probabilities)
    
    return self._solve_qsci_problem(configurations)
```

## 3. SingleTimeTE_QSCI: Fixed Evolution Time

Wrapper class for time evolution at a single, fixed evolution time.

### Class Implementation

```python
class SingleTimeTE_QSCI(TimeEvolvedQSCI):
    """Single-time TE-QSCI wrapper for testing compatibility.
    
    Provides a simplified interface for TE-QSCI at a fixed
    evolution time, ideal for systematic parameter studies.
    """
    
    def __init__(
        self, 
        hamiltonian: Operator, 
        sampler: ConcurrentSampler, 
        evolution_time: float, 
        num_states_pick_out: Optional[int] = None
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(
        self, 
        input_states: Sequence[CircuitQuantumState], 
        total_shots: int, 
        **kwargs
    ) -> QSCIResult:
        """Run single-time TE-QSCI at the fixed evolution time."""
        
        if len(input_states) != 1:
            raise ValueError("SingleTimeTE_QSCI expects exactly one initial state")
        
        return self.run_single_time(
            input_states[0], 
            self.evolution_time, 
            total_shots, 
            **kwargs
        )
```

### Scientific Rationale

**Time Evolution Theory:**
The time evolution operator generates configurations through its Taylor expansion:

```
|ψ(t)⟩ = e^(-iĤt)|ψ₀⟩ = Σₖ (-iĤt)ᵏ/k! |ψ₀⟩
```

The k-th order term naturally includes up to 2k-excitation configurations, providing a systematic exploration of the electronic configuration space.

**Optimal Evolution Time:**
```python
def estimate_optimal_evolution_time(
    hamiltonian: Operator,
    target_excitation_level: int = 2
) -> float:
    """Estimate optimal evolution time for target excitation level."""
    
    # Estimate Hamiltonian norm
    hamiltonian_norm = estimate_operator_norm(hamiltonian)
    
    # Time should be large enough to generate target excitations
    # but small enough to avoid over-rotation
    optimal_time = np.pi / (2 * hamiltonian_norm)
    
    return optimal_time
```

**Use Cases:**
- Parameter optimization studies
- Systematic analysis of evolution time effects
- Reproducible research with fixed parameters

## 4. TimeAverageTE_QSCI: Multi-Time Averaging

Wrapper class for averaging over multiple evolution times to improve sampling diversity.

### Class Implementation

```python
class TimeAverageTE_QSCI(TimeEvolvedQSCI):
    """Time-average TE-QSCI wrapper for testing compatibility.
    
    Averages configurations over multiple evolution times
    to reduce statistical fluctuations and improve sampling diversity.
    """
    
    def __init__(
        self, 
        hamiltonian: Operator, 
        sampler: ConcurrentSampler, 
        evolution_times: Sequence[float], 
        num_states_pick_out: Optional[int] = None
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_times = evolution_times
    
    def run(
        self, 
        input_states: Sequence[CircuitQuantumState], 
        total_shots: int, 
        **kwargs
    ) -> QSCIResult:
        """Run time-average TE-QSCI over the specified evolution times."""
        
        if len(input_states) != 1:
            raise ValueError("TimeAverageTE_QSCI expects exactly one initial state")
        
        shots_per_time = total_shots // len(self.evolution_times)
        
        return self.run_time_average(
            input_states[0], 
            self.evolution_times, 
            shots_per_time, 
            **kwargs
        )
```

### Statistical Benefits

**Variance Reduction:**
Time averaging reduces statistical fluctuations in configuration sampling:

```python
def analyze_variance_reduction(
    single_time_results: List[QSCIResult],
    time_average_result: QSCIResult
) -> float:
    """Analyze variance reduction from time averaging."""
    
    # Calculate variance in single-time results
    single_time_energies = [result.ground_state_energy for result in single_time_results]
    single_time_variance = np.var(single_time_energies)
    
    # Time-averaged result has lower variance
    # (theoretical reduction factor is 1/√N for N time points)
    n_times = len(single_time_results)
    theoretical_reduction = 1.0 / np.sqrt(n_times)
    
    return theoretical_reduction
```

**Evolution Time Selection:**
```python
def select_evolution_times(
    hamiltonian: Operator,
    n_times: int = 4,
    time_range: Tuple[float, float] = (0.5, 2.0)
) -> List[float]:
    """Select optimal evolution times for time averaging."""
    
    # Use logarithmic spacing for better coverage
    t_min, t_max = time_range
    times = np.logspace(np.log10(t_min), np.log10(t_max), n_times)
    
    return times.tolist()
```

**Use Cases:**
- Improved statistical reliability
- Exploration of multiple time scales
- Robust results across different evolution regimes

## 5. StateVectorTE_QSCI: Exact Simulation

Wrapper class for exact state vector calculation without sampling noise.

### Class Implementation

```python
class StateVectorTE_QSCI(TimeEvolvedQSCI):
    """State vector TE-QSCI wrapper for testing compatibility.
    
    Uses direct state vector calculation for exact results
    without sampling noise. Suitable for small systems and validation.
    """
    
    def __init__(
        self, 
        hamiltonian: Operator, 
        sampler: ConcurrentSampler, 
        evolution_time: float, 
        num_states_pick_out: Optional[int] = None
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(
        self, 
        input_states: Sequence[CircuitQuantumState], 
        total_shots: int, 
        **kwargs
    ) -> QSCIResult:
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
            initial_state, 
            self.evolution_time, 
            **kwargs
        )
```

### Exact Calculation Methods

**State Vector Computation:**
```python
def _calculate_exact_state_vector(
    self,
    evolution_circuit: QuantumCircuit
) -> np.ndarray:
    """Calculate exact state vector after time evolution."""
    
    from quri_parts.qulacs.simulator import evaluate_state_to_vector
    
    # Create quantum state from circuit
    n_qubits = evolution_circuit.qubit_count
    state = GeneralCircuitQuantumState(n_qubits, evolution_circuit)
    
    # Calculate exact state vector
    state_vector = evaluate_state_to_vector(state)
    
    return state_vector

def _extract_significant_configurations(
    self,
    state_vector: np.ndarray,
    threshold: float = 1e-6
) -> Mapping[int, int]:
    """Extract configurations with significant probability amplitudes."""
    
    configurations = {}
    
    for config, amplitude in enumerate(state_vector):
        probability = abs(amplitude) ** 2
        
        if probability > threshold:
            # Convert probability to "count" for QSCI interface
            # Use large reference count for high precision
            count = int(probability * 1e6)
            if count > 0:
                configurations[config] = count
    
    return configurations
```

**Use Cases:**
- Small system validation (≤ 10 qubits)
- Algorithm development and testing
- Exact benchmark results
- Noise-free performance analysis

## Algorithm Selection Guide

### Decision Matrix

| System Size | Goal | Recommended Algorithm | Key Parameters |
|-------------|------|----------------------|----------------|
| ≤ 6 qubits | Exact validation | StateVectorTE_QSCI | Complete subspace |
| 6-12 qubits | High accuracy | TimeAverageTE_QSCI | 4-6 time points |
| 12-20 qubits | Good performance | SingleTimeTE_QSCI | Optimized evolution time |
| > 20 qubits | Scalability test | VanillaQSCI | Multiple initial states |

### Configuration Examples

#### Small System (H₂ molecule)
```python
# H2 molecule - exact validation
algorithm = StateVectorTE_QSCI(
    hamiltonian=h2_hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_states_pick_out=16  # Complete 4-qubit space
)

result = algorithm.run([hf_state], total_shots=1000)
```

#### Medium System (H₆ chain)
```python
# H6 chain - balanced accuracy and performance
evolution_times = [0.5, 1.0, 1.5, 2.0]
algorithm = TimeAverageTE_QSCI(
    hamiltonian=h6_hamiltonian,
    sampler=sampler,
    evolution_times=evolution_times,
    num_states_pick_out=200
)

result = algorithm.run([hf_state], total_shots=2000)
```

#### Large System (optimization)
```python
# Large system - performance focused
algorithm = SingleTimeTE_QSCI(
    hamiltonian=large_hamiltonian,
    sampler=sampler,
    evolution_time=1.2,  # Optimized value
    num_states_pick_out=100
)

result = algorithm.run([hf_state], total_shots=5000)
```

## Performance Characteristics

### Computational Complexity

| Algorithm | Time Complexity | Space Complexity | Quantum Resources |
|-----------|----------------|------------------|-------------------|
| VanillaQSCI | O(R³) | O(R²) | Input state preparation |
| SingleTimeTE_QSCI | O(T·G + R³) | O(R²) | Time evolution circuit |
| TimeAverageTE_QSCI | O(N·T·G + R³) | O(R²) | N × time evolution |
| StateVectorTE_QSCI | O(2ⁿ + R³) | O(2ⁿ) | Exact simulation |

Where:
- R = num_states_pick_out (selected configurations)
- T = trotter_steps (time evolution accuracy)
- G = number of Hamiltonian terms
- N = number of evolution times
- n = number of qubits

### Scaling Recommendations

```python
def get_algorithm_recommendation(
    n_qubits: int,
    hamiltonian_terms: int,
    accuracy_target: float,
    time_budget: float
) -> Tuple[str, Dict[str, Any]]:
    """Get algorithm recommendation based on system characteristics."""
    
    if n_qubits <= 8 and accuracy_target < 1e-6:
        return "StateVectorTE_QSCI", {
            "evolution_time": 1.0,
            "num_states_pick_out": 2**n_qubits
        }
    
    elif n_qubits <= 15 and time_budget > 60:  # 1 minute
        return "TimeAverageTE_QSCI", {
            "evolution_times": [0.5, 1.0, 1.5, 2.0],
            "num_states_pick_out": min(500, 2**(n_qubits-2))
        }
    
    else:
        return "SingleTimeTE_QSCI", {
            "evolution_time": estimate_optimal_evolution_time(hamiltonian_terms),
            "num_states_pick_out": min(200, time_budget * 10)
        }
```

## Navigation

- **Previous**: [Performance](../design/performance.md) - Performance analysis and optimization
- **Related**: [Time Evolution](../design/time_evolution.md) - Time evolution implementation details
- **See Also**: [Architecture](../design/architecture.md) - System architecture overview
- **Examples**: Browse the `/examples` directory for usage demonstrations