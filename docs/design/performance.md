# Performance Considerations

The TE-QSCI implementation is designed for scalability and efficiency, with optimizations at multiple levels from quantum circuit compilation to classical post-processing. This page details the performance characteristics and optimization strategies.

## Computational Complexity Analysis

### QSCI Scaling

The computational complexity of QSCI algorithms depends on several key factors:

- **Configuration Selection**: O(N log N) where N is the total number of measured configurations
- **Hamiltonian Matrix Construction**: O(R²) where R is the number of selected configurations  
- **Eigenvalue Problem**: O(R³) for dense matrices, O(R²) for sparse matrices with iterative solvers

```python
def analyze_qsci_complexity(
    n_qubits: int,
    num_states_pick_out: int,
    hamiltonian_terms: int
) -> Dict[str, int]:
    """Analyze computational complexity for QSCI."""
    
    total_configurations = 2**n_qubits
    selected_configurations = min(num_states_pick_out, total_configurations)
    
    return {
        "configuration_selection": total_configurations * np.log(total_configurations),
        "hamiltonian_construction": selected_configurations**2 * hamiltonian_terms,
        "eigenvalue_solver": selected_configurations**3,
        "total_classical": (
            selected_configurations**3 + 
            selected_configurations**2 * hamiltonian_terms
        )
    }
```

### TE-QSCI Scaling

Time-evolved QSCI introduces additional complexity from quantum circuit simulation:

- **Time Evolution Circuit**: O(T × G) where T is Trotter steps, G is the number of Hamiltonian terms
- **Quantum Sampling**: O(S × D) where S is total shots, D is circuit depth
- **Multiple Time Points**: Linear scaling with the number of evolution times

```python
def analyze_te_qsci_complexity(
    n_qubits: int,
    hamiltonian_terms: int,
    evolution_time: float,
    trotter_steps: int,
    total_shots: int
) -> Dict[str, int]:
    """Analyze computational complexity for TE-QSCI."""
    
    # Circuit construction complexity
    circuit_gates = hamiltonian_terms * trotter_steps
    circuit_depth = estimate_circuit_depth(hamiltonian_terms, trotter_steps)
    
    # Sampling complexity (depends on backend)
    sampling_cost = total_shots * circuit_depth * n_qubits
    
    return {
        "circuit_construction": circuit_gates,
        "quantum_sampling": sampling_cost,
        "circuit_depth": circuit_depth,
        "time_evolution_overhead": circuit_gates / hamiltonian_terms
    }
```

## Memory Optimization

### Sparse Matrix Representations

For large molecular systems, sparse matrix representations provide significant memory savings:

```python
def _build_sparse_hamiltonian(
    self,
    configurations: Mapping[int, int]
) -> sparse.csr_matrix:
    """Build sparse Hamiltonian matrix for large systems."""
    
    config_list = list(configurations.keys())
    n_configs = len(config_list)
    
    # Pre-allocate arrays for sparse matrix construction
    max_nonzeros = n_configs * min(n_configs, 1000)  # Estimate sparsity
    row_indices = np.zeros(max_nonzeros, dtype=np.int32)
    col_indices = np.zeros(max_nonzeros, dtype=np.int32)
    data = np.zeros(max_nonzeros, dtype=np.complex128)
    
    nnz = 0  # Number of non-zeros
    
    for i, config_i in enumerate(config_list):
        for j, config_j in enumerate(config_list):
            matrix_element = self._compute_matrix_element(config_i, config_j)
            
            if abs(matrix_element) > 1e-12:  # Sparsity threshold
                if nnz >= max_nonzeros:
                    # Resize arrays if needed
                    row_indices = np.resize(row_indices, max_nonzeros * 2)
                    col_indices = np.resize(col_indices, max_nonzeros * 2)
                    data = np.resize(data, max_nonzeros * 2)
                    max_nonzeros *= 2
                
                row_indices[nnz] = i
                col_indices[nnz] = j
                data[nnz] = matrix_element
                nnz += 1
    
    # Trim arrays to actual size
    return sparse.csr_matrix(
        (data[:nnz], (row_indices[:nnz], col_indices[:nnz])), 
        shape=(n_configs, n_configs)
    )
```

**Memory Benefits:**
- Typical sparsity: 1-10% for molecular Hamiltonians
- Memory reduction: 10-100x for large systems
- Faster matrix-vector operations in iterative eigensolvers

### Efficient State Sampling

Adaptive sampling strategies reduce unnecessary quantum circuit evaluations:

```python
def _efficient_sampling_with_convergence(
    self,
    state: GeneralCircuitQuantumState,
    total_shots: int,
    convergence_threshold: float = 1e-3
) -> Mapping[int, int]:
    """Efficient sampling with convergence-based early termination."""
    
    # Adaptive sampling parameters
    initial_shots = min(1000, total_shots // 4)
    batch_size = 500
    remaining_shots = total_shots - initial_shots
    
    # Initial sampling
    counts = self._sample_state(state, initial_shots)
    previous_distribution = self._normalize_counts(counts)
    
    # Continue sampling until convergence or shot limit
    while remaining_shots > 0:
        # Sample additional batch
        current_batch = min(batch_size, remaining_shots)
        new_counts = self._sample_state(state, current_batch)
        
        # Merge counts
        for config, count in new_counts.items():
            counts[config] = counts.get(config, 0) + count
        
        # Check convergence
        current_distribution = self._normalize_counts(counts)
        if self._distributions_converged(
            previous_distribution, 
            current_distribution, 
            convergence_threshold
        ):
            break
        
        previous_distribution = current_distribution
        remaining_shots -= current_batch
    
    return counts

def _distributions_converged(
    self,
    dist1: Dict[int, float],
    dist2: Dict[int, float],
    threshold: float
) -> bool:
    """Check if two probability distributions have converged."""
    
    all_configs = set(dist1.keys()) | set(dist2.keys())
    
    # Calculate total variation distance
    tv_distance = 0.5 * sum(
        abs(dist1.get(config, 0) - dist2.get(config, 0))
        for config in all_configs
    )
    
    return tv_distance < threshold
```

## Parallelization Strategies

### Time Point Parallelization

For time-average TE-QSCI, different evolution times can be sampled in parallel:

```python
def run_time_average_parallel(
    self,
    initial_state: GeneralCircuitQuantumState,
    evolution_times: Sequence[float],
    total_shots: int,
    max_workers: Optional[int] = None
) -> QSCIResult:
    """Parallel time-average TE-QSCI implementation."""
    
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    if max_workers is None:
        max_workers = min(len(evolution_times), 4)  # Reasonable default
    
    shots_per_time = total_shots // len(evolution_times)
    
    def sample_evolution_time(t: float) -> Mapping[int, int]:
        """Sample a single evolution time."""
        evolved_state = self._create_time_evolved_state(initial_state, t)
        return self._sample_state(evolved_state, shots_per_time)
    
    # Execute sampling in parallel
    all_counts = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_time = {
            executor.submit(sample_evolution_time, t): t 
            for t in evolution_times
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_time):
            t = future_to_time[future]
            try:
                counts = future.result()
                
                # Merge counts from this time point
                for config, count in counts.items():
                    all_counts[config] = all_counts.get(config, 0) + count
                    
            except Exception as e:
                print(f"Error sampling time {t}: {e}")
                raise
    
    return self._solve_qsci_problem(all_counts)
```

### Configuration Processing Parallelization

Matrix element computation can be parallelized for large configuration spaces:

```python
def _compute_hamiltonian_parallel(
    self,
    configurations: List[int],
    max_workers: Optional[int] = None
) -> sparse.csr_matrix:
    """Parallel Hamiltonian matrix construction."""
    
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    
    n_configs = len(configurations)
    
    if max_workers is None:
        max_workers = min(mp.cpu_count(), 4)
    
    # Divide work into chunks
    chunk_size = max(1, n_configs // (max_workers * 4))
    row_chunks = [
        configurations[i:i + chunk_size] 
        for i in range(0, n_configs, chunk_size)
    ]
    
    def compute_matrix_chunk(row_configs: List[int]) -> List[Tuple[int, int, complex]]:
        """Compute matrix elements for a chunk of rows."""
        elements = []
        
        for i, config_i in enumerate(row_configs):
            for j, config_j in enumerate(configurations):
                element = self._compute_matrix_element(config_i, config_j)
                if abs(element) > 1e-12:
                    elements.append((i, j, element))
        
        return elements
    
    # Process chunks in parallel
    all_elements = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(compute_matrix_chunk, chunk) 
            for chunk in row_chunks
        ]
        
        row_offset = 0
        for future, chunk in zip(futures, row_chunks):
            elements = future.result()
            
            # Adjust row indices for global matrix
            adjusted_elements = [
                (row + row_offset, col, val) 
                for row, col, val in elements
            ]
            all_elements.extend(adjusted_elements)
            row_offset += len(chunk)
    
    # Build sparse matrix from elements
    if all_elements:
        rows, cols, data = zip(*all_elements)
        return sparse.csr_matrix(
            (data, (rows, cols)), 
            shape=(n_configs, n_configs)
        )
    else:
        return sparse.csr_matrix((n_configs, n_configs))
```

## Quantum Circuit Optimization

### Gate Count Reduction

Optimizations to reduce the number of quantum gates:

```python
def _optimize_pauli_evolution_circuit(
    self,
    circuit: QuantumCircuit
) -> QuantumCircuit:
    """Optimize Pauli evolution circuit for reduced gate count."""
    
    optimized_circuit = QuantumCircuit(circuit.qubit_count)
    
    # Group gates by type for optimization
    gate_groups = self._group_gates_by_type(circuit.gates)
    
    # Combine adjacent rotation gates
    for gate_type, gates in gate_groups.items():
        if gate_type in ["RX", "RY", "RZ"]:
            combined_gates = self._combine_rotation_gates(gates)
            for gate in combined_gates:
                optimized_circuit.add_gate(gate)
        else:
            # Add non-rotation gates as-is
            for gate in gates:
                optimized_circuit.add_gate(gate)
    
    return optimized_circuit

def _combine_rotation_gates(
    self,
    rotation_gates: List[Any]
) -> List[Any]:
    """Combine adjacent rotation gates on the same qubit."""
    
    # Group by qubit and rotation type
    qubit_rotations = {}
    for gate in rotation_gates:
        qubit = gate.target_indices[0]
        gate_type = gate.name
        key = (qubit, gate_type)
        
        if key not in qubit_rotations:
            qubit_rotations[key] = []
        qubit_rotations[key].append(gate)
    
    # Combine rotations for each (qubit, type) pair
    combined_gates = []
    for (qubit, gate_type), gates in qubit_rotations.items():
        total_angle = sum(gate.angle for gate in gates)
        
        # Only add gate if total angle is significant
        if abs(total_angle) > 1e-12:
            combined_gate = self._create_rotation_gate(gate_type, qubit, total_angle)
            combined_gates.append(combined_gate)
    
    return combined_gates
```

### Circuit Depth Optimization

Strategies to reduce circuit depth for better hardware performance:

```python
def _optimize_circuit_depth(
    self,
    circuit: QuantumCircuit
) -> QuantumCircuit:
    """Optimize circuit depth by reordering gates."""
    
    # Build dependency graph
    dependencies = self._build_gate_dependencies(circuit.gates)
    
    # Schedule gates to minimize depth
    scheduled_gates = self._schedule_gates_min_depth(dependencies)
    
    # Build optimized circuit
    optimized_circuit = QuantumCircuit(circuit.qubit_count)
    for gate in scheduled_gates:
        optimized_circuit.add_gate(gate)
    
    return optimized_circuit

def _schedule_gates_min_depth(
    self,
    dependencies: Dict[int, List[int]]
) -> List[Any]:
    """Schedule gates to minimize circuit depth using critical path method."""
    
    # Implement critical path scheduling
    # This is a simplified version - real implementation would be more complex
    
    ready_gates = [
        gate_id for gate_id, deps in dependencies.items() 
        if not deps
    ]
    scheduled = []
    completed = set()
    
    while ready_gates:
        # Schedule all ready gates in parallel
        current_layer = ready_gates[:]
        ready_gates = []
        
        for gate_id in current_layer:
            scheduled.append(self.gates[gate_id])
            completed.add(gate_id)
            
            # Check if any new gates become ready
            for other_gate_id, deps in dependencies.items():
                if (other_gate_id not in completed and 
                    other_gate_id not in ready_gates and
                    all(dep_id in completed for dep_id in deps)):
                    ready_gates.append(other_gate_id)
    
    return scheduled
```

## Performance Benchmarking

### Comprehensive Performance Analysis

```python
@dataclass
class PerformanceBenchmark:
    """Performance benchmark results for QSCI algorithms."""
    
    algorithm_name: str
    system_size: int
    execution_time: float
    memory_usage: float
    convergence_iterations: int
    final_accuracy: float
    
    def performance_score(self) -> float:
        """Calculate overall performance score."""
        # Weighted combination of metrics
        time_score = 1.0 / (1.0 + self.execution_time)
        memory_score = 1.0 / (1.0 + self.memory_usage / 1e6)  # MB
        accuracy_score = self.final_accuracy
        
        return (time_score + memory_score + accuracy_score) / 3.0

def benchmark_qsci_algorithms(
    hamiltonians: List[Operator],
    algorithms: List[str],
    shot_counts: List[int]
) -> List[PerformanceBenchmark]:
    """Comprehensive benchmarking of QSCI algorithms."""
    
    results = []
    
    for hamiltonian in hamiltonians:
        for algorithm_name in algorithms:
            for shots in shot_counts:
                
                # Setup algorithm
                algorithm = create_qsci_algorithm(
                    QSCIVariant[algorithm_name],
                    hamiltonian,
                    num_states_pick_out=200
                )
                
                # Benchmark execution
                start_time = time.time()
                start_memory = psutil.Process().memory_info().rss
                
                result = algorithm.run(
                    [self._create_initial_state(hamiltonian.qubit_count)],
                    shots
                )
                
                end_time = time.time()
                end_memory = psutil.Process().memory_info().rss
                
                # Calculate accuracy
                exact_energy = self._get_exact_ground_energy(hamiltonian)
                accuracy = 1.0 - abs(result.ground_state_energy - exact_energy) / abs(exact_energy)
                
                benchmark = PerformanceBenchmark(
                    algorithm_name=algorithm_name,
                    system_size=hamiltonian.qubit_count,
                    execution_time=end_time - start_time,
                    memory_usage=end_memory - start_memory,
                    convergence_iterations=result.execution_metadata.get('iterations', 1),
                    final_accuracy=accuracy
                )
                
                results.append(benchmark)
    
    return results
```

## Scaling Analysis

### Theoretical Scaling Behavior

```python
def analyze_scaling_behavior(
    system_sizes: List[int],
    algorithm_variants: List[str]
) -> Dict[str, Dict[str, Any]]:
    """Analyze theoretical and empirical scaling behavior."""
    
    scaling_analysis = {}
    
    for variant in algorithm_variants:
        
        # Theoretical complexity
        if variant == "VANILLA":
            time_complexity = "O(R³)"  # R = num_states_pick_out
            space_complexity = "O(R²)"
        elif variant in ["SINGLE_TIME_TE", "TIME_AVERAGE_TE"]:
            time_complexity = "O(T·G·S + R³)"  # T=trotter, G=gates, S=shots, R=configs
            space_complexity = "O(R²)"
        elif variant == "STATE_VECTOR_TE":
            time_complexity = "O(2^n)"  # n = number of qubits
            space_complexity = "O(2^n)"
        
        # Empirical measurements
        execution_times = []
        memory_usages = []
        
        for n_qubits in system_sizes:
            # Create test system
            hamiltonian = create_test_hamiltonian(n_qubits)
            
            # Measure performance
            performance = measure_algorithm_performance(hamiltonian, variant)
            execution_times.append(performance.execution_time)
            memory_usages.append(performance.memory_usage)
        
        # Fit scaling curves
        time_scaling = fit_polynomial(system_sizes, execution_times)
        memory_scaling = fit_polynomial(system_sizes, memory_usages)
        
        scaling_analysis[variant] = {
            "theoretical_time": time_complexity,
            "theoretical_space": space_complexity,
            "empirical_time_scaling": time_scaling,
            "empirical_memory_scaling": memory_scaling,
            "crossover_point": estimate_classical_quantum_crossover(variant)
        }
    
    return scaling_analysis
```

## Optimization Recommendations

### For Small Systems (< 10 qubits)
- Use StateVector TE-QSCI for exact results
- High Trotter steps for accurate time evolution
- Complete subspace sampling when feasible

### For Medium Systems (10-20 qubits)
- Single-time or time-average TE-QSCI
- Moderate Trotter steps (20-50)
- Sparse matrix representations
- Parallel time point sampling

### For Large Systems (> 20 qubits)
- Focus on most probable configurations
- Adaptive sampling strategies
- Aggressive circuit optimization
- Distributed computation for time averaging

```python
def get_optimization_strategy(
    n_qubits: int,
    available_memory: int,
    time_budget: float
) -> Dict[str, Any]:
    """Get recommended optimization strategy based on system constraints."""
    
    if n_qubits <= 10:
        return {
            "algorithm": "STATE_VECTOR_TE",
            "num_states_pick_out": 2**n_qubits,
            "trotter_steps": 50,
            "use_sparse": False,
            "parallel_sampling": False
        }
    elif n_qubits <= 20:
        return {
            "algorithm": "TIME_AVERAGE_TE",
            "num_states_pick_out": min(1000, 2**(n_qubits-2)),
            "trotter_steps": 30,
            "use_sparse": True,
            "parallel_sampling": True,
            "evolution_times": [0.5, 1.0, 1.5]
        }
    else:
        return {
            "algorithm": "SINGLE_TIME_TE",
            "num_states_pick_out": min(500, available_memory // (8 * n_qubits)),
            "trotter_steps": 20,
            "use_sparse": True,
            "parallel_sampling": True,
            "adaptive_sampling": True,
            "circuit_optimization": True
        }
```

## Navigation

- **Previous**: [Time Evolution](time_evolution.md) - Time evolution implementation details
- **Next**: [Algorithm Overview](../algorithms/overview.md) - Complete algorithm descriptions  
- **Related**: [Architecture](architecture.md) - System architecture and design
- **See Also**: [QURI Integration](quri_integration.md) - Ecosystem integration benefits