# API Reference

This section provides detailed API documentation for all QSCI modules and classes.

## Core Modules

### QSCI Algorithms
- **[qsci_algorithms](qsci_algorithms/)** - Core QSCI algorithm implementations
  - `VanillaQSCI` - Standard QSCI algorithm
  - `TimeEvolvedQSCI` - Base class for TE-QSCI variants
  - `SingleTimeTE_QSCI` - Single time evolution QSCI
  - `TimeAverageTE_QSCI` - Time-averaged QSCI
  - `StateVectorTE_QSCI` - State vector TE-QSCI

### Algorithm Interfaces
- **[qsci_algo_interface](qsci_algo_interface/)** - QURI Algo compatible interfaces
  - Algorithm wrappers for quri-algo integration
  - Standard interfaces for all QSCI variants
  - Result handling and data structures

### Probability Calculator
- **[probability_calculator](probability_calculator/)** - Probability computation utilities
  - State probability calculations
  - Measurement outcome processing
  - Statistical analysis tools

### QURI VM Analysis
- **[qsci_vm_analysis](qsci_vm_analysis/)** - Circuit analysis and resource estimation
  - LogicalCircuit analysis
  - ArchLogicalCircuit analysis with STAR architecture
  - Gate count and depth estimation
  - Quantum resource analysis

## Quick Navigation

| Module | Description | Key Classes |
|--------|-------------|-------------|
| [qsci_algorithms](qsci_algorithms/) | Core algorithm implementations | `VanillaQSCI`, `TimeEvolvedQSCI` |
| [qsci_algo_interface](qsci_algo_interface/) | QURI Algo interfaces | Algorithm wrappers and interfaces |
| [probability_calculator](probability_calculator/) | Probability utilities | Probability computation functions |
| [qsci_vm_analysis](qsci_vm_analysis/) | Circuit analysis | QURI VM integration tools |

## Usage Examples

### Basic Algorithm Usage

```python
from qsci_algorithms import VanillaQSCI
from quri_parts.core.operator import Operator
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Create QSCI algorithm instance
qsci = VanillaQSCI(
    hamiltonian=your_hamiltonian,
    sampler=create_qulacs_vector_concurrent_sampler(),
    num_states_pick_out=16
)

# Run algorithm
result = qsci.run(initial_states, total_shots=2000)
print(f"Ground state energy: {result.ground_state_energy}")
```

### Time-Evolved QSCI

```python
from qsci_algorithms import SingleTimeTE_QSCI

# Create TE-QSCI algorithm
te_qsci = SingleTimeTE_QSCI(
    hamiltonian=your_hamiltonian,
    sampler=sampler,
    evolution_time=1.0,
    num_trotter_steps=10
)

result = te_qsci.run(initial_states, total_shots=3000)
```

### QURI Algo Integration

```python
from qsci_algo_interface import QSCIGroundState

# Use with quri-algo
algorithm = QSCIGroundState(
    hamiltonian=hamiltonian,
    variant="vanilla",
    sampler=sampler
)

result = algorithm(initial_state)
```

### Circuit Analysis

```python
from qsci_vm_analysis import analyze_qsci_circuit

# Analyze circuit resources
analysis = analyze_qsci_circuit(
    circuit=your_circuit,
    architecture="star",
    analysis_level="arch_logical"
)

print(f"Gate count: {analysis.gate_count}")
print(f"Circuit depth: {analysis.circuit_depth}")
print(f"Estimated fidelity: {analysis.circuit_fidelity}")
```

## Type Annotations

All public functions and classes include comprehensive type annotations for better IDE support and documentation clarity.

## Error Handling

The API includes proper error handling with descriptive error messages. Common exceptions:

- `ValueError` - Invalid parameter values
- `TypeError` - Incorrect parameter types  
- `RuntimeError` - Algorithm execution errors
- `NotImplementedError` - Unimplemented features

## Performance Considerations

- Use `create_qulacs_vector_concurrent_sampler()` for best performance
- Set appropriate `num_states_pick_out` values (typically `2**n_qubits` for small systems)
- Consider time evolution parameters for TE-QSCI variants
- Use circuit analysis for resource planning

For detailed documentation of specific functions and classes, see the individual module pages linked above.