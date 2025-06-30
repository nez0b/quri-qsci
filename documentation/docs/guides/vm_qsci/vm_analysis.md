# VM Analysis and Resource Estimation

The VM-Aware QSCI module provides comprehensive quantum virtual machine analysis capabilities, enabling realistic assessment of quantum algorithm resource requirements and performance characteristics.

## Overview

VM analysis bridges the gap between ideal quantum algorithm design and realistic quantum hardware implementation by providing detailed resource estimation, performance metrics, and architecture-specific analysis.

## Analysis Levels

### LogicalCircuit Analysis
Analysis of the ideal quantum circuit without hardware-specific considerations:

- **Gate Count**: Total number of quantum gates
- **Circuit Depth**: Maximum number of sequential gate operations
- **Qubit Count**: Number of qubits required
- **Circuit Structure**: Gate types and connectivity patterns

### ArchLogicalCircuit Analysis  
Analysis considering specific quantum hardware architecture (STAR):

- **Transpiled Gate Count**: Gates after basis set conversion
- **Architecture Depth**: Circuit depth on target hardware
- **Circuit Latency**: Total execution time including gate times
- **Circuit Fidelity**: Expected fidelity after noise and errors
- **Resource Overhead**: Additional resources for error correction

## VM Interface Components

### QURIVMInterface

Central interface for VM analysis providing unified access to different analysis levels:

```python
from src.qsci_vm_analysis import QURIVMInterface, VMFactory

# Create VM interface with STAR and ideal VMs
vm_interface = QURIVMInterface(
    ideal_vm=VMFactory.create_ideal_vm(),
    star_vm=VMFactory.create_star_vm(error_rate=1e-3)
)

# Analyze circuit at logical level
logical_analysis = vm_interface.analyze_circuit_at_level(circuit, "LogicalCircuit")

# Analyze circuit at architecture level
arch_analysis = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")
```

### VMFactory

Factory for creating pre-configured VM instances:

```python
# Create ideal VM for baseline analysis
ideal_vm = VMFactory.create_ideal_vm()

# Create STAR VM with specific error rate
star_vm = VMFactory.create_star_vm(error_rate=1e-4)

# Create H6-optimized STAR VM
h6_vm = VMFactory.create_star_vm_for_h6(error_rate=1e-3)
```

### VMSampler

Quantum circuit sampling using VM instances:

```python
from src.qsci_vm_analysis import VMSampler

# Create sampler with STAR VM
sampler = VMSampler(star_vm)

# Sample quantum circuits
circuit_shot_pairs = [(circuit1, 1000), (circuit2, 2000)]
results = sampler(circuit_shot_pairs)
```

## Analysis Results

### VMAnalysisResult Structure

```python
@dataclass
class VMAnalysisResult:
    logical_analysis: Any      # LogicalCircuit analysis
    arch_analysis: Any         # ArchLogicalCircuit analysis 
    circuit_description: str   # Circuit metadata
    analysis_time: float       # Analysis timing
```

### Key Metrics

**Circuit Resource Metrics**:
- `gate_count`: Total quantum gates
- `circuit_depth`: Sequential gate depth
- `qubit_count`: Required qubits

**Architecture-Specific Metrics** (when available):
- `latency`: Circuit execution time
- `circuit_fidelity`: Expected output fidelity
- `error_budget`: Cumulative error estimates

## Practical Examples

### Basic Circuit Analysis

```python
from quri_parts.circuit import QuantumCircuit
from quri_parts.circuit.gate import RX, RY, CNOT
from src.qsci_vm_analysis import create_star_vm_for_h6

# Create sample circuit
circuit = QuantumCircuit(4)
circuit.add_gate(RX(0, 0.5))
circuit.add_gate(RY(1, 0.3))
circuit.add_gate(CNOT(0, 1))
circuit.add_gate(RX(2, 0.7))

# Create VM interface
vm_interface = create_star_vm_for_h6(error_rate=1e-3)

# Perform analysis
logical_result = vm_interface.analyze_circuit_at_level(circuit, "LogicalCircuit")
arch_result = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")

print("=== Logical Circuit Analysis ===")
print(f"Gates: {logical_result.gate_count}")
print(f"Depth: {logical_result.circuit_depth}")
print(f"Qubits: {logical_result.qubit_count}")

print("\n=== Architecture Circuit Analysis ===")
print(f"Gates: {arch_result.gate_count}")
print(f"Depth: {arch_result.circuit_depth}")
if hasattr(arch_result, 'latency') and arch_result.latency:
    print(f"Latency: {arch_result.latency}")
if hasattr(arch_result, 'circuit_fidelity') and arch_result.circuit_fidelity:
    print(f"Fidelity: {arch_result.circuit_fidelity}")
```

### Multi-Error Rate Analysis

```python
def analyze_error_impact(circuit, error_rates):
    """Analyze circuit performance across different error rates"""
    results = {}
    
    for error_rate in error_rates:
        vm_interface = create_star_vm_for_h6(error_rate)
        analysis = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")
        
        results[error_rate] = {
            'gates': analysis.gate_count,
            'depth': analysis.circuit_depth,
            'fidelity': getattr(analysis, 'circuit_fidelity', None),
            'latency': getattr(analysis, 'latency', None)
        }
    
    return results

# Analyze across error rates
error_rates = [1e-5, 1e-4, 1e-3, 1e-2]
results = analyze_error_impact(your_circuit, error_rates)

for error_rate, metrics in results.items():
    print(f"Error rate {error_rate}:")
    print(f"  Gates: {metrics['gates']}, Depth: {metrics['depth']}")
    if metrics['fidelity']:
        print(f"  Fidelity: {metrics['fidelity']:.4f}")
```

### VM Metrics Calculation

The VM-aware algorithms automatically calculate performance metrics:

```python
from src.qsci_vm_interface import create_vm_aware_te_qsci

# Create algorithm with VM analysis
algorithm = create_vm_aware_te_qsci(
    hamiltonian=hamiltonian,
    evolution_time=1.0,
    vm_interface=vm_interface,
    use_vm_sampling=True
)

# Run with automatic metrics calculation
result = algorithm.run()

# Access VM metrics
if hasattr(result, 'vm_overhead') and result.vm_overhead:
    print(f"VM Overhead: {result.vm_overhead:.2f}x")

if hasattr(result, 'fidelity_impact') and result.fidelity_impact:
    print(f"Fidelity Impact: {result.fidelity_impact:.4f}")
```

## Performance Considerations

### Analysis Caching

The VM analysis implements intelligent caching to avoid redundant computations:

```python
from src.qsci_vm_interface import VMAwareTimeEvolutionCircuit

# Circuit wrapper with analysis caching
wrapped_circuit = VMAwareTimeEvolutionCircuit(circuit, vm_interface)

# First analysis - performs computation
analysis1 = wrapped_circuit.analyze_at_level("LogicalCircuit")

# Second analysis - uses cached result  
analysis2 = wrapped_circuit.analyze_at_level("LogicalCircuit")
assert analysis1 is analysis2  # Same object reference
```

### Resource Optimization

For large circuits and molecular systems:

```python
# Use appropriate VM configurations for system size
def create_optimized_vm_interface(n_qubits, error_rate):
    if n_qubits <= 4:  # Small systems
        return QURIVMInterface(
            ideal_vm=VMFactory.create_ideal_vm(),
            star_vm=VMFactory.create_star_vm(error_rate)
        )
    elif n_qubits <= 12:  # Medium systems like H6
        return create_star_vm_for_h6(error_rate)
    else:  # Large systems
        # Could implement specialized configurations
        return create_star_vm_for_h6(error_rate)
```

### Batch Analysis

Analyze multiple circuits efficiently:

```python
def batch_analyze_circuits(circuits, vm_interface, level="ArchLogicalCircuit"):
    """Batch analysis of multiple circuits"""
    results = []
    
    for circuit in circuits:
        try:
            analysis = vm_interface.analyze_circuit_at_level(circuit, level)
            results.append({
                'circuit': circuit,
                'analysis': analysis,
                'gates': analysis.gate_count,
                'depth': analysis.circuit_depth
            })
        except Exception as e:
            print(f"Analysis failed for circuit: {e}")
            results.append({'circuit': circuit, 'error': str(e)})
    
    return results
```

## Advanced Analysis Features

### Circuit Comparison

Compare different quantum circuits or algorithms:

```python
def compare_circuits(circuit1, circuit2, vm_interface):
    """Compare resource requirements of two circuits"""
    analysis1 = vm_interface.analyze_circuit_at_level(circuit1, "ArchLogicalCircuit")
    analysis2 = vm_interface.analyze_circuit_at_level(circuit2, "ArchLogicalCircuit")
    
    comparison = {
        'gate_ratio': analysis2.gate_count / analysis1.gate_count,
        'depth_ratio': analysis2.circuit_depth / analysis1.circuit_depth,
    }
    
    if (hasattr(analysis1, 'latency') and analysis1.latency and 
        hasattr(analysis2, 'latency') and analysis2.latency):
        comparison['latency_ratio'] = analysis2.latency / analysis1.latency
    
    return comparison
```

### Resource Scaling Analysis

Study how resources scale with problem size:

```python
def analyze_scaling(hamiltonian_generator, qubit_counts, vm_interface):
    """Analyze resource scaling with system size"""
    scaling_data = []
    
    for n_qubits in qubit_counts:
        hamiltonian = hamiltonian_generator(n_qubits)
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=1.0,
            vm_interface=vm_interface
        )
        
        # Generate circuit for analysis
        circuit = algorithm._generate_authentic_time_evolution_circuit(n_qubits)
        analysis = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")
        
        scaling_data.append({
            'qubits': n_qubits,
            'gates': analysis.gate_count,
            'depth': analysis.circuit_depth,
            'hamiltonian_terms': len(hamiltonian.mapping)
        })
    
    return scaling_data
```

## Error Handling and Debugging

### VM Analysis Errors

Common issues and solutions:

```python
def robust_vm_analysis(circuit, vm_interface, level="ArchLogicalCircuit"):
    """Robust VM analysis with error handling"""
    try:
        return vm_interface.analyze_circuit_at_level(circuit, level)
    except Exception as e:
        print(f"VM analysis failed: {e}")
        
        # Try fallback to logical analysis
        if level == "ArchLogicalCircuit":
            print("Falling back to LogicalCircuit analysis")
            return vm_interface.analyze_circuit_at_level(circuit, "LogicalCircuit")
        else:
            raise
```

### Performance Monitoring

Monitor analysis performance:

```python
import time

def timed_analysis(circuit, vm_interface, level):
    """Analyze circuit with timing information"""
    start_time = time.time()
    analysis = vm_interface.analyze_circuit_at_level(circuit, level)
    end_time = time.time()
    
    return {
        'analysis': analysis,
        'analysis_time': end_time - start_time,
        'gates_per_second': analysis.gate_count / (end_time - start_time)
    }
```

## Integration with QSCI Algorithms

VM analysis integrates seamlessly with QSCI algorithms:

```python
# The VM-aware algorithm automatically performs analysis
result = algorithm.run()

# Access comprehensive results
print("=== QSCI Results ===")
print(f"Ground state energy: {result.qsci_result.ground_state_energy}")

print("\n=== VM Analysis ===")
if result.vm_analysis:
    logical = result.vm_analysis.logical_analysis
    arch = result.vm_analysis.arch_analysis
    print(f"Resource overhead: {arch.gate_count / logical.gate_count:.1f}x")
```

## Next Steps

- **[Molecular Systems](molecular_systems.md)** - Apply VM analysis to quantum chemistry
- **[Getting Started](getting_started.md)** - Basic usage examples
- **[API Reference](../../api/vm_qsci.md)** - Complete technical documentation