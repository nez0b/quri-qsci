# Getting Started with VM-Aware QSCI

This guide walks you through setting up and running your first VM-aware TE-QSCI calculation with authentic time evolution and quantum virtual machine analysis.

## Prerequisites

Ensure you have the required dependencies installed:

```bash
# Core dependencies
uv add quri-parts quri-vm

# Optional: For authentic time evolution
uv add quri-algo

# Development dependencies
uv add pytest pytest-cov
```

## Basic Setup

### 1. Import Core Components

```python
from quri_parts.core.operator import Operator, pauli_label
from src.qsci_vm_interface import create_vm_aware_te_qsci
from src.qsci_vm_analysis import create_star_vm_for_h6, VMFactory
```

### 2. Create a Simple Hamiltonian

For demonstration, let's create a simple 2-qubit Hamiltonian:

```python
# Create a simple H2-like Hamiltonian
hamiltonian = Operator({
    pauli_label("Z0"): -1.0523732,
    pauli_label("Z1"): -1.0523732,
    pauli_label("Z0 Z1"): -0.39793742,
    pauli_label("X0 X1"): -0.39793742,
    pauli_label("Y0 Y1"): -0.39793742,
})

print(f"Hamiltonian has {len(hamiltonian.mapping)} terms")
```

### 3. Configure VM Interface

Create a VM interface with STAR architecture:

```python
# Create VM interface with low error rate for testing
error_rate = 1e-5
vm_interface = create_star_vm_for_h6(error_rate)

print(f"VM interface configured with error rate: {error_rate}")
```

### 4. Create VM-Aware Algorithm

```python
# Create the VM-aware TE-QSCI algorithm
algorithm = create_vm_aware_te_qsci(
    hamiltonian=hamiltonian,
    evolution_time=1.0,
    trotter_steps=4,
    vm_interface=vm_interface,
    use_vm_sampling=True,
    sampling_level="ArchLogicalCircuit"  # Use architecture-aware analysis
)

print(f"Algorithm: {algorithm.name}")
```

### 5. Run the Algorithm

```python
# Execute the VM-aware calculation
result = algorithm.run()

# Access results
print("=== VM-Aware TE-QSCI Results ===")
print(f"Ground state energy: {result.qsci_result.ground_state_energy}")
print(f"Selected energies: {result.qsci_result.selected_energies}")

# VM Analysis results
if result.vm_analysis:
    logical_analysis = result.vm_analysis.logical_analysis
    arch_analysis = result.vm_analysis.arch_analysis
    
    print(f"\n=== Circuit Analysis ===")
    print(f"Logical circuit - Gates: {logical_analysis.gate_count}, Depth: {logical_analysis.circuit_depth}")
    print(f"Arch circuit - Gates: {arch_analysis.gate_count}, Depth: {arch_analysis.circuit_depth}")
    
    if hasattr(arch_analysis, 'latency') and arch_analysis.latency:
        print(f"Circuit latency: {arch_analysis.latency}")
    if hasattr(arch_analysis, 'circuit_fidelity') and arch_analysis.circuit_fidelity:
        print(f"Circuit fidelity: {arch_analysis.circuit_fidelity}")
```

## Complete Example

Here's a complete working example you can run:

```python
"""
VM-Aware QSCI Getting Started Example
"""
from quri_parts.core.operator import Operator, pauli_label
from src.qsci_vm_interface import create_vm_aware_te_qsci
from src.qsci_vm_analysis import create_star_vm_for_h6

def main():
    # 1. Create simple H2 Hamiltonian
    hamiltonian = Operator({
        pauli_label("Z0"): -1.0523732,
        pauli_label("Z1"): -1.0523732,
        pauli_label("Z0 Z1"): -0.39793742,
        pauli_label("X0 X1"): -0.39793742,
        pauli_label("Y0 Y1"): -0.39793742,
    })
    
    # 2. Configure VM with low error rate
    vm_interface = create_star_vm_for_h6(error_rate=1e-5)
    
    # 3. Create VM-aware algorithm
    algorithm = create_vm_aware_te_qsci(
        hamiltonian=hamiltonian,
        evolution_time=1.0,
        trotter_steps=4,
        vm_interface=vm_interface,
        use_vm_sampling=True
    )
    
    # 4. Run calculation
    print(f"Running {algorithm.name}...")
    result = algorithm.run()
    
    # 5. Display results
    print("=== Results ===")
    print(f"Ground state energy: {result.qsci_result.ground_state_energy:.6f}")
    print(f"Selected energies: {[f'{e:.6f}' for e in result.qsci_result.selected_energies]}")
    
    if result.vm_analysis:
        logical = result.vm_analysis.logical_analysis
        arch = result.vm_analysis.arch_analysis
        print(f"\nLogical Circuit: {logical.gate_count} gates, depth {logical.circuit_depth}")
        print(f"Arch Circuit: {arch.gate_count} gates, depth {arch.circuit_depth}")

if __name__ == "__main__":
    main()
```

## Understanding the Output

When you run the example, you'll see output similar to:

```
Running VMAwareSingleTimeTE_QSCI...
=== Results ===
Ground state energy: -1.857275
Selected energies: ['-1.857275', '-1.242458', '-0.882722', '-0.882722']

Logical Circuit: 24 gates, depth 12
Arch Circuit: 96 gates, depth 28
```

This shows:
- **Ground state energy**: The computed lowest energy eigenvalue
- **Selected energies**: Multiple energy levels computed by QSCI
- **Logical Circuit**: Gate count and depth of the ideal quantum circuit
- **Arch Circuit**: Resource requirements when mapped to STAR architecture

## Error Rate Impact

You can study the impact of quantum hardware errors by comparing different error rates:

```python
error_rates = [0.0, 1e-5, 1e-4, 1e-3]

for error_rate in error_rates:
    vm_interface = create_star_vm_for_h6(error_rate)
    algorithm = create_vm_aware_te_qsci(
        hamiltonian=hamiltonian,
        evolution_time=1.0,
        vm_interface=vm_interface
    )
    
    result = algorithm.run()
    print(f"Error rate {error_rate}: Energy = {result.qsci_result.ground_state_energy:.6f}")
```

## Next Steps

Now that you have a basic VM-aware QSCI calculation running:

1. **[Authentic Time Evolution](authentic_time_evolution.md)** - Learn about real vs dummy time evolution circuits
2. **[VM Analysis](vm_analysis.md)** - Dive deeper into VM metrics and resource estimation  
3. **[Molecular Systems](molecular_systems.md)** - Work with realistic quantum chemistry problems
4. **[API Reference](../../api/vm_qsci.md)** - Explore all available functionality

## Troubleshooting

### Common Issues

**ImportError for quri-algo**:
If you see import errors related to quri-algo, the algorithm will automatically fall back to simplified Trotter decomposition. For authentic time evolution, install quri-algo:
```bash
uv add quri-algo
```

**VM Creation Errors**:
Ensure quri-vm is properly installed and STAR device properties are available:
```bash
uv add quri-vm
```

**Performance Issues**:
For large molecular systems, start with smaller molecules (H2, H4) before moving to H6 or larger systems.