# State Conversion Guide

The `state_conversion` module provides utilities for converting quantum states between different representations. This is essential for integrating ffsim-generated states with the QSCI framework and other QURI Parts tools.

## Overview

The state conversion module handles:

- **Format Conversion**: ffsim fermionic states ↔ QURI Parts computational basis
- **Jordan-Wigner Mapping**: Proper qubit ordering and fermionic-to-qubit transformation
- **Validation Metrics**: Fidelity, overlap, and normalization checks
- **Multiple Methods**: Different conversion strategies for various use cases

## Understanding State Representations

### ffsim Fermionic States

ffsim represents quantum states in the fermionic Fock space:
- States are vectors in the space of fermionic configurations
- Dimension depends on the number of orbitals and electrons
- Natural for quantum chemistry calculations

### QURI Parts Computational Basis

QURI Parts uses computational basis states on qubits:
- States are superpositions of computational basis states (|00⟩, |01⟩, |10⟩, |11⟩, etc.)
- Each qubit can be |0⟩ or |1⟩
- Required for most quantum circuit operations

## Basic State Conversion

### Simple Conversion Example

```python
from ffsim_integration.molecular_systems import create_h2_molecule
from ffsim_integration.integration import create_lucj_ansatz
from ffsim_integration.state_conversion import ffsim_to_quri_state

# Create molecular system and ansatz
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)
lucj_result = create_lucj_ansatz(h2, n_reps=1, max_iterations=50)

# Convert state
n_qubits = 2 * h2.ffsim_mol_data.norb
nelec = h2.ffsim_mol_data.nelec

conversion_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=5000
)

print(f"Conversion Results:")
print(f"  - Input dimension: {len(lucj_result.state_vector)}")
print(f"  - Output dimension: 2^{n_qubits} = {2**n_qubits}")
print(f"  - Conversion fidelity: {conversion_result.metrics.fidelity:.4f}")
print(f"  - State vector norm: {conversion_result.metrics.state_vector_norm:.6f}")
```

## Conversion Methods

### 1. Sampling Circuit Method

The most robust method for general use:

```python
# Sampling circuit conversion
sampling_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=10000  # More shots = better accuracy
)

print(f"Sampling Circuit Results:")
print(f"  - Fidelity: {sampling_result.metrics.fidelity:.4f}")
print(f"  - Probability overlap: {sampling_result.metrics.probability_overlap:.4f}")
print(f"  - Max probability difference: {sampling_result.metrics.max_probability_difference:.4f}")
```

### 2. Direct Mapping Method

For cases where direct mapping is possible:

```python
# Direct mapping (when electron configuration is known)
direct_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="direct_mapping"
)

print(f"Direct Mapping Results:")
print(f"  - Exact conversion: {direct_result.metrics.fidelity == 1.0}")
print(f"  - State type: {type(direct_result.quri_state)}")
```

## Understanding Conversion Metrics

### Fidelity

Measures how accurately the state was converted:

```python
# Analyze conversion fidelity
fidelity = conversion_result.metrics.fidelity

if fidelity > 0.95:
    print(f"Excellent conversion: fidelity = {fidelity:.4f}")
elif fidelity > 0.85:
    print(f"Good conversion: fidelity = {fidelity:.4f}")
elif fidelity > 0.70:
    print(f"Acceptable conversion: fidelity = {fidelity:.4f}")
else:
    print(f"Poor conversion: fidelity = {fidelity:.4f}")
    print("Consider increasing total_shots or using a different method")
```

### Probability Overlap

Compares probability distributions:

```python
# Check probability overlap
overlap = conversion_result.metrics.probability_overlap
max_diff = conversion_result.metrics.max_probability_difference

print(f"Probability Analysis:")
print(f"  - Total overlap: {overlap:.4f}")
print(f"  - Max difference: {max_diff:.4f}")

if overlap > 0.90 and max_diff < 0.1:
    print("  - Probability distributions match well")
else:
    print("  - Significant probability differences detected")
```

## Advanced Conversion Options

### Controlling Conversion Parameters

```python
# High-accuracy conversion
high_accuracy_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=50000,  # Very high shot count
    validate_conversion=True  # Enable detailed validation
)

# Fast conversion for testing
fast_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=1000,   # Low shot count
    validate_conversion=False  # Skip validation
)

print(f"High accuracy fidelity: {high_accuracy_result.metrics.fidelity:.4f}")
print(f"Fast conversion fidelity: {fast_result.metrics.fidelity:.4f}")
```

### Working with UCJ Results

```python
# Convert UCJ results (similar to LUCJ)
from ffsim_integration.integration import create_ucj_ansatz

ucj_result = create_ucj_ansatz(h2, n_reps=1, max_iterations=30)

ucj_conversion = ffsim_to_quri_state(
    ffsim_state_vector=ucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=5000
)

print(f"UCJ Conversion:")
print(f"  - UCJ energy: {ucj_result.final_energy:.6f} Ha")
print(f"  - Conversion fidelity: {ucj_conversion.metrics.fidelity:.4f}")
```

## Direct Amplitude Mapping

For advanced users who need direct access to computational basis amplitudes:

```python
from ffsim_integration.state_conversion import _map_fermionic_amplitudes_directly

# Direct amplitude mapping
computational_amplitudes = _map_fermionic_amplitudes_directly(
    fermionic_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec
)

print(f"Direct Mapping Results:")
print(f"  - Input dimension: {len(lucj_result.state_vector)}")
print(f"  - Output dimension: {len(computational_amplitudes)}")
print(f"  - Norm preserved: {abs(np.linalg.norm(computational_amplitudes) - 1.0) < 1e-10}")

# Analyze the computational basis amplitudes
nonzero_indices = np.where(np.abs(computational_amplitudes) > 1e-10)[0]
print(f"  - Non-zero computational states: {len(nonzero_indices)}")

# Show the most significant computational basis states
sorted_indices = np.argsort(np.abs(computational_amplitudes))[::-1]
print("\\nMost significant computational basis states:")
for i in range(min(5, len(sorted_indices))):
    idx = sorted_indices[i]
    if np.abs(computational_amplitudes[idx]) > 1e-10:
        bitstring = format(idx, f'0{n_qubits}b')
        print(f"  |{bitstring}⟩: {computational_amplitudes[idx]:.6f}")
```

## Jordan-Wigner Mapping Details

### Understanding Qubit Ordering

The conversion uses interleaved Jordan-Wigner mapping:

```python
# For H₂ with 2 orbitals (4 qubits):
# Orbital 0: qubits 0 (α) and 1 (β)  
# Orbital 1: qubits 2 (α) and 3 (β)
# Pattern: α₀β₀α₁β₁

print("Jordan-Wigner Mapping for H₂:")
print("Orbital    α qubit    β qubit")
for orbital in range(h2.ffsim_mol_data.norb):
    alpha_qubit = 2 * orbital
    beta_qubit = 2 * orbital + 1
    print(f"{orbital:7d}    {alpha_qubit:7d}    {beta_qubit:7d}")
```

### Electron Configuration Analysis

```python
# Analyze which computational basis states correspond to the correct electron count
def analyze_electron_configurations(computational_amplitudes, n_qubits, target_electrons):
    print(f"Electron Configuration Analysis:")
    print(f"Target electrons: {target_electrons}")
    print("\\nComputational states with correct electron count:")
    
    significant_states = 0
    for idx in range(len(computational_amplitudes)):
        if np.abs(computational_amplitudes[idx]) > 1e-10:
            bitstring = format(idx, f'0{n_qubits}b')
            electron_count = bitstring.count('1')
            
            if electron_count == sum(target_electrons):
                print(f"  |{bitstring}⟩: {computational_amplitudes[idx]:.6f}")
                significant_states += 1
    
    print(f"\\nTotal significant states with correct electron count: {significant_states}")

# Run the analysis
target_nelec = sum(h2.ffsim_mol_data.nelec)
analyze_electron_configurations(computational_amplitudes, n_qubits, h2.ffsim_mol_data.nelec)
```

## Troubleshooting Conversion Issues

### Low Fidelity

```python
def improve_conversion_fidelity(ffsim_state, n_qubits, nelec):
    """Try different approaches to improve conversion fidelity."""
    
    methods_to_try = [
        {"method": "sampling_circuit", "shots": 1000},
        {"method": "sampling_circuit", "shots": 10000},
        {"method": "sampling_circuit", "shots": 50000},
        {"method": "direct_mapping", "shots": None}
    ]
    
    best_fidelity = 0.0
    best_result = None
    
    print("Trying different conversion approaches:")
    print("Method            Shots     Fidelity")
    
    for config in methods_to_try:
        try:
            if config["shots"] is not None:
                result = ffsim_to_quri_state(
                    ffsim_state, n_qubits, nelec,
                    conversion_method=config["method"],
                    total_shots=config["shots"]
                )
            else:
                result = ffsim_to_quri_state(
                    ffsim_state, n_qubits, nelec,
                    conversion_method=config["method"]
                )
            
            fidelity = result.metrics.fidelity
            shots_str = str(config["shots"]) if config["shots"] else "N/A"
            print(f"{config['method']:16s} {shots_str:>8s} {fidelity:>8.4f}")
            
            if fidelity > best_fidelity:
                best_fidelity = fidelity
                best_result = result
                
        except Exception as e:
            print(f"{config['method']:16s} {'Failed':>8s} {str(e)[:20]:>8s}")
    
    return best_result

# Use the troubleshooting function
best_conversion = improve_conversion_fidelity(
    lucj_result.state_vector, 
    n_qubits, 
    h2.ffsim_mol_data.nelec
)
```

### Memory Issues

For large systems, use chunked processing:

```python
# For large molecules, consider state vector compression
def efficient_conversion_for_large_systems(ffsim_state, n_qubits, nelec):
    """Convert large states efficiently."""
    
    # Use fewer shots for initial testing
    if n_qubits > 10:
        total_shots = 1000
        print(f"Large system detected ({n_qubits} qubits), using reduced shots")
    else:
        total_shots = 5000
    
    return ffsim_to_quri_state(
        ffsim_state, n_qubits, nelec,
        conversion_method="sampling_circuit",
        total_shots=total_shots,
        validate_conversion=(n_qubits <= 8)  # Skip validation for very large systems
    )
```

## Integration with QSCI

### Using Converted States with QSCI

```python
from qsci_algorithms import VanillaQSCI

# Convert state and use with QSCI
conversion_result = ffsim_to_quri_state(
    lucj_result.state_vector, n_qubits, nelec,
    conversion_method="sampling_circuit", total_shots=5000
)

# Create QSCI instance
qsci = VanillaQSCI(
    hamiltonian=h2.quri_hamiltonian,
    num_states_pick_out=8
)

# Run QSCI with converted state
qsci_result = qsci.run(conversion_result.quri_state)

print(f"Integrated Results:")
print(f"  - Original LUCJ energy: {lucj_result.final_energy:.6f} Ha")
print(f"  - QSCI refined energy: {qsci_result.ground_state_energy:.6f} Ha")
print(f"  - Conversion fidelity: {conversion_result.metrics.fidelity:.4f}")
```

## Best Practices

1. **Start with sampling_circuit method** - Most robust for general use
2. **Use adequate shot counts** - At least 5000 for production, 1000+ for testing
3. **Check fidelity** - Should be > 0.85 for reliable results
4. **Validate electron count** - Ensure computational states have correct number of electrons
5. **Monitor memory usage** - Large systems may need reduced shot counts

## Next Steps

- **[Getting Started](getting_started.md)**: See complete conversion examples in context
- **[Ansatz Creation](ansatz_creation.md)**: Create states to convert
- **[API Reference](../../api/ffsim_integration.md)**: Detailed function documentation