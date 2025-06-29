# Ansatz Creation Guide

The `integration` module provides functions for creating and optimizing variational quantum ansatz states using the ffsim library. This guide covers both LUCJ (Linear Unitary Coupled Cluster Jastrow) and UCJ (Unitary Coupled Cluster Jastrow) methods.

## Overview

Variational quantum ansatz are parameterized quantum circuits that can approximate ground states of molecular systems. The ffsim integration provides:

- **LUCJ**: Linear UCJ with reduced parameter count and faster optimization
- **UCJ**: Full unitary coupled cluster approach with comprehensive correlation treatment
- **Automatic optimization**: Built-in parameter optimization using classical methods
- **State vector extraction**: Direct access to optimized quantum states

## LUCJ (Linear Unitary Coupled Cluster Jastrow)

### Basic LUCJ Creation

```python
from ffsim_integration.molecular_systems import create_h2_molecule
from ffsim_integration.integration import create_lucj_ansatz

# Create molecular system
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

# Generate LUCJ ansatz
lucj_result = create_lucj_ansatz(
    mol_system=h2,
    n_reps=1,                    # Number of repetitions
    optimization_method="BFGS",  # Optimization algorithm
    max_iterations=50            # Maximum optimization steps
)

print(f"LUCJ Results:")
print(f"  - Final energy: {lucj_result.final_energy:.6f} Ha")
print(f"  - Optimization success: {lucj_result.optimization_success}")
print(f"  - Parameters: {len(lucj_result.optimized_parameters)}")
print(f"  - Iterations used: {lucj_result.n_iterations}")
```

### Understanding LUCJ Parameters

LUCJ uses fewer parameters than full UCJ, making optimization faster:

```python
# Compare parameter counts for different repetitions
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

for n_reps in [1, 2, 3]:
    lucj_result = create_lucj_ansatz(h2, n_reps=n_reps, max_iterations=10)
    print(f"n_reps={n_reps}: {len(lucj_result.optimized_parameters)} parameters")
```

### LUCJ Optimization Methods

Different optimization algorithms have different strengths:

```python
# Try different optimization methods
optimization_methods = ["BFGS", "Powell", "Nelder-Mead", "CG"]

for method in optimization_methods:
    try:
        lucj_result = create_lucj_ansatz(
            h2, 
            n_reps=1, 
            optimization_method=method,
            max_iterations=30
        )
        print(f"{method:12s}: {lucj_result.final_energy:.6f} Ha ({lucj_result.n_iterations} iter)")
    except Exception as e:
        print(f"{method:12s}: Failed - {e}")
```

## UCJ (Unitary Coupled Cluster Jastrow)

### Basic UCJ Creation

```python
from ffsim_integration.integration import create_ucj_ansatz

# Generate UCJ ansatz (more parameters, potentially more accurate)
ucj_result = create_ucj_ansatz(
    mol_system=h2,
    n_reps=1,
    optimization_method="BFGS",
    max_iterations=50
)

print(f"UCJ Results:")
print(f"  - Final energy: {ucj_result.final_energy:.6f} Ha")
print(f"  - Parameters: {len(ucj_result.optimized_parameters)}")
print(f"  - State vector norm: {np.linalg.norm(ucj_result.state_vector):.6f}")
```

### LUCJ vs UCJ Comparison

```python
import numpy as np

# Compare LUCJ and UCJ for the same system
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

lucj_result = create_lucj_ansatz(h2, n_reps=1, max_iterations=50)
ucj_result = create_ucj_ansatz(h2, n_reps=1, max_iterations=50)

print("LUCJ vs UCJ Comparison:")
print(f"Method    Energy (Ha)    Parameters    Iterations")
print(f"LUCJ      {lucj_result.final_energy:.6f}    {len(lucj_result.optimized_parameters):10d}    {lucj_result.n_iterations:10d}")
print(f"UCJ       {ucj_result.final_energy:.6f}    {len(ucj_result.optimized_parameters):10d}    {ucj_result.n_iterations:10d}")
print(f"HF        {h2.hartree_fock_energy:.6f}    {'N/A':>10s}    {'N/A':>10s}")
```

## Advanced Configuration

### Repetition Layers

The `n_reps` parameter controls the depth of the ansatz circuit:

```python
# Study the effect of repetition layers
bond_length = 0.74
n2 = create_n2_molecule(basis="sto-3g", bond_length=1.09, active_space=(6, 6))

print("Effect of repetition layers on N₂:")
print("n_reps    Energy (Ha)    Parameters    Time (approx)")

for n_reps in [1, 2, 3]:
    import time
    start_time = time.time()
    
    lucj_result = create_lucj_ansatz(
        n2, 
        n_reps=n_reps, 
        max_iterations=20  # Reduced for timing
    )
    
    elapsed = time.time() - start_time
    print(f"{n_reps:6d}    {lucj_result.final_energy:.6f}    {len(lucj_result.optimized_parameters):10d}    {elapsed:.1f}s")
```

### Optimization Convergence

Monitor optimization progress:

```python
# Create a custom optimization with detailed monitoring
def detailed_lucj_optimization(mol_system, n_reps=1):
    from ffsim_integration.integration import create_lucj_ansatz
    
    # Try different iteration limits to see convergence
    iteration_limits = [10, 25, 50, 100]
    
    print("Optimization Convergence Study:")
    print("Max Iter    Final Energy    Success    Parameters")
    
    for max_iter in iteration_limits:
        result = create_lucj_ansatz(
            mol_system,
            n_reps=n_reps,
            optimization_method="BFGS",
            max_iterations=max_iter
        )
        
        print(f"{max_iter:8d}    {result.final_energy:.6f}    {str(result.optimization_success):7s}    {len(result.optimized_parameters):10d}")

# Run the convergence study
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)
detailed_lucj_optimization(h2, n_reps=1)
```

## Working with Ansatz Results

### Extracting State Information

```python
# Analyze the resulting quantum state
lucj_result = create_lucj_ansatz(h2, n_reps=1, max_iterations=50)

print("State Vector Analysis:")
print(f"  - Dimension: {len(lucj_result.state_vector)}")
print(f"  - Norm: {np.linalg.norm(lucj_result.state_vector):.6f}")
print(f"  - Max amplitude: {np.max(np.abs(lucj_result.state_vector)):.6f}")
print(f"  - Non-zero elements: {np.sum(np.abs(lucj_result.state_vector) > 1e-10)}")

# Look at the largest amplitude components
amplitudes = np.abs(lucj_result.state_vector)
sorted_indices = np.argsort(amplitudes)[::-1]

print("\\nLargest amplitude components:")
for i in range(min(5, len(sorted_indices))):
    idx = sorted_indices[i]
    print(f"  State {idx:2d}: amplitude = {amplitudes[idx]:.6f}")
```

### Energy Analysis

```python
# Compare energies across different methods
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

# Reference energies
hf_energy = h2.hartree_fock_energy
fci_energy = h2.fci_energy

# Ansatz energies
lucj_result = create_lucj_ansatz(h2, n_reps=1, max_iterations=50)
lucj_energy = lucj_result.final_energy

print("Energy Comparison:")
print(f"Method           Energy (Ha)    vs HF (Ha)    % Recovery")
print(f"Hartree-Fock     {hf_energy:.6f}    {0.0:.6f}    {0.0:8.1f}%")
print(f"LUCJ             {lucj_energy:.6f}    {lucj_energy - hf_energy:.6f}    {100 * (lucj_energy - hf_energy) / (fci_energy - hf_energy):8.1f}%")
print(f"FCI              {fci_energy:.6f}    {fci_energy - hf_energy:.6f}    {100.0:8.1f}%")
```

## Best Practices

### Parameter Selection Guidelines

**For Development/Testing:**
```python
# Fast configuration for testing
test_result = create_lucj_ansatz(
    mol_system=h2,
    n_reps=1,                    # Start with 1 repetition
    optimization_method="Powell", # Often faster than BFGS
    max_iterations=20            # Limit iterations for speed
)
```

**For Production/Accuracy:**
```python
# Accurate configuration for final results
production_result = create_lucj_ansatz(
    mol_system=h2,
    n_reps=2,                    # More repetitions for accuracy
    optimization_method="BFGS",  # Usually most reliable
    max_iterations=100           # Allow full convergence
)
```

### Error Handling

```python
# Robust ansatz creation with error handling
def safe_lucj_creation(mol_system, **kwargs):
    try:
        result = create_lucj_ansatz(mol_system, **kwargs)
        
        if not result.optimization_success:
            print("Warning: Optimization did not converge")
            print(f"Final energy: {result.final_energy:.6f} Ha")
            print(f"Iterations used: {result.n_iterations}")
        
        return result
        
    except Exception as e:
        print(f"Error creating LUCJ ansatz: {e}")
        print("Try reducing max_iterations or changing optimization_method")
        return None

# Use the safe function
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)
result = safe_lucj_creation(h2, n_reps=1, max_iterations=50)
```

## Performance Tips

1. **Start Small**: Begin with `n_reps=1` and small iteration limits
2. **Choose Method Wisely**: BFGS is usually best, but Powell can be faster
3. **Monitor Convergence**: Check `optimization_success` flag
4. **Scale Gradually**: Increase complexity only when needed
5. **Use Validation**: Always compare with Hartree-Fock and FCI when available

## Next Steps

- **[State Conversion](state_conversion.md)**: Convert ansatz states for use with QSCI
- **[Getting Started](getting_started.md)**: See complete workflow examples  
- **[Molecular Systems](molecular_systems.md)**: Create different molecular systems to test
- **[API Reference](../../api/ffsim_integration.md)**: Detailed function documentation