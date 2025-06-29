# Molecular Systems Guide

The `molecular_systems` module provides utilities for creating and configuring molecular systems for quantum chemistry calculations. This guide covers the available functions and their parameters.

## Overview

The molecular systems module supports:
- H₂ (hydrogen) molecule creation
- N₂ (nitrogen) molecule creation with active space configuration
- Reference energy lookup and validation
- Automatic Hartree-Fock and FCI calculations

## H₂ Molecule Creation

### Basic Usage

```python
from ffsim_integration.molecular_systems import create_h2_molecule

# Create H₂ with default parameters
h2 = create_h2_molecule()
print(f"HF energy: {h2.hartree_fock_energy:.6f} Ha")
```

### Customizing Parameters

```python
# Create H₂ with specific parameters
h2_custom = create_h2_molecule(
    basis="6-31g",           # Basis set
    bond_length=0.8,         # Bond length in Angstroms
    charge=0,                # Molecular charge
    spin=0                   # Spin multiplicity (2S)
)

print(f"Molecule: {h2_custom.name}")
print(f"Basis: {h2_custom.basis}")
print(f"Bond length: {h2_custom.bond_length} Å")
print(f"Active space: {h2_custom.active_space.n_active_ele} electrons, {h2_custom.active_space.n_active_orb} orbitals")
```

### Available Basis Sets

Common basis sets for H₂:
- `"sto-3g"`: Minimal basis, fast calculations
- `"6-31g"`: Split-valence basis, better accuracy
- `"6-31g*"`: With polarization functions
- `"cc-pvdz"`: Correlation-consistent basis

## N₂ Molecule Creation

### Basic N₂ System

```python
from ffsim_integration.molecular_systems import create_n2_molecule

# Create N₂ with active space configuration
n2 = create_n2_molecule(
    basis="sto-3g",
    bond_length=1.09,        # Equilibrium bond length
    active_space=(6, 6)      # (n_electrons, n_orbitals) in active space
)

print(f"N₂ system created:")
print(f"  - HF energy: {n2.hartree_fock_energy:.6f} Ha")
print(f"  - Active space: {n2.active_space.n_active_ele} electrons, {n2.active_space.n_active_orb} orbitals")
print(f"  - Total qubits needed: {2 * n2.ffsim_mol_data.norb}")
```

### Active Space Configuration

The active space determines which orbitals and electrons are included in the calculation:

```python
# Small active space (faster, less accurate)
n2_small = create_n2_molecule(
    basis="sto-3g",
    bond_length=1.09,
    active_space=(4, 4)      # 4 electrons in 4 orbitals
)

# Larger active space (slower, more accurate)
n2_large = create_n2_molecule(
    basis="sto-3g", 
    bond_length=1.09,
    active_space=(8, 8)      # 8 electrons in 8 orbitals
)

print(f"Small active space qubits: {2 * n2_small.ffsim_mol_data.norb}")
print(f"Large active space qubits: {2 * n2_large.ffsim_mol_data.norb}")
```

## Reference Energies

### Looking Up Reference Values

```python
from ffsim_integration.molecular_systems import get_reference_energies

# Get reference energies for comparison
h2_refs = get_reference_energies("H2", "6-31g")
print("H₂ Reference Energies:")
print(f"  - Hartree-Fock: {h2_refs['hartree_fock']:.6f} Ha")
print(f"  - FCI: {h2_refs['fci']:.6f} Ha")
print(f"  - Experimental: {h2_refs['experimental']:.6f} Ha")

n2_refs = get_reference_energies("N2", "sto-3g")
print("\\nN₂ Reference Energies:")
print(f"  - Hartree-Fock: {n2_refs['hartree_fock']:.6f} Ha")
print(f"  - FCI: {n2_refs['fci']:.6f} Ha")
```

## System Validation

### Checking System Consistency

```python
from ffsim_integration.molecular_systems import validate_molecular_system

# Validate a molecular system
is_valid = validate_molecular_system(h2)
print(f"H₂ system valid: {is_valid}")

# The validation checks:
# - Energy calculations completed successfully
# - Active space configuration is reasonable
# - ffsim and QURI Parts data are consistent
# - Basis set is supported
```

## Understanding the MolecularSystem Object

Each created molecular system contains several important attributes:

```python
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

# Basic properties
print(f"Name: {h2.name}")
print(f"Basis: {h2.basis}")
print(f"Bond length: {h2.bond_length}")
print(f"Charge: {h2.charge}")
print(f"Spin: {h2.spin}")

# Energy information
print(f"\\nEnergies:")
print(f"  - Hartree-Fock: {h2.hartree_fock_energy:.6f} Ha")
print(f"  - FCI: {h2.fci_energy:.6f} Ha")
print(f"  - Correlation: {h2.fci_energy - h2.hartree_fock_energy:.6f} Ha")

# Active space details
print(f"\\nActive Space:")
print(f"  - Electrons: {h2.active_space.n_active_ele}")
print(f"  - Orbitals: {h2.active_space.n_active_orb}")

# ffsim-specific data
print(f"\\nffsim Data:")
print(f"  - Norb: {h2.ffsim_mol_data.norb}")
print(f"  - Nelec: {h2.ffsim_mol_data.nelec}")

# QURI Parts Hamiltonian
print(f"\\nQURI Parts:")
print(f"  - Hamiltonian terms: {len(h2.quri_hamiltonian)}")
```

## Performance Considerations

### Choosing Parameters

**For Quick Testing:**
```python
# Fast configuration for development/testing
h2_fast = create_h2_molecule(basis="sto-3g", bond_length=0.74)
n2_fast = create_n2_molecule(basis="sto-3g", active_space=(4, 4))
```

**For Accurate Results:**
```python
# More accurate configuration for production
h2_accurate = create_h2_molecule(basis="6-31g", bond_length=0.74)
n2_accurate = create_n2_molecule(basis="6-31g", active_space=(6, 6))
```

### Memory and Time Scaling

- **Basis set size**: Larger basis sets require more memory and time
- **Active space size**: Scales exponentially with number of orbitals
- **Bond length**: Affects convergence; use literature values when possible

## Examples for Different Applications

### Potential Energy Curves

```python
import numpy as np

# Generate H₂ potential energy curve
bond_lengths = np.linspace(0.5, 2.0, 10)
energies = []

for r in bond_lengths:
    h2 = create_h2_molecule(basis="sto-3g", bond_length=r)
    energies.append(h2.hartree_fock_energy)

# Plot or analyze the curve
print("Bond lengths and HF energies:")
for r, e in zip(bond_lengths, energies):
    print(f"  {r:.2f} Å: {e:.6f} Ha")
```

### Basis Set Comparison

```python
# Compare different basis sets for H₂
basis_sets = ["sto-3g", "6-31g", "6-31g*"]
bond_length = 0.74

print("Basis set comparison for H₂:")
for basis in basis_sets:
    h2 = create_h2_molecule(basis=basis, bond_length=bond_length)
    print(f"  {basis:8s}: {h2.hartree_fock_energy:.6f} Ha")
```

## Next Steps

- **[Ansatz Creation](ansatz_creation.md)**: Use these molecular systems to create LUCJ/UCJ ansatz states
- **[Getting Started](getting_started.md)**: See complete workflow examples
- **[API Reference](../../api/ffsim_integration.md)**: Detailed function documentation