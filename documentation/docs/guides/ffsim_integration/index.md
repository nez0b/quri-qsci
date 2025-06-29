# ffsim Integration Overview

The `ffsim_integration` module provides a powerful bridge between the [ffsim](https://github.com/qiskit-community/ffsim) quantum simulation library and the QSCI (Quantum Selected Configuration Interaction) framework.

## Why Use ffsim with QSCI?

The integration combines the strengths of both frameworks:

- **ffsim**: Efficient fermionic quantum simulation with optimized state vector operations
- **QSCI**: Selected configuration interaction methods for quantum chemistry problems
- **Unified Workflow**: Seamless conversion between different state representations

## Key Benefits

✅ **Efficient Ansatz Creation**: Generate LUCJ and UCJ ansatz states using ffsim's optimized routines  
✅ **Flexible State Conversion**: Convert between ffsim and QURI Parts state formats  
✅ **Molecular System Support**: Built-in support for H₂ and N₂ molecules  
✅ **QSCI Integration**: Direct integration with QSCI algorithms for enhanced performance  

## Module Components

The `ffsim_integration` module consists of four main components:

### 🧪 Molecular Systems (`molecular_systems.py`)
Create and configure molecular systems for quantum chemistry calculations:
- H₂ molecule creation with customizable parameters
- N₂ molecule setup with active space configuration
- Hartree-Fock and FCI reference energy calculations

### 🔬 Ansatz Creation (`integration.py`)
Generate variational quantum ansatz states:
- **LUCJ (Linear Unitary Coupled Cluster Jastrow)**: Efficient parameterized ansatz
- **UCJ (Unitary Coupled Cluster Jastrow)**: Full unitary coupled cluster approach
- Parameter optimization and energy minimization

### 🔄 State Conversion (`state_conversion.py`)
Convert between different quantum state representations:
- ffsim fermionic states ↔ QURI Parts computational basis
- Jordan-Wigner mapping with proper qubit ordering
- Fidelity and overlap calculations

### 🎯 High-Level Interface (`qsci_interface.py`)
Complete workflow orchestration:
- End-to-end LUCJ + QSCI pipelines
- Automated molecule → ansatz → QSCI workflows
- Performance benchmarking and result analysis

## Installation

The ffsim integration is available as an optional dependency:

```bash
pip install quri-qsci[ffsim]
```

## Quick Start

Here's a simple example to get you started:

```python
from ffsim_integration.molecular_systems import create_h2_molecule
from ffsim_integration.integration import create_lucj_ansatz

# Create H₂ molecule
h2 = create_h2_molecule(basis="sto-3g", bond_length=0.74)

# Generate LUCJ ansatz
lucj_result = create_lucj_ansatz(h2, n_reps=1, max_iterations=50)

print(f"LUCJ energy: {lucj_result.final_energy:.6f} Ha")
```

## Next Steps

- 📚 **[Getting Started](getting_started.md)**: Complete tutorial with working examples
- 🧪 **[Molecular Systems](molecular_systems.md)**: Detailed guide to molecule creation
- 🔬 **[Ansatz Creation](ansatz_creation.md)**: Understanding LUCJ and UCJ methods
- 🔄 **[State Conversion](state_conversion.md)**: Working with different state formats
- 📖 **[API Reference](../../api/ffsim_integration.md)**: Complete function documentation

## Support

For questions and issues related to the ffsim integration:
- Check the [API documentation](../../api/ffsim_integration.md) for detailed function references
- Review the example notebooks in the `examples/ffsim_integration/` directory
- Report issues on the [GitHub repository](https://github.com/nez0b/quri-qsci)