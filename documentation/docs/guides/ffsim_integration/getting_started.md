# Getting Started with ffsim Integration

This tutorial provides a complete, step-by-step guide to using the `ffsim_integration` module with the QSCI framework. By the end of this tutorial, you'll be able to create molecular systems, generate ansatz states, and run QSCI calculations.

## Prerequisites

Ensure you have installed the optional ffsim dependency:

```bash
pip install quri-qsci[ffsim]
```

## Complete Example: H₂ Molecule with LUCJ + QSCI

Let's walk through a complete workflow using the H₂ molecule:

### Step 1: Create a Molecular System

```python
from ffsim_integration.molecular_systems import create_h2_molecule

# Create H₂ molecule with STO-3G basis at equilibrium bond length
h2_system = create_h2_molecule(
    basis="sto-3g",
    bond_length=0.74,  # Angstroms
    charge=0,
    spin=0
)

print(f"Created H₂ system:")
print(f"  - Basis: {h2_system.basis}")
print(f"  - Bond length: {h2_system.bond_length} Å")
print(f"  - HF energy: {h2_system.hartree_fock_energy:.6f} Ha")
print(f"  - Active space: {h2_system.active_space.n_active_ele} electrons, {h2_system.active_space.n_active_orb} orbitals")
```

### Step 2: Generate LUCJ Ansatz

```python
from ffsim_integration.integration import create_lucj_ansatz

# Create LUCJ (Linear Unitary Coupled Cluster Jastrow) ansatz
lucj_result = create_lucj_ansatz(
    mol_system=h2_system,
    n_reps=1,                      # Number of repetitions
    optimization_method="BFGS",    # Optimization algorithm
    max_iterations=50              # Maximum optimization steps
)

print(f"\\nLUCJ Ansatz Results:")
print(f"  - Final energy: {lucj_result.final_energy:.6f} Ha")
print(f"  - Energy improvement: {h2_system.hartree_fock_energy - lucj_result.final_energy:.6f} Ha")
print(f"  - Optimization success: {lucj_result.optimization_success}")
print(f"  - State vector dimension: {len(lucj_result.state_vector)}")
```

### Step 3: Convert State for QSCI

```python
from ffsim_integration.state_conversion import ffsim_to_quri_state

# Convert ffsim state vector to QURI Parts format
n_qubits = 2 * h2_system.ffsim_mol_data.norb
nelec = h2_system.ffsim_mol_data.nelec

quri_state_result = ffsim_to_quri_state(
    ffsim_state_vector=lucj_result.state_vector,
    n_qubits=n_qubits,
    nelec=nelec,
    conversion_method="sampling_circuit",
    total_shots=5000
)

print(f"\\nState Conversion Results:")
print(f"  - Conversion fidelity: {quri_state_result.metrics.fidelity:.4f}")
print(f"  - State vector norm: {quri_state_result.metrics.state_vector_norm:.6f}")
print(f"  - QURI Parts state: {type(quri_state_result.quri_state)}")
```

### Step 4: Run QSCI Analysis

```python
from qsci_algorithms import VanillaQSCI

# Create QSCI algorithm instance
qsci = VanillaQSCI(
    hamiltonian=h2_system.quri_hamiltonian,
    num_states_pick_out=8  # Number of states to select
)

# Run QSCI with the converted state
qsci_result = qsci.run(quri_state_result.quri_state)

print(f"\\nQSCI Results:")
print(f"  - Ground state energy: {qsci_result.ground_state_energy:.6f} Ha")
print(f"  - Selected states: {len(qsci_result.selected_states)}")
print(f"  - Improvement over HF: {h2_system.hartree_fock_energy - qsci_result.ground_state_energy:.6f} Ha")
```

### Step 5: Complete Workflow with High-Level Interface

For convenience, you can run the entire workflow with a single function:

```python
from ffsim_integration.qsci_interface import run_lucj_qsci

# Run complete LUCJ + QSCI workflow
workflow_result = run_lucj_qsci(
    molecule="H2",
    ansatz_type="LUCJ",
    subspace_sizes=[4, 8, 12],     # Different QSCI subspace sizes to try
    basis="sto-3g",
    bond_length=0.74,
    n_reps=1,
    max_optimization_iterations=50,
    total_shots=5000,
    conversion_method="sampling_circuit"
)

print(f"\\nComplete Workflow Results:")
print(f"  - Best QSCI energy: {workflow_result.best_energy:.6f} Ha")
print(f"  - Best subspace size: R = {workflow_result.best_subspace_size}")
print(f"  - LUCJ ansatz energy: {workflow_result.ansatz_energy:.6f} Ha")
print(f"  - Conversion fidelity: {workflow_result.conversion_metrics.fidelity:.4f}")
```

## Expected Output

When you run this complete example, you should see output similar to:

```
Created H₂ system:
  - Basis: sto-3g
  - Bond length: 0.74 Å
  - HF energy: -1.116759 Ha
  - Active space: 2 electrons, 2 orbitals

LUCJ Ansatz Results:
  - Final energy: -1.117763 Ha
  - Energy improvement: 0.001004 Ha
  - Optimization success: True
  - State vector dimension: 4

State Conversion Results:
  - Conversion fidelity: 0.8500
  - State vector norm: 1.000000
  - QURI Parts state: <class 'quri_parts.core.state.ComputationalBasisSuperposition'>

QSCI Results:
  - Ground state energy: -1.117500 Ha
  - Selected states: 4
  - Improvement over HF: 0.000741 Ha

Complete Workflow Results:
  - Best QSCI energy: -1.117500 Ha
  - Best subspace size: R = 8
  - LUCJ ansatz energy: -1.117763 Ha
  - Conversion fidelity: 0.8500
```

## Understanding the Results

- **Energy Improvement**: The LUCJ ansatz provides a better energy than Hartree-Fock
- **Conversion Fidelity**: Measures how accurately the state was converted (closer to 1.0 is better)
- **QSCI Enhancement**: Further refinement of the energy using selected configuration interaction

## Next Steps

Now that you've completed the basic workflow, explore:

- **[Molecular Systems](molecular_systems.md)**: Learn about N₂ molecules and active space configuration
- **[Ansatz Creation](ansatz_creation.md)**: Understand the differences between LUCJ and UCJ methods
- **[State Conversion](state_conversion.md)**: Deep dive into state format conversion details
- **[API Reference](../../api/ffsim_integration.md)**: Complete function documentation

## Troubleshooting

**Import Errors**: Ensure ffsim is installed with `pip install quri-qsci[ffsim]`

**Optimization Failures**: Try reducing `max_iterations` or changing the `optimization_method` to "Powell" or "Nelder-Mead"

**Low Fidelity**: Increase `total_shots` for better state conversion accuracy

**Performance Issues**: For larger molecules, consider using smaller active spaces or fewer optimization iterations