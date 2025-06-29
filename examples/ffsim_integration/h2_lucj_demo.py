#!/usr/bin/env python3
"""
H2 LUCJ ansatz demonstration using quri-qsci[ffsim].

This example shows how to use LUCJ (Local Unitary Coupled Cluster Jastrow) ansatz
with QSCI for the H2 molecule.

Requires: pip install quri-qsci[ffsim]
"""

try:
    import ffsim
except ImportError:
    print("Error: ffsim is required for this example.")
    print("Install with: pip install quri-qsci[ffsim]")
    exit(1)

try:
    from ffsim_integration.qsci_interface import run_lucj_qsci
except ImportError:
    print("Error: Could not import ffsim_integration module.")
    print("Make sure quri-qsci[ffsim] is properly installed.")
    exit(1)


def main():
    """Main H2 LUCJ demonstration."""
    print("H2 LUCJ Ansatz + QSCI Demonstration")
    print("=" * 50)
    
    # Run H2 LUCJ + QSCI workflow
    result = run_lucj_qsci(
        molecule="H2",
        ansatz_type="LUCJ",
        subspace_sizes=[10, 25, 50],
        basis="sto-3g",
        bond_length=0.74,
        n_reps=1,
        max_optimization_iterations=50,
        total_shots=5000,
        conversion_method="sampling_circuit"
    )
    
    print("\n🎯 KEY INSIGHTS:")
    print(f"✓ LUCJ ansatz successfully optimized H2 molecule")
    print(f"✓ Energy improvement over HF: {result.hartree_fock_energy - result.ansatz_energy:.6f} Ha")
    print(f"✓ State conversion fidelity: {result.conversion_metrics.fidelity:.4f}")
    print(f"✓ Best QSCI subspace size: R = {result.best_subspace_size}")
    print(f"✓ LUCJ uses locality constraints for better scalability")
    
    return result


if __name__ == "__main__":
    result = main()