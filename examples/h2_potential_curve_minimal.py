#!/usr/bin/env python3
"""
Minimal H2 Potential Energy Curve using QSCI methods.

This is a simplified version that focuses on generating the comparison plot
with reduced computational requirements for quick demonstration.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# PySCF and quantum chemistry imports
from pyscf import scf, gto
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner

# QURI Parts core imports
from quri_parts.core.operator import Operator
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit

# QURI Parts sampler imports
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Our QSCI implementations
from src import (
    create_qsci_algorithm, QSCIVariant,
    VanillaQSCIAlgorithm, SingleTimeTeQSCIAlgorithm, TimeAverageTeQSCIAlgorithm
)


def get_h2_hamiltonian_and_hf_state(r: float):
    """Generate H2 Hamiltonian and Hartree-Fock state for given bond length."""
    # Create H2 molecule
    mole = gto.M(atom=f"H 0 0 0; H 0 0 {r}", basis="sto-3g", verbose=0)
    mf = scf.RHF(mole)
    mf.run(verbose=0)
    
    # Generate molecular orbital integrals
    active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
    
    # Create qubit Hamiltonian using Jordan-Wigner mapping
    hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
        active_space, mo_eint_set, fermion_qubit_mapping=jordan_wigner
    )
    
    # Create Hartree-Fock state
    n_qubits = 2 * active_space.n_active_orb  # Should be 4 for H2
    n_electrons = active_space.n_active_ele    # Should be 2 for H2
    
    hf_circuit = QuantumCircuit(n_qubits)
    for i in range(n_electrons):
        hf_circuit.add_X_gate(i)
    
    hf_state = GeneralCircuitQuantumState(n_qubits, hf_circuit)
    
    return hamiltonian, hf_state, mf.e_tot


def run_qsci_calculation(variant, hamiltonian, hf_state, shots=1000, num_states=16):
    """Run QSCI calculation for given variant."""
    
    sampler = create_qulacs_vector_concurrent_sampler()
    
    if variant == "vanilla":
        algorithm = create_qsci_algorithm(
            QSCIVariant.VANILLA,
            hamiltonian,
            sampler=sampler,
            num_states_pick_out=num_states
        )
        result = algorithm.run(
            input_states=[hf_state],
            total_shots=shots
        )
    elif variant == "single_time_te":
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=sampler,
            evolution_time=1.0,
            num_states_pick_out=num_states
        )
        result = algorithm.run(
            initial_state=hf_state,
            total_shots=shots
        )
    elif variant == "time_average_te":
        algorithm = create_qsci_algorithm(
            QSCIVariant.TIME_AVERAGE_TE,
            hamiltonian,
            sampler=sampler,
            evolution_times=[0.5, 1.0, 1.5],
            num_states_pick_out=num_states
        )
        shots_per_time = shots // 3
        result = algorithm.run(
            initial_state=hf_state,
            shots_per_time=shots_per_time
        )
    else:
        raise ValueError(f"Unknown variant: {variant}")
    
    return result.ground_state_energy


def main():
    """Run minimal H2 potential curve study."""
    
    print("H2 Molecule Potential Energy Curve - Minimal Demo")
    print("=" * 60)
    
    # Denser sampling for smoother curves
    bond_lengths = np.linspace(0.5, 2.0, 20)  # Increased to 20 points for smoother curves
    shots = 1000  # Reduced shots for speed
    num_states = 16  # Reduced states for speed
    
    print(f"Bond lengths: {len(bond_lengths)} points from {bond_lengths[0]:.1f} to {bond_lengths[-1]:.1f} Å")
    print(f"QSCI shots: {shots}")
    print(f"States to pick: {num_states}")
    print()
    
    # Storage for results
    results = {
        'bond_lengths': bond_lengths,
        'hf_energies': [],
        'vanilla_energies': [],
        'single_time_te_energies': [],
        'time_average_te_energies': []
    }
    
    # Run calculations
    for i, r in enumerate(bond_lengths):
        print(f"Bond length {r:.2f} Å ({i+1}/{len(bond_lengths)}):")
        
        # Get system
        hamiltonian, hf_state, hf_energy = get_h2_hamiltonian_and_hf_state(r)
        results['hf_energies'].append(hf_energy)
        print(f"  HF energy: {hf_energy:.6f} Ha")
        
        # Run QSCI variants
        for variant_name, variant_key in [
            ("Vanilla QSCI", "vanilla"),
            ("Single-time TE-QSCI", "single_time_te"),
            ("Time-avg TE-QSCI", "time_average_te")
        ]:
            try:
                start_time = time.time()
                energy = run_qsci_calculation(variant_key, hamiltonian, hf_state, shots, num_states)
                calc_time = time.time() - start_time
                
                results[f'{variant_key}_energies'].append(energy)
                print(f"  {variant_name}: {energy:.6f} Ha ({calc_time:.1f}s)")
                
            except Exception as e:
                print(f"  {variant_name}: Failed ({e})")
                results[f'{variant_key}_energies'].append(np.nan)
        
        print()
    
    # Create plot
    print("Creating potential energy curve plot...")
    
    # Set up the plot style
    sns.set_theme("talk")
    plt.rcParams["figure.figsize"] = (12, 8)
    plt.rcParams["font.size"] = 12
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Main plot
    bond_lengths = results['bond_lengths']
    
    # Plot all methods with enhanced visibility for overlapping curves
    methods_to_plot = [
        ('Hartree-Fock', results['hf_energies'], 'gray', '--', 'Reference HF', 6, 2),
        ('Vanilla QSCI', results['vanilla_energies'], 'green', '-s', f'Vanilla QSCI (R={num_states})', 7, 2),
        ('Single-time TE-QSCI', results['single_time_te_energies'], 'red', '-^', 'Single-time TE-QSCI (t=1.0)', 9, 3),  # Larger markers, thicker line
        ('Time-average TE-QSCI', results['time_average_te_energies'], 'purple', '-d', 'Time-avg TE-QSCI', 6, 2)
    ]
    
    for method_name, energies, color, style, label, markersize, linewidth in methods_to_plot:
        energies_array = np.array(energies)
        valid_mask = ~np.isnan(energies_array)
        
        if np.any(valid_mask):
            if method_name == 'Single-time TE-QSCI':
                # Special handling for single-time TE-QSCI to make it visible even when overlapping
                ax1.plot(bond_lengths[valid_mask], energies_array[valid_mask], 
                        style, color=color, label=label, markersize=markersize, linewidth=linewidth,
                        markeredgecolor='black', markeredgewidth=1.5, zorder=10)  # Black edge and higher z-order
            else:
                ax1.plot(bond_lengths[valid_mask], energies_array[valid_mask], 
                        style, color=color, label=label, markersize=markersize, linewidth=linewidth)
    
    ax1.set_title("H$_2$ Potential Energy Curves: QSCI Methods Comparison", fontsize=16, fontweight='bold')
    ax1.set_xlabel("Bond Length (Å)", fontsize=12)
    ax1.set_ylabel("Energy (Hartree)", fontsize=12)
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Find reasonable y-limits
    all_energies = []
    for energies in [results['hf_energies'], results['vanilla_energies'], 
                     results['single_time_te_energies'], results['time_average_te_energies']]:
        valid_energies = [e for e in energies if not np.isnan(e)]
        all_energies.extend(valid_energies)
    
    if all_energies:
        min_energy = min(all_energies)
        max_energy = max(all_energies)
        energy_range = max_energy - min_energy
        ax1.set_ylim(min_energy - 0.1*energy_range, max_energy + 0.1*energy_range)
    
    # Difference plot (relative to HF)
    hf_energies = np.array(results['hf_energies'])
    
    methods_for_diff = [
        ('Vanilla QSCI', results['vanilla_energies'], 'green', '-s', 7, 2),
        ('Single-time TE-QSCI', results['single_time_te_energies'], 'red', '-^', 9, 3),  # Larger markers
        ('Time-average TE-QSCI', results['time_average_te_energies'], 'purple', '-d', 6, 2)
    ]
    
    for method_name, energies, color, style, markersize, linewidth in methods_for_diff:
        energies_array = np.array(energies)
        diff = (energies_array - hf_energies) * 1000  # Convert to milliHartree
        valid_mask = ~np.isnan(diff)
        
        if np.any(valid_mask):
            if method_name == 'Single-time TE-QSCI':
                # Special handling for single-time TE-QSCI to make it visible even when overlapping
                ax2.plot(bond_lengths[valid_mask], diff[valid_mask], 
                        style, color=color, label=f'{method_name} - HF', markersize=markersize, linewidth=linewidth,
                        markeredgecolor='black', markeredgewidth=1.5, zorder=10)  # Black edge and higher z-order
            else:
                ax2.plot(bond_lengths[valid_mask], diff[valid_mask], 
                        style, color=color, label=f'{method_name} - HF', markersize=markersize, linewidth=linewidth)
    
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.7, label='HF reference')
    ax2.set_title("Energy Differences Relative to Hartree-Fock", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Bond Length (Å)", fontsize=12)
    ax2.set_ylabel("Energy Difference (mHa)", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    # Create figures directory if it doesn't exist
    figures_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    filename = os.path.join(figures_dir, 'h2_potential_curve_comparison.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved to: {filename}")
    
    # Also save to docs/images for documentation
    docs_images_dir = os.path.join(os.path.dirname(__file__), '..', 'documentation', 'docs', 'images')
    os.makedirs(docs_images_dir, exist_ok=True)
    docs_filename = os.path.join(docs_images_dir, 'h2_potential_curve_comparison.png')
    plt.savefig(docs_filename, dpi=300, bbox_inches='tight')
    print(f"✓ Plot also saved to docs: {docs_filename}")
    
    plt.show()
    
    # Save results
    results_filename = os.path.join(figures_dir, 'h2_potential_curve_results.npz')
    np.savez(results_filename, **results)
    print(f"✓ Results saved to: {results_filename}")
    
    # Print summary
    print("\nResults Summary:")
    print("=" * 50)
    
    # Find minimum energies
    methods = {
        'Hartree-Fock': results['hf_energies'],
        'Vanilla QSCI': results['vanilla_energies'],
        'Single-time TE-QSCI': results['single_time_te_energies'],
        'Time-average TE-QSCI': results['time_average_te_energies']
    }
    
    print(f"{'Method':<20} {'Min Energy (Ha)':<15} {'Eq. Bond (Å)':<15}")
    print("-" * 50)
    
    for method_name, energies in methods.items():
        energies_array = np.array(energies)
        valid_mask = ~np.isnan(energies_array)
        
        if np.any(valid_mask):
            min_energy = np.nanmin(energies_array)
            min_idx = np.nanargmin(energies_array)
            eq_bond = bond_lengths[min_idx]
            
            print(f"{method_name:<20} {min_energy:<15.6f} {eq_bond:<15.2f}")
        else:
            print(f"{method_name:<20} {'N/A':<15} {'N/A':<15}")
    
    print(f"\n✅ H2 potential curve analysis completed successfully!")
    print(f"Generated plot: {filename}")


if __name__ == "__main__":
    main()