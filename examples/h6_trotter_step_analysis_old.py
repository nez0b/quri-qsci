
"""
H6 molecule TE-QSCI Trotter step size analysis.

This script investigates the effect of the Trotter step size on the accuracy of
TE-QSCI calculations for the H6 molecule. It plots the energy error
(E_QSCI - E_exact) against the evolution time for various Trotter step sizes.
"""

import numpy as np
import matplotlib.pyplot as plt


# QURI Parts imports
from quri_parts.core.operator import Operator
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Quantum chemistry imports
from pyscf import gto, scf, fci
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner

# Our QSCI implementations
from quri_qsci import (
    create_qsci_algorithm,
    QSCIVariant,
)

def create_h6_hamiltonian_and_hf_state():
    """Creates the H6 Hamiltonian, HF state, and calculates the exact FCI energy."""
    molecule_geometry = "H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3; H 0 0 4; H 0 0 5"
    geometry_list = []
    for part in molecule_geometry.split(";"):
        part = part.strip()
        if part:
            atom_coords = part.split()
            if len(atom_coords) == 4:
                atom = atom_coords[0]
                x, y, z = atom_coords[1:4]
                geometry_list.append([atom, (float(x), float(y), float(z))])

    mole = gto.M(
        atom=geometry_list,
        basis="sto-3g",
        charge=0,
        spin=0,
        verbose=0
    )
    mf = scf.RHF(mole).run()

    active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
    hamiltonian, _ = get_qubit_mapped_hamiltonian(
        active_space,
        mo_eint_set,
        fermion_qubit_mapping=jordan_wigner
    )

    # Calculate exact FCI energy
    cisolver = fci.FCI(mole, mf.mo_coeff)
    exact_fci_energy, _ = cisolver.kernel()

    # Create Hartree-Fock state
    n_qubits = 2 * active_space.n_active_orb
    n_electrons = active_space.n_active_ele
    hf_circuit = QuantumCircuit(n_qubits)
    for i in range(n_electrons):
        hf_circuit.add_X_gate(i)
    hf_state = GeneralCircuitQuantumState(n_qubits, hf_circuit)

    return hamiltonian, hf_state, exact_fci_energy

def run_trotter_analysis():
    """
    Performs the TE-QSCI calculation for different Trotter step sizes and plots the results.
    """
    print("Setting up H6 molecule...")
    hamiltonian, hf_state, exact_energy = create_h6_hamiltonian_and_hf_state()
    print(f"H6 Hamiltonian created. Exact FCI energy: {exact_energy:.6f} Ha")

    # Parameters for the study
    evolution_times = np.linspace(0, 3.0, 10)
    trotter_step_sizes = [0.01, 0.1, 0.2, 0.5]
    subspace_size_R = 90
    total_shots = 10000

    sampler = create_qulacs_vector_concurrent_sampler()
    results = {dt: [] for dt in trotter_step_sizes}

    print("\nStarting TE-QSCI calculations for different Trotter step sizes...")
    for dt in trotter_step_sizes:
        print(f"  Trotter step size (dt): {dt}")
        energy_errors = []
        for t in evolution_times:
            if t > 0:
                trotter_steps = max(1, int(round(t / dt)))
            else:
                trotter_steps = 1  # For t=0, any number of steps is fine

            print(f"    t = {t:.2f}, trotter_steps = {trotter_steps}")

            algorithm = create_qsci_algorithm(
                QSCIVariant.SINGLE_TIME_TE,
                hamiltonian=hamiltonian,
                sampler=sampler,
                num_states_pick_out=subspace_size_R,
                evolution_time=t,
                trotter_steps=trotter_steps,
                time_evolution_method="trotter"
            )

            qsci_result = algorithm.run(hf_state, total_shots)
            energy_error = qsci_result.ground_state_energy - exact_energy
            energy_errors.append(energy_error)

        results[dt] = energy_errors

    # Plotting the results
    print("\nPlotting results...")
    plt.figure(figsize=(10, 7))
    for dt in trotter_step_sizes:
        plt.scatter(evolution_times, results[dt], label=f"dt = {dt}")

    plt.xlabel("Evolution Time (t)")
    plt.ylabel("Energy Error (E_QSCI - E_exact) [Ha]")
    plt.title(f"TE-QSCI Trotter Error Analysis for H6/STO-3G (R={subspace_size_R})")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    plot_filename = "h6_trotter_analysis.png"
    plt.savefig(plot_filename)
    print(f"Plot saved to {plot_filename}")
    plt.show()

if __name__ == "__main__":
    run_trotter_analysis()
