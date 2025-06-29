#!/usr/bin/env python3
"""
H6 Table III Reproduction: Time-Average TE-QSCI Energy Analysis

This script reproduces Table III from the paper "Quantum-selected configuration 
interaction with time-evolved state" (arXiv:2412.13839) for H6 molecule using 
STO-3G basis set.

The analysis tests time-average TE-QSCI over different time intervals and 
calculates E_QSCI - E_exact energy differences for:
- t = 1.4 (single time point)
- t ∈ [1.0, 2.0] (time interval)
- t ∈ [0.5, 1.5] (time interval) 
- t ∈ [1.0, 2.5] (time interval)
- t ∈ [0.5, 2.5] (time interval)

Uses exact time evolution (not Trotter approximation) for high accuracy.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

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
from src import (
    create_qsci_algorithm,
    QSCIVariant,
    TimeEvolvedQSCI
)


class H6TableIIIAnalysis:
    """
    Analysis class for reproducing Table III from TE-QSCI paper.
    
    Tests time-average TE-QSCI with exact time evolution for different
    time intervals and calculates energy differences relative to exact FCI.
    """
    
    def __init__(self):
        """Initialize H6 analysis with Table III parameters."""
        self.molecule_geometry = "H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3; H 0 0 4; H 0 0 5"
        self.basis_set = "sto-3g"
        self.n_electrons = 6
        self.expected_qubits = 12
        
        # Time intervals from Table III
        self.time_intervals = [
            (1.4, 1.4),      # Single time point: t = 1.4
            (1.0, 2.0),      # Time interval: t ∈ [1.0, 2.0]
            (0.5, 1.5),      # Time interval: t ∈ [0.5, 1.5]
            (1.0, 2.5),      # Time interval: t ∈ [1.0, 2.5]
            (0.5, 2.5)       # Time interval: t ∈ [0.5, 2.5]
        ]
        
        # QSCI parameters (estimated based on existing H6 studies)
        self.subspace_size_R = 50   # Reduced for faster execution
        self.shots_per_time = 1000  # Reduced for faster execution
        self.num_time_points = 11   # Reduced for faster execution
        
        print("H6 TABLE III REPRODUCTION - TIME-AVERAGE TE-QSCI ANALYSIS")
        print("=" * 80)
        print(f"Target: Reproduce Table III energy differences for H6/STO-3G")
        print(f"Method: Exact time evolution with time-average TE-QSCI")
        print()
        
        # Initialize molecular system and references
        self.hamiltonian = None
        self.hf_state = None
        self.exact_fci_energy = None
        self.sampler = None
        
    def setup_h6_system(self):
        """
        Set up H6 molecular system and calculate exact reference energy.
        
        Returns:
            Tuple of (hamiltonian, hf_state, exact_fci_energy)
        """
        print("✓ Setting up H6 molecular system...")
        
        # Parse geometry string into list format
        geometry_list = []
        for part in self.molecule_geometry.split(";"):
            part = part.strip()
            if part:
                atom_coords = part.split()
                if len(atom_coords) == 4:
                    atom = atom_coords[0]
                    x, y, z = atom_coords[1:4]
                    geometry_list.append([atom, (float(x), float(y), float(z))])
        
        # Create PySCF molecule object
        mole = gto.M(
            atom=geometry_list,
            basis=self.basis_set,
            charge=0,
            spin=0,
            verbose=0
        )
        
        # Run Hartree-Fock calculation
        mf = scf.RHF(mole).run()
        
        print(f"  • System: Linear H6 chain with {self.basis_set} basis")
        print(f"  • {self.n_electrons} electrons, {self.expected_qubits} qubits expected")
        print(f"  • Hartree-Fock energy: {mf.e_tot:.6f} Ha")
        
        # Generate molecular orbital integrals and qubit Hamiltonian
        active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
        hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
            active_space,
            mo_eint_set,
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Verify system dimensions
        n_qubits_actual = 2 * active_space.n_active_orb
        n_electrons_actual = active_space.n_active_ele
        
        assert n_qubits_actual == self.expected_qubits, f"Expected {self.expected_qubits} qubits, got {n_qubits_actual}"
        assert n_electrons_actual == self.n_electrons, f"Expected {self.n_electrons} electrons, got {n_electrons_actual}"
        
        # Calculate exact FCI energy (reference for energy differences)
        print("  • Calculating exact FCI energy...")
        cisolver = fci.FCI(mole, mf.mo_coeff)
        exact_fci_energy, fci_civec = cisolver.kernel()
        
        print(f"  • Exact FCI energy: {exact_fci_energy:.8f} Ha")
        
        # Create Hartree-Fock initial state
        hf_circuit = QuantumCircuit(n_qubits_actual)
        for i in range(n_electrons_actual):
            hf_circuit.add_X_gate(i)
        hf_state = GeneralCircuitQuantumState(n_qubits_actual, hf_circuit)
        
        # Store system data
        self.hamiltonian = hamiltonian
        self.hf_state = hf_state
        self.exact_fci_energy = exact_fci_energy
        self.sampler = create_qulacs_vector_concurrent_sampler()
        
        print(f"  • Hamiltonian: {len(hamiltonian)} terms")
        print(f"  • State space: 2^{n_qubits_actual} = {2**n_qubits_actual:,} states")
        print()
        
        return hamiltonian, hf_state, exact_fci_energy
    
    def run_time_average_te_qsci(self, time_start, time_end):
        """
        Run time-average TE-QSCI for a given time interval using exact evolution.
        
        Args:
            time_start: Start time of interval
            time_end: End time of interval
            
        Returns:
            Tuple of (qsci_energy, num_selected_states, energy_difference)
        """
        print(f"Running time-average TE-QSCI for interval [{time_start}, {time_end}]...")
        
        # Create time points for averaging
        if time_start == time_end:
            # Single time point
            evolution_times = [time_start]
            interval_label = f"t = {time_start}"
        else:
            # Multiple time points for averaging
            evolution_times = np.linspace(time_start, time_end, self.num_time_points)
            interval_label = f"t ∈ [{time_start}, {time_end}]"
        
        print(f"  • Time points: {len(evolution_times)} points")
        print(f"  • Evolution times: {evolution_times[0]:.2f} to {evolution_times[-1]:.2f}")
        
        # Create Time-Average TE-QSCI algorithm with exact time evolution
        te_qsci_algorithm = create_qsci_algorithm(
            QSCIVariant.TIME_AVERAGE_TE,
            self.hamiltonian,
            sampler=self.sampler,
            evolution_times=evolution_times,
            num_states_pick_out=self.subspace_size_R,
            time_evolution_method="exact"  # Use exact time evolution
        )
        
        print(f"  • QSCI subspace size: R = {self.subspace_size_R}")
        print(f"  • Shots per time point: {self.shots_per_time}")
        
        # Run the algorithm
        start_time = time.time()
        
        try:
            result = te_qsci_algorithm.run(
                initial_state=self.hf_state,
                shots_per_time=self.shots_per_time
            )
            
            execution_time = time.time() - start_time
            
            # Extract results
            qsci_energy = result.ground_state_energy
            num_selected_states = result.subspace_dimension
            energy_difference = qsci_energy - self.exact_fci_energy
            
            print(f"  • E_QSCI: {qsci_energy:.8f} Ha")
            print(f"  • E_exact: {self.exact_fci_energy:.8f} Ha") 
            print(f"  • E_QSCI - E_exact: {energy_difference:.8f} Ha ({energy_difference*1000:.3f} mHa)")
            print(f"  • Selected states: {num_selected_states}")
            print(f"  • Execution time: {execution_time:.1f}s")
            print()
            
            return qsci_energy, num_selected_states, energy_difference, interval_label
            
        except Exception as e:
            print(f"  ✗ Error during TE-QSCI calculation: {e}")
            print()
            return None, None, None, interval_label
    
    def run_full_analysis(self):
        """
        Run the complete Table III analysis for all time intervals.
        
        Returns:
            List of result dictionaries with analysis data
        """
        print("STARTING FULL TABLE III ANALYSIS")
        print("-" * 50)
        
        # Setup molecular system
        self.setup_h6_system()
        
        # Run analysis for each time interval
        results = []
        
        for i, (t_start, t_end) in enumerate(self.time_intervals, 1):
            print(f"Analysis {i}/{len(self.time_intervals)}: ", end="")
            
            qsci_energy, num_states, energy_diff, interval_label = self.run_time_average_te_qsci(t_start, t_end)
            
            if qsci_energy is not None:
                result = {
                    "interval": interval_label,
                    "time_start": t_start,
                    "time_end": t_end,
                    "E_QSCI": qsci_energy,
                    "E_exact": self.exact_fci_energy,
                    "energy_difference_Ha": energy_diff,
                    "energy_difference_mHa": energy_diff * 1000,
                    "num_selected_states": num_states,
                    "success": True
                }
            else:
                result = {
                    "interval": interval_label,
                    "time_start": t_start,
                    "time_end": t_end,
                    "success": False
                }
            
            results.append(result)
        
        return results
    
    def format_table_iii_results(self, results):
        """
        Format results in Table III style for publication.
        
        Args:
            results: List of result dictionaries from analysis
        """
        print("TABLE III REPRODUCTION: TIME-AVERAGE TE-QSCI ENERGY DIFFERENCES")
        print("=" * 80)
        print(f"System: H6 linear chain, {self.basis_set.upper()} basis")
        print(f"Method: Time-average TE-QSCI with exact time evolution")
        print(f"Reference: Exact FCI energy = {self.exact_fci_energy:.8f} Ha")
        print()
        
        # Table header
        print(f"{'Time Interval':<18} | {'N_det':>6} | {'E_QSCI - E_exact (mHa)':>25}")
        print("-" * 60)
        
        # Table data
        successful_results = [r for r in results if r.get("success", False)]
        
        for result in successful_results:
            interval = result["interval"]
            n_det = result["num_selected_states"]
            energy_diff_mha = result["energy_difference_mHa"]
            
            print(f"{interval:<18} | {n_det:>6} | {energy_diff_mha:>25.3f}")
        
        print()
        
        # Summary statistics
        if successful_results:
            energy_diffs = [r["energy_difference_mHa"] for r in successful_results]
            n_dets = [r["num_selected_states"] for r in successful_results]
            
            print("SUMMARY STATISTICS:")
            print(f"  • Successful calculations: {len(successful_results)}/{len(results)}")
            print(f"  • Energy differences range: {min(energy_diffs):.3f} to {max(energy_diffs):.3f} mHa")
            print(f"  • Average |energy difference|: {np.mean(np.abs(energy_diffs)):.3f} mHa")
            print(f"  • Selected states range: {min(n_dets)} to {max(n_dets)}")
            print()
        
        # Failed calculations
        failed_results = [r for r in results if not r.get("success", False)]
        if failed_results:
            print("FAILED CALCULATIONS:")
            for result in failed_results:
                print(f"  ✗ {result['interval']}")
            print()


def main():
    """Main function to run Table III reproduction analysis."""
    
    analysis = H6TableIIIAnalysis()
    
    try:
        # Run complete analysis
        results = analysis.run_full_analysis()
        
        # Format and display results
        analysis.format_table_iii_results(results)
        
        print("✅ Table III reproduction completed successfully!")
        
        return results
        
    except Exception as e:
        print(f"✗ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    results = main()