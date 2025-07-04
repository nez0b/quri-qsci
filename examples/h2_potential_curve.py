"""
H2 Potential Energy Curve Calculation using QSCI and TE-QSCI

This example calculates the H2 molecule potential energy curve using different methods:
1. VQE (reference from quantum chemistry introduction)
2. Vanilla QSCI 
3. Single-time TE-QSCI
4. Time-average TE-QSCI

The results are compared and plotted on the same graph to demonstrate the
effectiveness of TE-QSCI methods compared to VQE for quantum chemistry.

Based on: 2_quantum-chemistry/0_introduction/0_introduction.ipynb
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time


# PySCF and quantum chemistry imports
from pyscf import scf, gto
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.ansatz import KUpCCGSD
from quri_parts.openfermion.transforms import jordan_wigner

# QURI Parts core imports
from quri_parts.core.operator import Operator
from quri_parts.core.state import apply_circuit, ParametricCircuitQuantumState
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit

# QURI Parts estimator and VQE imports
from quri_parts.qulacs.estimator import (
    create_qulacs_vector_parametric_estimator,
    create_qulacs_vector_concurrent_parametric_estimator
)
from quri_parts.core.estimator.gradient import create_parameter_shift_gradient_estimator
from quri_parts.algo.optimizer import OptimizerStatus, LBFGS

# QURI Parts sampler imports
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler

# Our QSCI implementations
from quri_qsci import (
    create_qsci_algorithm, QSCIVariant,
    VanillaQSCIAlgorithm, SingleTimeTeQSCIAlgorithm, TimeAverageTeQSCIAlgorithm
)


class H2PotentialCurveStudy:
    """Comprehensive H2 potential energy curve study using multiple methods."""
    
    def __init__(self):
        """Initialize the H2 potential curve study."""
        self.bond_lengths = np.linspace(0.5, 2.5, 15)  # Reduced range for demonstration
        self.results = {
            'bond_lengths': self.bond_lengths,
            'vqe_energies': [],
            'vanilla_qsci_energies': [],
            'single_time_te_qsci_energies': [],
            'time_average_te_qsci_energies': [],
            'hartree_fock_energies': []  # Reference HF energies
        }
        
        # QSCI parameters
        self.qsci_shots = 2000  # Reduced shots for faster execution
        self.num_states_pick_out = 32  # Reduced states for demonstration
        self.te_evolution_time = 1.0
        self.te_evolution_times = [0.5, 1.0, 1.5]
        
        # Create sampler for QSCI methods
        self.sampler = create_qulacs_vector_concurrent_sampler()
        
        print("H2 Potential Energy Curve Study")
        print("=" * 50)
        print(f"Bond lengths: {self.bond_lengths[0]:.1f} - {self.bond_lengths[-1]:.1f} Å ({len(self.bond_lengths)} points)")
        print(f"QSCI shots: {self.qsci_shots}")
        print(f"States to pick: {self.num_states_pick_out}")
        print(f"TE evolution time: {self.te_evolution_time}")
        print(f"Time-average times: {self.te_evolution_times}")
        print()
    
    def get_h2_hamiltonian_and_hf_state(self, r: float) -> tuple[Operator, GeneralCircuitQuantumState, float]:
        """
        Generate H2 Hamiltonian and Hartree-Fock state for given bond length.
        
        Args:
            r: Bond length in Angstrom
            
        Returns:
            Tuple of (hamiltonian, hf_state, hf_energy)
        """
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
    
    def get_h2_vqe_ansatz(self, r: float) -> tuple[Operator, ParametricCircuitQuantumState]:
        """
        Generate H2 Hamiltonian and VQE ansatz state (from original notebook).
        
        Args:
            r: Bond length in Angstrom
            
        Returns:
            Tuple of (hamiltonian, ansatz_state)
        """
        mole = gto.M(atom=f"H 0 0 0; H 0 0 {r}", basis="sto-3g", verbose=0)
        mf = scf.RHF(mole).run(verbose=0)
        
        hamiltonian, mapping = get_qubit_mapped_hamiltonian(
            *get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
        )

        ansatz_circuit = KUpCCGSD(2*mole.nao, k=1, fermion_qubit_mapping=mapping)
        hf_state = mapping.state_mapper([0, 1])
        ansatz_state = apply_circuit(ansatz_circuit, hf_state)
        
        return hamiltonian, ansatz_state
    
    def run_vqe_single_point(self, r: float) -> float:
        """Run VQE for a single bond length (from original notebook)."""
        
        parametric_estimator = create_qulacs_vector_parametric_estimator()
        cp_estimator = create_qulacs_vector_concurrent_parametric_estimator()
        gradient_estimator = create_parameter_shift_gradient_estimator(cp_estimator)
        
        hamiltonian, ansatz_state = self.get_h2_vqe_ansatz(r)
        
        # Define cost and gradient functions
        cost = lambda param: parametric_estimator(hamiltonian, ansatz_state, param).value.real
        grad = lambda param: np.array(gradient_estimator(hamiltonian, ansatz_state, param).values).real
        
        # VQE optimization
        init_param = np.random.random(ansatz_state.parametric_circuit.parameter_count)
        optimizer = LBFGS()
        op_state = optimizer.get_init_state(init_param)
        
        while True:
            op_state = optimizer.step(op_state, cost, grad) 
            if op_state.status == OptimizerStatus.CONVERGED or op_state.status == OptimizerStatus.FAILED:
                break
        
        return op_state.cost
    
    def run_vanilla_qsci_single_point(self, r: float) -> float:
        """Run Vanilla QSCI for a single bond length."""
        
        hamiltonian, hf_state, hf_energy = self.get_h2_hamiltonian_and_hf_state(r)
        
        # Create Vanilla QSCI algorithm
        algorithm = create_qsci_algorithm(
            QSCIVariant.VANILLA,
            hamiltonian,
            sampler=self.sampler,
            num_states_pick_out=self.num_states_pick_out
        )
        
        # Run algorithm
        result = algorithm.run(
            input_states=[hf_state],
            total_shots=self.qsci_shots
        )
        
        return result.ground_state_energy
    
    def run_single_time_te_qsci_single_point(self, r: float) -> float:
        """Run Single-time TE-QSCI for a single bond length."""
        
        hamiltonian, hf_state, hf_energy = self.get_h2_hamiltonian_and_hf_state(r)
        
        # Create Single-time TE-QSCI algorithm
        algorithm = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            sampler=self.sampler,
            evolution_time=self.te_evolution_time,
            num_states_pick_out=self.num_states_pick_out
        )
        
        # Run algorithm
        result = algorithm.run(
            initial_state=hf_state,
            total_shots=self.qsci_shots
        )
        
        return result.ground_state_energy
    
    def run_time_average_te_qsci_single_point(self, r: float) -> float:
        """Run Time-average TE-QSCI for a single bond length."""
        
        hamiltonian, hf_state, hf_energy = self.get_h2_hamiltonian_and_hf_state(r)
        
        # Create Time-average TE-QSCI algorithm
        algorithm = create_qsci_algorithm(
            QSCIVariant.TIME_AVERAGE_TE,
            hamiltonian,
            sampler=self.sampler,
            evolution_times=self.te_evolution_times,
            num_states_pick_out=self.num_states_pick_out
        )
        
        # Run algorithm
        shots_per_time = self.qsci_shots // len(self.te_evolution_times)
        result = algorithm.run(
            initial_state=hf_state,
            shots_per_time=shots_per_time
        )
        
        return result.ground_state_energy
    
    def run_full_potential_curve_study(self):
        """Run the complete potential energy curve study with all methods."""
        
        print("Running H2 potential energy curve calculations...")
        print("-" * 50)
        
        for i, r in enumerate(self.bond_lengths):
            print(f"Bond length {r:.2f} Å ({i+1}/{len(self.bond_lengths)}):")
            
            # Get Hartree-Fock reference energy
            _, _, hf_energy = self.get_h2_hamiltonian_and_hf_state(r)
            self.results['hartree_fock_energies'].append(hf_energy)
            print(f"  HF energy: {hf_energy:.6f} Ha")
            
            # Run VQE
            try:
                start_time = time.time()
                vqe_energy = self.run_vqe_single_point(r)
                vqe_time = time.time() - start_time
                self.results['vqe_energies'].append(vqe_energy)
                print(f"  VQE energy: {vqe_energy:.6f} Ha ({vqe_time:.1f}s)")
            except Exception as e:
                print(f"  VQE failed: {e}")
                self.results['vqe_energies'].append(np.nan)
            
            # Run Vanilla QSCI
            try:
                start_time = time.time()
                vanilla_energy = self.run_vanilla_qsci_single_point(r)
                vanilla_time = time.time() - start_time
                self.results['vanilla_qsci_energies'].append(vanilla_energy)
                print(f"  Vanilla QSCI: {vanilla_energy:.6f} Ha ({vanilla_time:.1f}s)")
            except Exception as e:
                print(f"  Vanilla QSCI failed: {e}")
                self.results['vanilla_qsci_energies'].append(np.nan)
            
            # Run Single-time TE-QSCI
            try:
                start_time = time.time()
                st_te_energy = self.run_single_time_te_qsci_single_point(r)
                st_te_time = time.time() - start_time
                self.results['single_time_te_qsci_energies'].append(st_te_energy)
                print(f"  Single-time TE-QSCI: {st_te_energy:.6f} Ha ({st_te_time:.1f}s)")
            except Exception as e:
                print(f"  Single-time TE-QSCI failed: {e}")
                self.results['single_time_te_qsci_energies'].append(np.nan)
            
            # Run Time-average TE-QSCI
            try:
                start_time = time.time()
                ta_te_energy = self.run_time_average_te_qsci_single_point(r)
                ta_te_time = time.time() - start_time
                self.results['time_average_te_qsci_energies'].append(ta_te_energy)
                print(f"  Time-avg TE-QSCI: {ta_te_energy:.6f} Ha ({ta_te_time:.1f}s)")
            except Exception as e:
                print(f"  Time-average TE-QSCI failed: {e}")
                self.results['time_average_te_qsci_energies'].append(np.nan)
            
            print()
    
    def analyze_results(self):
        """Analyze and print summary of results."""
        print("Results Summary:")
        print("=" * 70)
        
        # Find minimum energies (equilibrium bond lengths)
        methods = {
            'Hartree-Fock': self.results['hartree_fock_energies'],
            'VQE': self.results['vqe_energies'], 
            'Vanilla QSCI': self.results['vanilla_qsci_energies'],
            'Single-time TE-QSCI': self.results['single_time_te_qsci_energies'],
            'Time-average TE-QSCI': self.results['time_average_te_qsci_energies']
        }
        
        print(f"{'Method':<20} {'Min Energy (Ha)':<15} {'Eq. Bond (Å)':<15} {'Diff from VQE (mHa)':<20}")
        print("-" * 70)
        
        vqe_min = np.nanmin(self.results['vqe_energies'])
        
        for method_name, energies in methods.items():
            energies_array = np.array(energies)
            valid_mask = ~np.isnan(energies_array)
            
            if np.any(valid_mask):
                min_energy = np.nanmin(energies_array)
                min_idx = np.nanargmin(energies_array)
                eq_bond = self.bond_lengths[min_idx]
                diff_from_vqe = (min_energy - vqe_min) * 1000  # Convert to milliHartree
                
                print(f"{method_name:<20} {min_energy:<15.6f} {eq_bond:<15.2f} {diff_from_vqe:<20.2f}")
            else:
                print(f"{method_name:<20} {'N/A':<15} {'N/A':<15} {'N/A':<20}")
        
        print()
    
    def plot_potential_curves(self, save_plot=True):
        """Plot and compare potential energy curves."""
        
        # Set up the plot style (matching the original notebook)
        sns.set_theme("talk")
        plt.rcParams["figure.figsize"] = (16, 10)
        plt.rcParams["font.size"] = 14
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
        
        # Main plot
        bond_lengths = self.results['bond_lengths']
        
        # Plot all methods
        methods_to_plot = [
            ('Hartree-Fock', self.results['hartree_fock_energies'], 'gray', '--', 'Reference HF'),
            ('VQE', self.results['vqe_energies'], 'blue', '-o', 'VQE (1-UpCCGSD)'),
            ('Vanilla QSCI', self.results['vanilla_qsci_energies'], 'green', '-s', f'Vanilla QSCI (R={self.num_states_pick_out})'),
            ('Single-time TE-QSCI', self.results['single_time_te_qsci_energies'], 'red', '-^', f'Single-time TE-QSCI (t={self.te_evolution_time})'),
            ('Time-average TE-QSCI', self.results['time_average_te_qsci_energies'], 'purple', '-d', f'Time-avg TE-QSCI (t={self.te_evolution_times})')
        ]
        
        for method_name, energies, color, style, label in methods_to_plot:
            energies_array = np.array(energies)
            valid_mask = ~np.isnan(energies_array)
            
            if np.any(valid_mask):
                ax1.plot(bond_lengths[valid_mask], energies_array[valid_mask], 
                        style, color=color, label=label, markersize=6, linewidth=2)
        
        ax1.set_title("H₂ Potential Energy Curves: VQE vs QSCI Methods", fontsize=18, fontweight='bold')
        ax1.set_xlabel("Bond Length (Å)", fontsize=14)
        ax1.set_ylabel("Energy (Hartree)", fontsize=14)
        ax1.legend(fontsize=12, loc='upper right')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-1.2, -0.9)  # Focus on relevant energy range
        
        # Difference plot (relative to VQE)
        vqe_energies = np.array(self.results['vqe_energies'])
        
        methods_for_diff = [
            ('Vanilla QSCI', self.results['vanilla_qsci_energies'], 'green', '-s'),
            ('Single-time TE-QSCI', self.results['single_time_te_qsci_energies'], 'red', '-^'),
            ('Time-average TE-QSCI', self.results['time_average_te_qsci_energies'], 'purple', '-d')
        ]
        
        for method_name, energies, color, style in methods_for_diff:
            energies_array = np.array(energies)
            diff = (energies_array - vqe_energies) * 1000  # Convert to milliHartree
            valid_mask = ~np.isnan(diff)
            
            if np.any(valid_mask):
                ax2.plot(bond_lengths[valid_mask], diff[valid_mask], 
                        style, color=color, label=f'{method_name} - VQE', markersize=6, linewidth=2)
        
        ax2.axhline(y=0, color='blue', linestyle='--', alpha=0.7, label='VQE reference')
        ax2.set_title("Energy Differences Relative to VQE", fontsize=16, fontweight='bold')
        ax2.set_xlabel("Bond Length (Å)", fontsize=14)
        ax2.set_ylabel("Energy Difference (mHa)", fontsize=14)
        ax2.legend(fontsize=12)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plot:
            # Create figures directory if it doesn't exist
            figures_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
            os.makedirs(figures_dir, exist_ok=True)
            
            filename = os.path.join(figures_dir, 'h2_potential_curve_comparison.png')
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {filename}")
        
        plt.show()
    
    def save_results(self, filename='h2_potential_curve_results.npz'):
        """Save results to file for later analysis."""
        # Create figures directory if it doesn't exist
        figures_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
        os.makedirs(figures_dir, exist_ok=True)
        
        filepath = os.path.join(figures_dir, filename)
        np.savez(filepath, **self.results)
        print(f"Results saved to: {filepath}")
    
    def run_complete_study(self):
        """Run the complete H2 potential curve study."""
        start_time = time.time()
        
        # Run calculations
        self.run_full_potential_curve_study()
        
        # Analyze results
        self.analyze_results()
        
        # Plot results
        self.plot_potential_curves()
        
        # Save results
        self.save_results()
        
        total_time = time.time() - start_time
        print(f"\nTotal computation time: {total_time:.1f} seconds")
        print("Study completed successfully!")


def main():
    """Main function to run the H2 potential curve study."""
    
    print("H₂ Molecule Potential Energy Curve Study")
    print("Comparing VQE, Vanilla QSCI, and TE-QSCI methods")
    print("=" * 60)
    print()
    
    # Create and run study
    study = H2PotentialCurveStudy()
    study.run_complete_study()


if __name__ == "__main__":
    main()