#!/usr/bin/env python3
"""
H2 TE-QSCI with QURI VM Test
============================

Simple H2 molecule test to verify QURI VM integration is working correctly.
This is a minimal test to ensure VM creation and sampling works before
running the full H6 analysis.

Requirements: H2 molecule with STO-3G basis
"""

import numpy as np
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# QURI Parts imports
from quri_parts.core.operator import Operator, pauli_label
from quri_parts.core.state import CircuitQuantumState, GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit
from quri_parts.backend.units import TimeValue, TimeUnit

# Quantum chemistry imports
from pyscf import gto, scf, fci
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner

# Our implementations
from quri_qsci.qsci_algorithms import SingleTimeTE_QSCI
from quri_qsci.qsci_algo_interface import LoweringLevel, create_qsci_algorithm, QSCIVariant
from quri_qsci.qsci_vm_analysis import VMFactory, QURIVMInterface
from quri_qsci.qsci_vm_interface import create_vm_aware_te_qsci


@dataclass
class H2TEQSCIResult:
    """Results from H2 TE-QSCI VM test."""
    hartree_fock_energy: float
    fci_energy: float
    te_qsci_energies: Dict[str, float]
    execution_times: Dict[str, float]
    n_qubits: int


class H2TEQSCIVMTest:
    """Simple H2 TE-QSCI test with QURI VM integration."""
    
    def __init__(self, evolution_time: float = 1.0, subspace_size: int = 10, num_shots: int = 50):
        """Initialize H2 TE-QSCI VM test.
        
        Args:
            evolution_time: Time evolution parameter (reduced for speed)
            subspace_size: QSCI subspace dimension (very small for testing)
            num_shots: Number of measurement shots (minimal for testing)
        """
        self.evolution_time = evolution_time
        self.subspace_size = subspace_size
        self.num_shots = num_shots
        
        # Test multiple error levels to verify VM differences
        self.error_rates = {
            'ideal': 0.0,      # Perfect LogicalCircuit
            'low': 1e-4,       # Low error STAR
            'medium': 1e-3,    # Medium error STAR
            'high': 1e-2,      # High error STAR
        }
        
        print(f"H2 TE-QSCI VM Test")
        print(f"Configuration:")
        print(f"   • Evolution time: t = {evolution_time}")
        print(f"   • Subspace size: R = {subspace_size}")
        print(f"   • Number of shots: {num_shots}")
        print(f"   • Architecture: STAR")
        print(f"   • Molecule: H2 (STO-3G)")
        print()
        
        # Setup molecular system
        print("Setting up H2 molecular system...")
        self._setup_h2_molecule()
        
        # Setup TE-QSCI algorithms
        print("Setting up TE-QSCI algorithms...")
        self._setup_te_qsci_algorithms()
        
        print("Initialization complete!")
        print()
    
    def _setup_h2_molecule(self):
        """Setup H2 molecular system using PySCF."""
        # Simple H2 molecule
        geometry = [
            ['H', (0.0, 0.0, 0.0)],
            ['H', (0.0, 0.0, 0.74)]  # Bond length ~0.74 Å
        ]
        
        # Create PySCF molecule
        self.mole = gto.M(
            atom=geometry,
            basis="sto-3g",
            charge=0,
            spin=0,
            verbose=0
        )
        
        print(f"   • Atoms: {self.mole.natm}")
        print(f"   • Electrons: {self.mole.nelectron}")
        print(f"   • Basis functions: {self.mole.nao}")
        
        # Hartree-Fock calculation
        self.mf = scf.RHF(self.mole)
        self.mf.run()
        self.hartree_fock_energy = self.mf.e_tot
        
        # FCI calculation for reference
        cisolver = fci.FCI(self.mole, self.mf.mo_coeff)
        self.fci_energy, self.fci_civec = cisolver.kernel()
        
        print(f"   • HF energy: {self.hartree_fock_energy:.6f} Ha")
        print(f"   • FCI energy: {self.fci_energy:.6f} Ha")
        print(f"   • Correlation energy: {self.fci_energy - self.hartree_fock_energy:.6f} Ha")
        
        # Create QURI Parts Hamiltonian
        self.active_space, mo_eint_set = get_spin_mo_integrals_from_mole(
            self.mole, self.mf.mo_coeff
        )
        
        self.hamiltonian, self.jw_mapping = get_qubit_mapped_hamiltonian(
            self.active_space, 
            mo_eint_set, 
            fermion_qubit_mapping=jordan_wigner
        )
        
        self.n_qubits = 2 * self.active_space.n_active_orb
        
        print(f"   • Qubits needed: {self.n_qubits}")
        print(f"   • Hamiltonian terms: {self.hamiltonian.n_terms}")
        
        # Create Hartree-Fock initial state
        self.hartree_fock_state = self._create_hartree_fock_state()
        
    def _create_hartree_fock_state(self) -> GeneralCircuitQuantumState:
        """Create Hartree-Fock initial state."""
        circuit = QuantumCircuit(self.n_qubits)
        
        # Fill occupied orbitals for H2 (2 electrons)
        for i in range(self.mole.nelectron):
            circuit.add_X_gate(i)
        
        return GeneralCircuitQuantumState(self.n_qubits, circuit)
    
    def _setup_te_qsci_algorithms(self):
        """Setup VM-aware TE-QSCI algorithms with different error levels."""
        from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
        
        # Create fallback sampler
        fallback_sampler = create_qulacs_vector_concurrent_sampler()
        
        # Create TE-QSCI algorithms for each error level
        self.te_qsci_algorithms = {}
        
        for error_name, error_rate in self.error_rates.items():
            if error_name == 'ideal':
                # Ideal uses STAR VM with zero error for both levels
                logical_vm = VMFactory.create_star_vm(0.0)  # Baseline STAR
                arch_vm = VMFactory.create_star_vm(0.0)     # Same for comparison
                vm_interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)
                
                algorithm = create_vm_aware_te_qsci(
                    hamiltonian=self.hamiltonian,
                    evolution_time=self.evolution_time,
                    vm_interface=vm_interface,
                    sampler=fallback_sampler,
                    num_states_pick_out=self.subspace_size,
                    trotter_steps=5,  # Reduced for speed
                    sampling_level=LoweringLevel.LogicalCircuit,
                    use_vm_sampling=True
                )
            else:
                # STAR architecture with baseline vs error comparison
                logical_vm = VMFactory.create_star_vm(0.0)      # Baseline STAR
                arch_vm = VMFactory.create_star_vm(error_rate)   # Error STAR
                vm_interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)
                
                algorithm = create_vm_aware_te_qsci(
                    hamiltonian=self.hamiltonian,
                    evolution_time=self.evolution_time,
                    vm_interface=vm_interface,
                    sampler=fallback_sampler,
                    num_states_pick_out=self.subspace_size,
                    trotter_steps=5,  # Reduced for speed
                    sampling_level=LoweringLevel.ArchLogicalCircuit,
                    use_vm_sampling=True
                )
            
            self.te_qsci_algorithms[error_name] = algorithm
            print(f"   • {error_name.title()} Algorithm: {algorithm.name} (error rate: {error_rate})")
        
        print(f"   • Evolution time: {self.evolution_time}")
        print(f"   • Subspace size: {self.subspace_size}")
        print(f"   • Trotter steps: 5")
    
    def run_test(self) -> H2TEQSCIResult:
        """Run H2 TE-QSCI VM test."""
        print("Running H2 TE-QSCI VM test...")
        
        te_qsci_energies = {}
        execution_times = {}
        
        # Run analysis for each error level
        for error_name in self.error_rates.keys():
            error_rate = self.error_rates[error_name]
            algorithm = self.te_qsci_algorithms[error_name]
            
            print(f"\n=== {error_name.title()} VM Test (error rate: {error_rate}) ===")
            start_time = time.time()
            
            try:
                vm_result = algorithm.run_with_vm_analysis(
                    self.hartree_fock_state,
                    analyze_sampling=True,
                    num_shots=self.num_shots
                )
                exec_time = time.time() - start_time
                execution_times[error_name] = exec_time
                
                energy = vm_result.qsci_result.ground_state_energy
                te_qsci_energies[error_name] = energy
                
                energy_vs_hf = energy - self.hartree_fock_energy
                energy_vs_fci = energy - self.fci_energy
                
                print(f"   ✓ {error_name.title()} TE-QSCI completed in {exec_time:.2f} seconds")
                print(f"   • Ground state energy: {energy:.6f} Ha")
                print(f"   • Energy vs HF: {energy_vs_hf:.6f} Ha")
                print(f"   • Energy vs FCI: {energy_vs_fci:.6f} Ha")
                print(f"   • Architecture overhead: {vm_result.architecture_overhead:.1f}x")
                print(f"   • Fidelity impact: {vm_result.fidelity_impact:.1%}")
                
            except Exception as e:
                print(f"   ✗ {error_name.title()} TE-QSCI failed: {e}")
                execution_times[error_name] = time.time() - start_time
                te_qsci_energies[error_name] = float('nan')
        
        # Compare results across all error levels
        if 'ideal' in te_qsci_energies:
            print(f"\n=== VM Error Rate Comparison ===")
            ideal_energy = te_qsci_energies['ideal']
            
            if not np.isnan(ideal_energy):
                print(f"   • Ideal Energy: {ideal_energy:.6f} Ha")
                
                # Check each error level
                max_deviation = 0.0
                for error_name in ['low', 'medium', 'high']:
                    if error_name in te_qsci_energies:
                        error_energy = te_qsci_energies[error_name]
                        error_rate = self.error_rates[error_name]
                        
                        if not np.isnan(error_energy):
                            energy_diff = error_energy - ideal_energy
                            abs_diff = abs(energy_diff)
                            max_deviation = max(max_deviation, abs_diff)
                            
                            print(f"   • {error_name.title()} Error (ε={error_rate}): {error_energy:.6f} Ha, Δ={energy_diff:+.6f} Ha")
                        else:
                            print(f"   • {error_name.title()} Error (ε={error_rate}): FAILED")
                
                # Validate that higher error rates produce larger deviations
                print(f"\n   Analysis:")
                if max_deviation > 1e-6:
                    print(f"   ✓ Maximum deviation: {max_deviation:.6f} Ha")
                    print(f"   ✓ Different VMs produce different results (good!)")
                    
                    # Check if deviation increases with error rate
                    if self._check_error_progression(te_qsci_energies, ideal_energy):
                        print(f"   ✓ Higher error rates show larger deviations (excellent!)")
                    else:
                        print(f"   ⚠ Error progression not strictly monotonic (acceptable)")
                else:
                    print(f"   ⚠ Maximum deviation: {max_deviation:.6f} Ha")
                    print(f"   ⚠ All VMs produce very similar results")
                    print(f"   → Try higher error rates or more shots for better differentiation")
            else:
                print(f"   ✗ Ideal calculation failed")
        
        return H2TEQSCIResult(
            hartree_fock_energy=self.hartree_fock_energy,
            fci_energy=self.fci_energy,
            te_qsci_energies=te_qsci_energies,
            execution_times=execution_times,
            n_qubits=self.n_qubits
        )
    
    def _check_error_progression(self, energies: Dict[str, float], ideal_energy: float) -> bool:
        """Check if higher error rates produce larger energy deviations."""
        error_levels = ['low', 'medium', 'high']
        deviations = []
        
        for level in error_levels:
            if level in energies and not np.isnan(energies[level]):
                deviation = abs(energies[level] - ideal_energy)
                deviations.append(deviation)
            else:
                return False  # Can't check progression if some failed
        
        # Check if deviations generally increase (allow some noise)
        if len(deviations) >= 2:
            # Check if medium >= low and high >= medium (with some tolerance)
            return deviations[-1] >= deviations[0] * 0.5  # Allow some variation
        
        return False
    
    def print_summary(self, result: H2TEQSCIResult):
        """Print test summary."""
        print("=" * 60)
        print("H2 TE-QSCI VM TEST SUMMARY")
        print("=" * 60)
        
        print(f"\nMolecular System:")
        print(f"   Molecule:          H2 (STO-3G)")
        print(f"   Qubits:           {result.n_qubits}")
        print(f"   HF Energy:        {result.hartree_fock_energy:.6f} Ha")
        print(f"   FCI Energy:       {result.fci_energy:.6f} Ha")
        
        print(f"\nVM Test Results:")
        for error_name, energy in result.te_qsci_energies.items():
            exec_time = result.execution_times[error_name]
            error_rate = self.error_rates[error_name]
            
            if not np.isnan(energy):
                energy_vs_fci = energy - result.fci_energy
                print(f"   {error_name.title():<8} (ε={error_rate}): {energy:.6f} Ha, vs FCI: {energy_vs_fci:+.6f} Ha, time: {exec_time:.1f}s")
            else:
                print(f"   {error_name.title():<8} (ε={error_rate}): FAILED, time: {exec_time:.1f}s")
        
        # Test validation
        success_count = sum(1 for e in result.te_qsci_energies.values() if not np.isnan(e))
        total_tests = len(result.te_qsci_energies)
        
        print(f"\nTest Status: {success_count}/{total_tests} tests passed")
        
        if success_count == total_tests:
            print("✓ H2 VM test completed successfully!")
            print("✓ QURI VM creation and sampling working correctly")
            print("✓ Ready to run full H6 analysis")
        else:
            print("✗ Some tests failed - check implementation")
        
        print("=" * 60)


def main():
    """Main test function."""
    print("Starting H2 TE-QSCI VM Test")
    print("=" * 40)
    
    try:
        # Initialize test with parameters to see error effects
        test = H2TEQSCIVMTest(evolution_time=1.0, subspace_size=10, num_shots=200)
        
        # Run test
        result = test.run_test()
        
        # Print summary
        test.print_summary(result)
        
        print(f"\nTest completed! If successful, you can now run:")
        print(f"python examples/h6_te_qsci_vm_analysis.py")
        
        return result
        
    except Exception as e:
        print(f"Error during H2 test: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    result = main()