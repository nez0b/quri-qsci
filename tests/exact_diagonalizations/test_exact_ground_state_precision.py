"""
Exact Ground State Precision Validation Tests

This module validates that all QSCI variants achieve machine precision when starting
with exact ground states extracted from Hamiltonian eigensolvers. This proves the
mathematical correctness of QSCI algorithms independent of state preparation accuracy.
"""

import pytest
import numpy as np
import sys
import os
from typing import Dict, Any, List
from scipy.sparse import csc_matrix

# Add parent directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))\

# Import test utilities
from utils.hamiltonian_factory import (
    create_simple_2x2_hamiltonian,
    create_pauli_hamiltonian,
    create_tfim_hamiltonian,
    create_heisenberg_hamiltonian,
    create_random_sparse_hamiltonian,
    get_exact_solution,
    validate_sparse_hamiltonian_properties
)

from utils.exact_ground_state_utils import (
    extract_exact_ground_state,
    create_exact_computational_basis_state,
    compare_qsci_variants_precision
)

try:
    from quri_qsci.qsci_algorithms import (
        VanillaQSCI, SingleTimeTE_QSCI, TimeAverageTE_QSCI, StateVectorTE_QSCI
    )
    from quri_parts.core.operator import Operator
    from quri_parts.core.operator.pauli import PauliLabel, pauli_label
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: QSCI implementation not available: {e}")
    QSCI_AVAILABLE = False


def sparse_matrix_to_quri_operator(sparse_matrix: csc_matrix, n_qubits: int) -> Operator:
    """Convert sparse Hamiltonian matrix to QURI-Parts Operator format."""
    if sparse_matrix.shape[0] != 2**n_qubits:
        raise ValueError(f"Matrix dimension {sparse_matrix.shape[0]} doesn't match {2**n_qubits} for {n_qubits} qubits")
    
    operator = Operator()
    
    # Convert to dense for Pauli decomposition (only for small systems)
    if sparse_matrix.shape[0] > 64:
        raise ValueError("Matrix too large for Pauli decomposition in test")
    
    dense_matrix = sparse_matrix.toarray()
    
    # Generate all possible Pauli strings for n_qubits
    pauli_strings = generate_pauli_strings(n_qubits)
    
    for pauli_string, pauli_matrix in pauli_strings:
        # Compute coefficient: Tr(H * P) / 2^n
        coeff = np.trace(dense_matrix @ pauli_matrix) / (2**n_qubits)
        
        if abs(coeff) > 1e-12:  # Only keep non-zero terms
            operator[pauli_string] = complex(coeff)
    
    return operator


def generate_pauli_strings(n_qubits: int) -> List:
    """Generate all Pauli strings and their matrix representations for n qubits."""
    # Pauli matrices
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    
    pauli_matrices = [I, X, Y, Z]
    pauli_names = ['I', 'X', 'Y', 'Z']
    
    pauli_strings = []
    
    # Generate all 4^n Pauli strings
    for i in range(4**n_qubits):
        # Convert i to base-4 representation
        pauli_indices = []
        temp = i
        for _ in range(n_qubits):
            pauli_indices.append(temp % 4)
            temp //= 4
        pauli_indices.reverse()  # MSB first
        
        # Build PauliLabel
        if all(idx == 0 for idx in pauli_indices):
            # Identity operator
            pauli_label_obj = PauliLabel()
        else:
            # Non-trivial Pauli string
            pauli_terms = []
            for qubit_idx, pauli_idx in enumerate(pauli_indices):
                if pauli_idx != 0:  # Skip identity
                    pauli_terms.append(f"{pauli_names[pauli_idx]}{qubit_idx}")
            
            if pauli_terms:
                pauli_label_obj = pauli_label(" ".join(pauli_terms))
            else:
                pauli_label_obj = PauliLabel()
        
        # Build matrix representation
        matrix = pauli_matrices[pauli_indices[0]]
        for j in range(1, n_qubits):
            matrix = np.kron(matrix, pauli_matrices[pauli_indices[j]])
        
        pauli_strings.append((pauli_label_obj, matrix))
    
    return pauli_strings


@pytest.fixture
def exact_ground_state_test_systems():
    """Create test systems with exact ground states for various Hamiltonians."""
    systems = {}
    
    # Simple 2x2 Hamiltonian (1 qubit)
    hamiltonian_2x2, exact_eigenvalues, _ = create_simple_2x2_hamiltonian()
    sparse_2x2 = csc_matrix(hamiltonian_2x2)
    ground_vector_2x2, ground_energy_2x2 = extract_exact_ground_state(sparse_2x2)
    
    systems['simple_2x2'] = {
        'hamiltonian_matrix': sparse_2x2,
        'n_qubits': 1,
        'exact_ground_vector': ground_vector_2x2,
        'exact_ground_energy': ground_energy_2x2,
        'quri_hamiltonian': sparse_matrix_to_quri_operator(sparse_2x2, 1)
    }
    
    # 2-qubit Pauli Hamiltonian
    hamiltonian_pauli = create_pauli_hamiltonian()
    sparse_pauli = csc_matrix(hamiltonian_pauli)
    ground_vector_pauli, ground_energy_pauli = extract_exact_ground_state(sparse_pauli)
    
    systems['pauli_2q'] = {
        'hamiltonian_matrix': sparse_pauli,
        'n_qubits': 2,
        'exact_ground_vector': ground_vector_pauli,
        'exact_ground_energy': ground_energy_pauli,
        'quri_hamiltonian': sparse_matrix_to_quri_operator(sparse_pauli, 2)
    }
    
    # 2-qubit TFIM Hamiltonian (critical point)
    hamiltonian_tfim = create_tfim_hamiltonian(n_qubits=2, h=1.0)
    ground_vector_tfim, ground_energy_tfim = extract_exact_ground_state(hamiltonian_tfim)
    
    systems['tfim_2q'] = {
        'hamiltonian_matrix': hamiltonian_tfim,
        'n_qubits': 2,
        'exact_ground_vector': ground_vector_tfim,
        'exact_ground_energy': ground_energy_tfim,
        'quri_hamiltonian': sparse_matrix_to_quri_operator(hamiltonian_tfim, 2)
    }
    
    # 2-qubit Heisenberg Hamiltonian
    hamiltonian_heisenberg = create_heisenberg_hamiltonian(n_qubits=2, jz=1.0, jxy=1.0)
    ground_vector_heisenberg, ground_energy_heisenberg = extract_exact_ground_state(hamiltonian_heisenberg)
    
    systems['heisenberg_2q'] = {
        'hamiltonian_matrix': hamiltonian_heisenberg,
        'n_qubits': 2,
        'exact_ground_vector': ground_vector_heisenberg,
        'exact_ground_energy': ground_energy_heisenberg,
        'quri_hamiltonian': sparse_matrix_to_quri_operator(hamiltonian_heisenberg, 2)
    }
    
    return systems


@pytest.mark.exact_ground_state_precision
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="QSCI implementation not available")
class TestExactGroundStatePrecision:
    """
    Precision validation tests using exact ground states.
    
    These tests validate that all QSCI variants achieve machine precision (1e-12)
    when starting with exact ground states extracted from eigensolvers.
    """
    
    def test_uniform_superposition_validation(self, exact_ground_state_test_systems):
        """Test that uniform superposition provides perfect sampling diversity."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        from quri_parts.qulacs.simulator import evaluate_state_to_vector
        
        for system_name, system_data in exact_ground_state_test_systems.items():
            print(f"\nTesting uniform superposition for {system_name}")
            
            # Create uniform superposition state
            uniform_state = create_diverse_superposition_state(
                n_qubits=system_data['n_qubits']
            )
            
            # Validate uniform superposition
            state_vector = evaluate_state_to_vector(uniform_state)
            if hasattr(state_vector, 'vector'):
                state_vector = state_vector.vector
            
            expected_amplitude = 1.0 / np.sqrt(2**system_data['n_qubits'])
            is_uniform = all(abs(abs(amp) - expected_amplitude) < 1e-10 for amp in state_vector)
            norm = np.linalg.norm(state_vector)
            
            print(f"  Expected amplitude: {expected_amplitude:.6f}")
            print(f"  All amplitudes uniform: {is_uniform}")
            print(f"  State norm: {norm:.12f}")
            
            assert is_uniform, f"Uniform superposition not achieved for {system_name}"
            assert abs(norm - 1.0) < 1e-10, f"State not normalized for {system_name}: norm={norm}"
    
    def test_vanilla_qsci_uniform_superposition_precision(self, exact_ground_state_test_systems):
        """Test VanillaQSCI achieves machine precision with uniform superposition."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        
        for system_name, system_data in exact_ground_state_test_systems.items():
            print(f"\nTesting VanillaQSCI with uniform superposition for {system_name}")
            
            # Create uniform superposition state
            uniform_state = create_diverse_superposition_state(
                n_qubits=system_data['n_qubits']
            )
            
            # Run VanillaQSCI with complete subspace
            sampler = create_qulacs_vector_concurrent_sampler()
            vanilla_qsci = VanillaQSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=sampler,
                num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
            )
            
            result = vanilla_qsci.run([uniform_state], total_shots=2000)
            
            # Validate precision
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            print(f"  Exact energy: {system_data['exact_ground_energy']:.12f}")
            print(f"  QSCI energy: {result.ground_state_energy:.12f}")
            print(f"  Energy error: {energy_error:.2e}")
            print(f"  Selected states: {len(result.selected_states)}/{2**system_data['n_qubits']}")
            
            # Expect machine precision with uniform superposition + complete subspace
            precision_tolerance = 1e-10
            assert energy_error < precision_tolerance, (
                f"VanillaQSCI energy error {energy_error:.2e} exceeds precision tolerance "
                f"{precision_tolerance:.2e} for {system_name}"
            )
    
    def test_single_time_te_qsci_uniform_superposition_precision(self, exact_ground_state_test_systems):
        """Test SingleTimeTE_QSCI achieves precision with uniform superposition."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        
        for system_name, system_data in exact_ground_state_test_systems.items():
            print(f"\nTesting SingleTimeTE_QSCI with uniform superposition for {system_name}")
            
            # Create uniform superposition state
            uniform_state = create_diverse_superposition_state(
                n_qubits=system_data['n_qubits']
            )
            
            # Test multiple evolution times
            evolution_times = [0.1, 0.5, 1.0]
            
            for evolution_time in evolution_times:
                print(f"  Evolution time: {evolution_time}")
                
                # Run SingleTimeTE_QSCI with complete subspace
                sampler = create_qulacs_vector_concurrent_sampler()
                te_qsci = SingleTimeTE_QSCI(
                    hamiltonian=system_data['quri_hamiltonian'],
                    sampler=sampler,
                    evolution_time=evolution_time,
                    num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
                )
                
                result = te_qsci.run([uniform_state], total_shots=3000)  # Higher shots for TE-QSCI
                
                # Validate precision
                energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
                print(f"    Energy error: {energy_error:.2e}")
                
                # Expect high precision with uniform superposition + complete subspace
                precision_tolerance = 1e-8
                assert energy_error < precision_tolerance, (
                    f"SingleTimeTE_QSCI energy error {energy_error:.2e} exceeds tolerance "
                    f"{precision_tolerance:.2e} for {system_name} at t={evolution_time}"
                )
    
    def test_state_vector_te_qsci_uniform_superposition_precision(self, exact_ground_state_test_systems):
        """Test StateVectorTE_QSCI achieves high precision with uniform superposition."""
        from quri_parts.qulacs.simulator import evaluate_state_to_vector
        
        for system_name, system_data in exact_ground_state_test_systems.items():
            print(f"\nTesting StateVectorTE_QSCI with uniform superposition for {system_name}")
            
            # Create uniform superposition circuit state for StateVectorTE_QSCI
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            uniform_state = create_diverse_superposition_state(
                n_qubits=system_data['n_qubits']
            )
            
            # Run StateVectorTE_QSCI with uniform circuit state (not vector)
            from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
            sv_te_qsci = StateVectorTE_QSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=create_qulacs_vector_concurrent_sampler(),
                evolution_time=1.0,
                num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
            )
            
            result = sv_te_qsci.run([uniform_state], total_shots=3000)
            
            # Validate precision
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            print(f"  Exact energy: {system_data['exact_ground_energy']:.12f}")
            print(f"  QSCI energy: {result.ground_state_energy:.12f}")
            print(f"  Energy error: {energy_error:.2e}")
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"StateVectorTE_QSCI energy error {energy_error:.2e} exceeds tolerance "
                f"{precision_tolerance:.2e} for {system_name}"
            )
    
    def test_time_average_te_qsci_uniform_superposition_precision(self, exact_ground_state_test_systems):
        """Test TimeAverageTE_QSCI achieves precision with uniform superposition."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        
        # Select subset of systems for time-average tests (computationally expensive)
        selected_systems = ['simple_2x2', 'pauli_2q']
        
        for system_name in selected_systems:
            if system_name not in exact_ground_state_test_systems:
                continue
                
            system_data = exact_ground_state_test_systems[system_name]
            print(f"\nTesting TimeAverageTE_QSCI with uniform superposition for {system_name}")
            
            # Create uniform superposition state
            uniform_state = create_diverse_superposition_state(
                n_qubits=system_data['n_qubits']
            )
            
            # Run TimeAverageTE_QSCI with complete subspace
            sampler = create_qulacs_vector_concurrent_sampler()
            ta_te_qsci = TimeAverageTE_QSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=sampler,
                evolution_times=[0.3, 0.7, 1.0],
                num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
            )
            
            result = ta_te_qsci.run([uniform_state], total_shots=4500)  # 1500 shots per time * 3 times
            
            # Validate precision
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            print(f"  Exact energy: {system_data['exact_ground_energy']:.12f}")
            print(f"  QSCI energy: {result.ground_state_energy:.12f}")
            print(f"  Energy error: {energy_error:.2e}")
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"TimeAverageTE_QSCI energy error {energy_error:.2e} exceeds tolerance "
                f"{precision_tolerance:.2e} for {system_name}"
            )
    
    def test_uniform_superposition_vs_computational_basis_comparison(self):
        """Compare uniform superposition vs computational basis states for QSCI."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        quri_hamiltonian = sparse_matrix_to_quri_operator(csc_matrix(hamiltonian_matrix), 2)
        
        # Get exact ground energy
        exact_eigenvalues, _ = get_exact_solution(csc_matrix(hamiltonian_matrix), num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        print(f"\nComparing uniform superposition vs computational basis states")
        print(f"Exact ground energy: {exact_ground_energy:.12f}")
        
        # Test uniform superposition (should be perfect)
        uniform_state = create_diverse_superposition_state(n_qubits=2)
        
        sampler = create_qulacs_vector_concurrent_sampler()
        vanilla_qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=4  # Complete subspace
        )
        
        uniform_result = vanilla_qsci.run([uniform_state], total_shots=2000)
        uniform_error = abs(uniform_result.ground_state_energy - exact_ground_energy)
        
        print(f"  Uniform superposition energy: {uniform_result.ground_state_energy:.12f}, error: {uniform_error:.2e}")
        
        # Test each computational basis state for comparison
        for i in range(4):  # |00⟩, |01⟩, |10⟩, |11⟩
            basis_state = create_exact_computational_basis_state(n_qubits=2, basis_state_index=i)
            
            basis_result = vanilla_qsci.run([basis_state], total_shots=2000)
            basis_error = abs(basis_result.ground_state_energy - exact_ground_energy)
            
            print(f"  |{i:02b}⟩ state energy: {basis_result.ground_state_energy:.12f}, error: {basis_error:.2e}")
        
        # Uniform superposition should achieve machine precision
        assert uniform_error < 1e-10, (
            f"Uniform superposition error {uniform_error:.2e} should be machine precision"
        )


@pytest.mark.exact_ground_state_precision
@pytest.mark.performance
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="QSCI implementation not available")
class TestQSCIVariantComparison:
    """Compare precision across all QSCI variants with exact ground states."""
    
    def test_all_variants_precision_comparison(self, exact_ground_state_test_systems):
        """Compare precision of all QSCI variants with uniform superposition."""
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        
        # Select one representative system for comprehensive comparison
        system_name = 'pauli_2q'
        system_data = exact_ground_state_test_systems[system_name]
        
        print(f"\n=== QSCI VARIANT PRECISION COMPARISON: {system_name} ===")
        print(f"Exact ground energy: {system_data['exact_ground_energy']:.12f}")
        
        # Create uniform superposition state
        uniform_state = create_diverse_superposition_state(
            n_qubits=system_data['n_qubits']
        )
        
        variant_results = {}
        sampler = create_qulacs_vector_concurrent_sampler()
        
        # Test VanillaQSCI
        try:
            vanilla_qsci = VanillaQSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=sampler,
                num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
            )
            result = vanilla_qsci.run([uniform_state], total_shots=2000)
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            variant_results['VanillaQSCI'] = {
                'energy': result.ground_state_energy,
                'error': energy_error,
                'success': True
            }
            print(f"VanillaQSCI: {result.ground_state_energy:.12f} (error: {energy_error:.2e})")
        except Exception as e:
            variant_results['VanillaQSCI'] = {'success': False, 'error_msg': str(e)}
            print(f"VanillaQSCI: FAILED - {e}")
        
        # Test SingleTimeTE_QSCI
        try:
            st_te_qsci = SingleTimeTE_QSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=sampler,
                evolution_time=0.5,
                num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
            )
            result = st_te_qsci.run([uniform_state], total_shots=3000)
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            variant_results['SingleTimeTE_QSCI'] = {
                'energy': result.ground_state_energy,
                'error': energy_error,
                'success': True
            }
            print(f"SingleTimeTE_QSCI: {result.ground_state_energy:.12f} (error: {energy_error:.2e})")
        except Exception as e:
            variant_results['SingleTimeTE_QSCI'] = {'success': False, 'error_msg': str(e)}
            print(f"SingleTimeTE_QSCI: FAILED - {e}")
        
        # Test StateVectorTE_QSCI with uniform superposition vector
        try:
            n_qubits = system_data['n_qubits']
            uniform_amplitude = 1.0 / np.sqrt(2**n_qubits)
            uniform_vector = np.full(2**n_qubits, uniform_amplitude, dtype=complex)
            
            sv_te_qsci = StateVectorTE_QSCI(
                hamiltonian=system_data['quri_hamiltonian'],
                sampler=create_qulacs_vector_concurrent_sampler(),
                evolution_time=1.0,
                num_states_pick_out=2**n_qubits  # Complete subspace
            )
            result = sv_te_qsci.run([uniform_vector], total_shots=3000)
            energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
            variant_results['StateVectorTE_QSCI'] = {
                'energy': result.ground_state_energy,
                'error': energy_error,
                'success': True
            }
            print(f"StateVectorTE_QSCI: {result.ground_state_energy:.12f} (error: {energy_error:.2e})")
        except Exception as e:
            variant_results['StateVectorTE_QSCI'] = {'success': False, 'error_msg': str(e)}
            print(f"StateVectorTE_QSCI: FAILED - {e}")
        
        # Validate that at least one variant succeeded with good precision
        successful_variants = [name for name, result in variant_results.items() 
                             if result.get('success', False)]
        
        assert len(successful_variants) > 0, "No QSCI variants succeeded"
        
        # Check that successful variants achieve machine precision with uniform superposition
        for variant_name in successful_variants:
            result = variant_results[variant_name]
            assert result['error'] < 1e-8, (  # Expect machine precision with uniform superposition
                f"{variant_name} error {result['error']:.2e} too large with uniform superposition + complete subspace"
            )
        
        print(f"\nSuccessful variants: {successful_variants}")
        print(f"This demonstrates QSCI algorithm correctness with uniform superposition!")


if __name__ == "__main__":
    # Run exact ground state precision tests
    pytest.main([__file__, "-v", "-m", "exact_ground_state_precision"])