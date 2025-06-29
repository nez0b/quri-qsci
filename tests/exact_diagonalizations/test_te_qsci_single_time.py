"""
Single-Time TE-QSCI Verification Tests using Exact Diagonalization

This module contains verification tests that compare single-time TE-QSCI results 
against exact diagonalization for small, well-defined quantum systems. These tests
validate the single-time evolution variant of TE-QSCI algorithms.
"""

import pytest
import numpy as np
import sys
import os
from typing import List, Tuple

# Add parent directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

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

try:
    from src.qsci_algorithms import SingleTimeTE_QSCI
    from quri_parts.core.operator import Operator
    from quri_parts.core.operator.pauli import PauliLabel, pauli_label
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Single-time TE-QSCI implementation not available: {e}")
    QSCI_AVAILABLE = False

from scipy.sparse import csc_matrix


def sparse_matrix_to_quri_operator(sparse_matrix: csc_matrix, 
                                   n_qubits: int) -> Operator:
    """
    Convert a sparse Hamiltonian matrix to QURI-Parts Operator format.
    
    This function converts a 2^n x 2^n sparse matrix representation
    of a quantum Hamiltonian into the Pauli operator format used
    by QURI-Parts.
    
    Args:
        sparse_matrix: Sparse matrix representation of Hamiltonian
        n_qubits: Number of qubits
        
    Returns:
        QURI-Parts Operator object
    """
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


def generate_pauli_strings(n_qubits: int) -> List[Tuple[PauliLabel, np.ndarray]]:
    """
    Generate all Pauli strings and their matrix representations for n qubits.
    
    Args:
        n_qubits: Number of qubits
        
    Returns:
        List of (PauliLabel, matrix) tuples
    """
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


@pytest.mark.te_qsci_single_time
@pytest.mark.exact_diagonalization
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="Single-time TE-QSCI implementation not available")
class TestSingleTimeTE_QSCI_ExactVerification:
    """
    Verification tests comparing single-time TE-QSCI results against exact diagonalization.
    
    These tests provide mathematical validation by comparing single-time TE-QSCI ground state
    energies and eigenvalues against exact solutions computed using
    scipy.sparse.linalg.eigsh.
    """
    
    def test_simple_2x2_hamiltonian_single_time_te_qsci_exact(self, tolerance_settings, evolution_time_parameters):
        """Test single-time TE-QSCI against exact diagonalization for simple 2x2 case."""
        # Create simple test Hamiltonian
        hamiltonian_matrix, exact_eigenvalues, exact_eigenvectors = create_simple_2x2_hamiltonian()
        
        # Convert to QURI operator
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=1)
        
        # Test multiple evolution times
        for evolution_time in evolution_time_parameters[:2]:  # Test first 2 evolution times
            # Create diverse state
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits=1, theta=np.pi/4)
            
            # Run single-time TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2
            )
            
            result = te_qsci.run([diverse_state], total_shots=1000)
            
            # Validate ground state energy
            exact_ground_energy = np.min(exact_eigenvalues)
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            assert energy_error < tolerance_settings['qsci_verification_atol'], (
                f"Single-time TE-QSCI (t={evolution_time}) ground state energy error "
                f"{energy_error:.6e} exceeds tolerance {tolerance_settings['qsci_verification_atol']:.6e}. "
                f"TE-QSCI: {result.ground_state_energy:.6f}, "
                f"Exact: {exact_ground_energy:.6f}"
            )
    
    def test_pauli_hamiltonian_single_time_te_qsci_exact(self, tolerance_settings, evolution_time_parameters):
        """Test single-time TE-QSCI against exact diagonalization for 2-qubit Pauli Hamiltonian."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        
        # Get exact solution
        sparse_ham = csc_matrix(hamiltonian_matrix)
        exact_eigenvalues, exact_eigenvectors = get_exact_solution(sparse_ham, num_states=3)
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Test multiple evolution times
        for evolution_time in evolution_time_parameters:
            # Create uniform superposition state (already updated in conftest.py)
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/6)
            
            # Run single-time TE-QSCI with complete subspace coverage
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=4  # Complete subspace for 2-qubit system
            )
            
            result = te_qsci.run([diverse_state], total_shots=2000)  # Higher shots for better convergence
            
            # Validate ground state energy
            exact_ground_energy = exact_eigenvalues[0]
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"2-qubit Pauli Hamiltonian single-time TE-QSCI (t={evolution_time}) "
                f"ground state energy error {energy_error:.6e} "
                f"exceeds tolerance {precision_tolerance:.6e}"
            )
    
    @pytest.mark.parametrize("n_qubits,h_field", [
        (2, 0.5),  # Ferromagnetic phase
        (2, 1.0),  # Critical point
        (2, 2.0),  # Paramagnetic phase
    ])
    def test_tfim_hamiltonian_single_time_te_qsci_exact(self, n_qubits, h_field, tolerance_settings, evolution_time_parameters):
        """Test single-time TE-QSCI against exact diagonalization for TFIM."""
        # Create TFIM Hamiltonian
        hamiltonian_matrix = create_tfim_hamiltonian(n_qubits, h_field)
        
        # Validate Hamiltonian properties
        ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
        assert ham_props['is_valid_hamiltonian'], "TFIM Hamiltonian should be valid"
        
        # Get exact solution
        exact_eigenvalues, exact_eigenvectors = get_exact_solution(
            hamiltonian_matrix, num_states=min(3, 2**n_qubits)
        )
        
        # Convert to QURI operator  
        quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)
        
        # Test evolution time parameter sweep
        for evolution_time in evolution_time_parameters[:3]:  # Test first 3 evolution times
            # Create uniform superposition state (already updated in conftest.py)
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/8)
            
            # Run single-time TE-QSCI with complete subspace coverage
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**n_qubits  # Complete subspace coverage
            )
            
            result = te_qsci.run([diverse_state], total_shots=3000)  # Higher shots for TE-QSCI
            
            # Validate ground state energy
            exact_ground_energy = exact_eigenvalues[0]
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"TFIM ({n_qubits} qubits, h={h_field}) single-time TE-QSCI (t={evolution_time}) "
                f"ground state energy error {energy_error:.6e} exceeds tolerance "
                f"{precision_tolerance:.6e}. TE-QSCI: {result.ground_state_energy:.6f}, "
                f"Exact: {exact_ground_energy:.6f}"
            )
    
    @pytest.mark.parametrize("n_qubits,jz,jxy", [
        (2, 1.0, 1.0),  # Isotropic Heisenberg
        (2, 2.0, 1.0),  # Easy-axis anisotropic
        (2, 0.0, 1.0),  # XX model
    ])
    def test_heisenberg_hamiltonian_single_time_te_qsci_exact(self, n_qubits, jz, jxy, tolerance_settings, evolution_time_parameters):
        """Test single-time TE-QSCI against exact diagonalization for Heisenberg model."""
        # Create Heisenberg Hamiltonian
        hamiltonian_matrix = create_heisenberg_hamiltonian(n_qubits, jz, jxy)
        
        # Validate Hamiltonian properties
        ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
        assert ham_props['is_valid_hamiltonian'], "Heisenberg Hamiltonian should be valid"
        
        # Get exact solution
        exact_eigenvalues, exact_eigenvectors = get_exact_solution(
            hamiltonian_matrix, num_states=min(3, 2**n_qubits)
        )
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)
        
        # Test evolution time parameter sweep
        for evolution_time in evolution_time_parameters[:2]:  # Test first 2 evolution times
            # Create uniform superposition state (already updated in conftest.py)
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/5)
            
            # Run single-time TE-QSCI with complete subspace coverage
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**n_qubits  # Complete subspace coverage
            )
            
            result = te_qsci.run([diverse_state], total_shots=3000)  # Higher shots for TE-QSCI
            
            # Validate ground state energy
            exact_ground_energy = exact_eigenvalues[0]
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"Heisenberg ({n_qubits} qubits, Jz={jz}, Jxy={jxy}) single-time TE-QSCI "
                f"(t={evolution_time}) ground state energy error {energy_error:.6e} exceeds "
                f"tolerance {precision_tolerance:.6e}. "
                f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
            )
    
    def test_evolution_time_parameter_sweep(self, tolerance_settings, evolution_time_parameters):
        """Test evolution time dependency for single-time TE-QSCI."""
        # Use simple 2-qubit Pauli Hamiltonian for parameter sweep
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact ground state energy
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Create uniform superposition state (already updated in conftest.py)
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/4)
        
        energy_errors = []
        
        # Test all evolution times
        for evolution_time in evolution_time_parameters:
            # Run single-time TE-QSCI with complete subspace coverage
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=4  # Complete subspace for 2-qubit system
            )
            
            result = te_qsci.run([diverse_state], total_shots=3000)  # Higher shots for TE-QSCI
            
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            energy_errors.append(energy_error)
            
            # Expect high precision with uniform superposition + complete subspace
            precision_tolerance = 1e-8
            assert energy_error < precision_tolerance, (
                f"Evolution time sweep: t={evolution_time} energy error {energy_error:.6e} "
                f"exceeds tolerance {precision_tolerance:.6e}"
            )
        
        # Validate that we get consistent results across evolution times
        # (For a stable algorithm, energy errors shouldn't vary too wildly)
        energy_error_std = np.std(energy_errors)
        assert energy_error_std < 0.3, (
            f"Energy error standard deviation {energy_error_std:.6e} too large across "
            f"evolution times. Errors: {energy_errors}"
        )
    
    def test_random_sparse_hamiltonian_single_time_te_qsci_stress(self, tolerance_settings, evolution_time_parameters):
        """Stress test single-time TE-QSCI with random sparse Hamiltonians."""
        for trial in range(2):  # Multiple random trials
            # Create random sparse Hamiltonian
            n_qubits = 2  # Keep small for exact diagonalization
            hamiltonian_matrix = create_random_sparse_hamiltonian(
                n_qubits, density=0.3, random_seed=42 + trial
            )
            
            # Validate Hamiltonian properties
            ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
            assert ham_props['is_valid_hamiltonian'], f"Random Hamiltonian {trial} should be valid"
            
            # Get exact solution
            exact_eigenvalues, _ = get_exact_solution(hamiltonian_matrix, num_states=2)
            
            # Convert to QURI operator
            quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)
            
            # Test with one evolution time
            evolution_time = evolution_time_parameters[1]  # Use middle evolution time
            
            # Create diverse state
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/7)
            
            # Run single-time TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = SingleTimeTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=3
            )
            
            result = te_qsci.run([diverse_state], total_shots=1500)
            
            # Validate ground state energy
            exact_ground_energy = exact_eigenvalues[0]
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use loose tolerance for random Hamiltonians
            verification_tolerance = max(tolerance_settings['qsci_loose_atol'], 0.3)
            assert energy_error < verification_tolerance, (
                f"Random sparse Hamiltonian trial {trial} single-time TE-QSCI "
                f"(t={evolution_time}) ground state energy error {energy_error:.6e} "
                f"exceeds verification tolerance {verification_tolerance:.6e}. "
                f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
            )


if __name__ == "__main__":
    # Run single-time TE-QSCI verification tests
    pytest.main([__file__, "-v", "-m", "te_qsci_single_time"])