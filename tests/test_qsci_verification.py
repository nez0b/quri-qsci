"""
QSCI Verification Tests using Exact Diagonalization

This module contains verification tests that compare QSCI results against
exact diagonalization for small, well-defined quantum systems. These tests
provide mathematical validation and serve as ground truth for the QSCI
implementation.

The tests use the hamiltonian_factory utilities to create various quantum
model Hamiltonians and validate QSCI results against scipy.sparse.linalg.eigsh.
"""

import pytest
import numpy as np
import sys
import os
from typing import List, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import test utilities
from .utils.hamiltonian_factory import (
    create_hermitian_matrix,
    create_simple_2x2_hamiltonian,
    create_degenerate_hamiltonian,
    create_identity_hamiltonian,
    create_diagonal_hamiltonian,
    create_pauli_hamiltonian,
    create_tfim_hamiltonian,
    create_heisenberg_hamiltonian,
    create_random_sparse_hamiltonian,
    get_exact_solution,
    validate_sparse_hamiltonian_properties
)

try:
    from quri_qsci.qsci_algorithms import VanillaQSCI
    from quri_parts.core.operator import Operator
    from quri_parts.core.operator.pauli import PauliLabel, pauli_label
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: QSCI implementation not available: {e}")
    QSCI_AVAILABLE = False

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh


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
        
    Note:
        This is a simplified conversion that may not be efficient
        for large systems. For production use, consider using
        dedicated Pauli decomposition algorithms.
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


def create_diverse_superposition_state(n_qubits: int) -> GeneralCircuitQuantumState:
    """
    Create a quantum state with diverse measurement outcomes for QSCI testing.

    This creates a superposition state that will yield multiple different
    computational basis states when measured, ensuring QSCI has diverse
    input data to work with.

    Args:
        n_qubits: Number of qubits

    Returns:
        Quantum state with measurement diversity
    """
    circuit = QuantumCircuit(n_qubits)

    # Create a uniform superposition by applying a Hadamard gate to each qubit
    for i in range(n_qubits):
        circuit.add_H_gate(i)

    return GeneralCircuitQuantumState(n_qubits, circuit)


@pytest.mark.verification
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="QSCI implementation not available")
class TestExactDiagonalizationVerification:
    """
    Verification tests comparing QSCI results against exact diagonalization.
    
    These tests provide mathematical validation by comparing QSCI ground state
    energies and eigenvalues against exact solutions computed using
    scipy.sparse.linalg.eigsh.
    """
    
    def test_simple_2x2_hamiltonian_exact(self, tolerance_settings):
        """Test QSCI against exact diagonalization for simple 2x2 case."""
        # Create simple test Hamiltonian
        hamiltonian_matrix, exact_eigenvalues, exact_eigenvectors = create_simple_2x2_hamiltonian()
        
        # Convert to QURI operator
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=1)
        
        # Create diverse state
        diverse_state = create_diverse_superposition_state(n_qubits=1)
        
        # Run QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=2
        )
        
        result = qsci.run([diverse_state], total_shots=1000)
        
        # Validate ground state energy
        exact_ground_energy = np.min(exact_eigenvalues)
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        assert energy_error < tolerance_settings['eigenvalue_atol'], (
            f"Ground state energy error {energy_error:.6e} exceeds tolerance "
            f"{tolerance_settings['eigenvalue_atol']:.6e}. "
            f"QSCI: {result.ground_state_energy:.6f}, "
            f"Exact: {exact_ground_energy:.6f}"
        )
    
    def test_pauli_hamiltonian_exact(self, tolerance_settings):
        """Test QSCI against exact diagonalization for 2-qubit Pauli Hamiltonian."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        
        # Get exact solution
        sparse_ham = csc_matrix(hamiltonian_matrix)
        exact_eigenvalues, exact_eigenvectors = get_exact_solution(sparse_ham, num_states=3)
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Create diverse state
        diverse_state = create_diverse_superposition_state(n_qubits=2)
        
        # Run QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=4  # Use full subspace
        )
        
        result = qsci.run([diverse_state], total_shots=1500)
        
        # Validate ground state energy
        exact_ground_energy = exact_eigenvalues[0]
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        assert energy_error < tolerance_settings['eigenvalue_atol'], (
            f"2-qubit Pauli Hamiltonian ground state energy error {energy_error:.6e} "
            f"exceeds tolerance {tolerance_settings['eigenvalue_atol']:.6e}"
        )
        
        # Validate multiple eigenvalues if available
        if len(result.eigenvalues) > 1:
            for i, qsci_eigenval in enumerate(result.eigenvalues[:len(exact_eigenvalues)]):
                exact_eigenval = exact_eigenvalues[i]
                eigenval_error = abs(qsci_eigenval - exact_eigenval)
                
                assert eigenval_error < tolerance_settings['eigenvalue_atol'], (
                    f"Eigenvalue {i} error {eigenval_error:.6e} exceeds tolerance. "
                    f"QSCI: {qsci_eigenval:.6f}, Exact: {exact_eigenval:.6f}"
                )
    
    @pytest.mark.parametrize("n_qubits,h_field", [
        (2, 0.5),  # Ferromagnetic phase
        (2, 1.0),  # Critical point
        (2, 2.0),  # Paramagnetic phase
        (3, 1.0),  # 3-qubit critical
    ])
    def test_tfim_hamiltonian_exact(self, n_qubits, h_field, tolerance_settings):
        """Test QSCI against exact diagonalization for TFIM."""
        # Create TFIM Hamiltonian
        hamiltonian_matrix = create_tfim_hamiltonian(n_qubits, h_field)

        # Validate Hamiltonian properties
        ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
        assert ham_props['is_valid_hamiltonian'], "TFIM Hamiltonian should be valid"

        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(
            hamiltonian_matrix, num_states=min(3, 2**n_qubits)
        )
        exact_ground_energy = exact_eigenvalues[0]

        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)

        # Create uniform superposition state for complete subspace sampling
        diverse_state = create_diverse_superposition_state(n_qubits)

        # Run QSCI with full subspace selection
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=2**n_qubits  # Use full subspace
        )

        result = qsci.run([diverse_state], total_shots=2000)

        # Validate ground state energy with high precision
        energy_error = abs(result.ground_state_energy - exact_ground_energy)

        # Use a high precision tolerance for validation
        precision_tolerance = 1e-9
        assert energy_error < precision_tolerance, (
            f"TFIM ({n_qubits} qubits, h={h_field}) ground state energy error "
            f"{energy_error:.6e} exceeds precision tolerance {precision_tolerance:.6e}. "
            f"QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
        )
    
    @pytest.mark.parametrize("n_qubits,jz,jxy", [
        (2, 1.0, 1.0),  # Isotropic Heisenberg
        (2, 2.0, 1.0),  # Easy-axis anisotropic
        (2, 0.0, 1.0),  # XX model
        (3, 1.0, 1.0),  # 3-qubit isotropic
    ])
    def test_heisenberg_hamiltonian_exact(self, n_qubits, jz, jxy, tolerance_settings):
        """Test QSCI against exact diagonalization for Heisenberg model."""
        # Create Heisenberg Hamiltonian
        hamiltonian_matrix = create_heisenberg_hamiltonian(n_qubits, jz, jxy)

        # Validate Hamiltonian properties
        ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
        assert ham_props['is_valid_hamiltonian'], "Heisenberg Hamiltonian should be valid"

        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(
            hamiltonian_matrix, num_states=min(3, 2**n_qubits)
        )
        exact_ground_energy = exact_eigenvalues[0]

        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)

        # Create uniform superposition state for complete subspace sampling
        diverse_state = create_diverse_superposition_state(n_qubits)

        # Run QSCI with full subspace selection
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=2**n_qubits  # Use full subspace
        )

        result = qsci.run([diverse_state], total_shots=2000)

        # Validate ground state energy with high precision
        energy_error = abs(result.ground_state_energy - exact_ground_energy)

        # Use a high precision tolerance for validation
        precision_tolerance = 1e-9
        assert energy_error < precision_tolerance, (
            f"Heisenberg ({n_qubits} qubits, Jz={jz}, Jxy={jxy}) ground state energy error "
            f"{energy_error:.6e} exceeds precision tolerance {precision_tolerance:.6e}. "
            f"QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
        )
    
    def test_degenerate_eigenvalue_case(self, tolerance_settings):
        """Test QSCI with degenerate eigenvalues (edge case)."""
        # Create Hamiltonian with known degeneracies
        hamiltonian_matrix, known_eigenvalues = create_degenerate_hamiltonian(size=4)

        # This is a 4x4 matrix, so we can treat it as a 2-qubit system
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)

        # Create uniform superposition state for complete subspace sampling
        diverse_state = create_diverse_superposition_state(n_qubits=2)

        # Run QSCI with full subspace selection
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=4  # Use full subspace
        )

        result = qsci.run([diverse_state], total_shots=1500)

        # Validate ground state energy with high precision
        exact_ground_energy = np.min(known_eigenvalues)
        energy_error = abs(result.ground_state_energy - exact_ground_energy)

        # Use a high precision tolerance for validation
        precision_tolerance = 1e-9
        assert energy_error < precision_tolerance, (
            f"Degenerate case ground state energy error {energy_error:.6e} "
            f"exceeds precision tolerance {precision_tolerance:.6e}"
        )
    
    def test_random_sparse_hamiltonian_stress(self, tolerance_settings):
        """Stress test QSCI with random sparse Hamiltonians."""
        for trial in range(3):  # Multiple random trials
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
            exact_ground_energy = exact_eigenvalues[0]

            # Convert to QURI operator
            quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits)

            # Create uniform superposition state for complete subspace sampling
            diverse_state = create_diverse_superposition_state(n_qubits)

            # Run QSCI with full subspace selection
            sampler = create_qulacs_vector_concurrent_sampler()
            qsci = VanillaQSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                num_states_pick_out=2**n_qubits  # Use full subspace
            )

            result = qsci.run([diverse_state], total_shots=1500)

            # Validate ground state energy with high precision
            energy_error = abs(result.ground_state_energy - exact_ground_energy)

            # Use a high precision tolerance for validation
            precision_tolerance = 1e-9
            assert energy_error < precision_tolerance, (
                f"Random sparse Hamiltonian trial {trial} ground state energy error "
                f"{energy_error:.6e} exceeds precision tolerance {precision_tolerance:.6e}. "
                f"QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
            )


@pytest.mark.verification
class TestHamiltonianFactory:
    """Test the Hamiltonian factory functions themselves."""
    
    def test_tfim_hamiltonian_properties(self):
        """Test TFIM Hamiltonian creation and properties."""
        hamiltonian = create_tfim_hamiltonian(n_qubits=3, h=1.0)
        props = validate_sparse_hamiltonian_properties(hamiltonian)
        
        assert props['is_valid_hamiltonian'], "TFIM should be a valid Hamiltonian"
        assert props['is_sparse'], "TFIM should be sparse"
        assert hamiltonian.shape == (8, 8), "3-qubit TFIM should be 8x8"
    
    def test_heisenberg_hamiltonian_properties(self):
        """Test Heisenberg Hamiltonian creation and properties."""
        hamiltonian = create_heisenberg_hamiltonian(n_qubits=2, jz=1.0, jxy=1.0)
        props = validate_sparse_hamiltonian_properties(hamiltonian)
        
        assert props['is_valid_hamiltonian'], "Heisenberg should be a valid Hamiltonian"
        assert props['is_sparse'], "Heisenberg should be sparse"
        assert hamiltonian.shape == (4, 4), "2-qubit Heisenberg should be 4x4"
    
    def test_random_sparse_hamiltonian_properties(self):
        """Test random sparse Hamiltonian creation and properties."""
        hamiltonian = create_random_sparse_hamiltonian(n_qubits=2, density=0.2)
        props = validate_sparse_hamiltonian_properties(hamiltonian)
        
        assert props['is_valid_hamiltonian'], "Random sparse should be a valid Hamiltonian"
        assert props['is_sparse'], "Random sparse should be sparse"
        # With density=0.2 and Hermitian symmetrization, expect roughly 40-60% sparsity
        assert 0.3 < props['sparsity'] < 0.7, f"Expected 30-70% sparsity, got {props['sparsity']:.2f}"
    
    def test_exact_solution_function(self):
        """Test the exact solution computation function."""
        # Simple 2x2 case
        matrix = create_simple_2x2_hamiltonian()[0]
        sparse_matrix = csc_matrix(matrix)
        
        eigenvalues, eigenvectors = get_exact_solution(sparse_matrix, num_states=2)
        
        assert len(eigenvalues) == 2, "Should return 2 eigenvalues"
        assert eigenvectors.shape == (2, 2), "Should return 2x2 eigenvector matrix"
        assert eigenvalues[0] <= eigenvalues[1], "Eigenvalues should be sorted"


@pytest.mark.verification 
class TestEdgeCasesVerification:
    """Test edge cases in verification scenarios."""
    
    def test_identity_hamiltonian_trivial_case(self):
        """Test QSCI with identity Hamiltonian (trivial case)."""
        if not QSCI_AVAILABLE:
            pytest.skip("QSCI not available")
            
        # Identity Hamiltonian - all eigenvalues should be 1
        identity_matrix = create_identity_hamiltonian(4)  # 2-qubit
        sparse_ham = csc_matrix(identity_matrix)
        
        # Exact solution
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=2)
        assert np.allclose(exact_eigenvalues, 1.0), "Identity eigenvalues should all be 1"
        
        # Convert to QURI operator (should be just identity term)
        quri_hamiltonian = Operator()
        quri_hamiltonian[PauliLabel()] = 1.0  # Identity operator
        
        # Create state
        diverse_state = create_diverse_superposition_state(n_qubits=2)
        
        # Run QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=2
        )
        
        result = qsci.run([diverse_state], total_shots=1000)
        
        # Ground state energy should be close to 1
        energy_error = abs(result.ground_state_energy - 1.0)
        assert energy_error < 1e-6, f"Identity case energy error: {energy_error}"
    
    def test_diagonal_hamiltonian_known_eigenvalues(self):
        """Test QSCI with diagonal Hamiltonian (known eigenvalues)."""
        if not QSCI_AVAILABLE:
            pytest.skip("QSCI not available")
            
        # Diagonal Hamiltonian with known eigenvalues
        known_eigenvalues = np.array([-2.0, -1.0, 0.5, 1.5])
        diagonal_matrix = create_diagonal_hamiltonian(known_eigenvalues)
        sparse_ham = csc_matrix(diagonal_matrix)
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Create diverse state
        diverse_state = create_diverse_superposition_state(n_qubits=2)
        
        # Run QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        qsci = VanillaQSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            num_states_pick_out=4  # Use full subspace
        )
        
        result = qsci.run([diverse_state], total_shots=1500)
        
        # Ground state should be minimum eigenvalue
        exact_ground = np.min(known_eigenvalues)
        energy_error = abs(result.ground_state_energy - exact_ground)
        assert energy_error < 1e-6, (
            f"Diagonal case ground state error: {energy_error}. "
            f"QSCI: {result.ground_state_energy}, Expected: {exact_ground}"
        )


if __name__ == "__main__":
    # Run verification tests
    pytest.main([__file__, "-v", "-m", "verification"])