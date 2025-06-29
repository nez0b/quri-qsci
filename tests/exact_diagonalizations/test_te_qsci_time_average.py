"""
Time-Average TE-QSCI Verification Tests using Exact Diagonalization

This module contains verification tests that compare time-average TE-QSCI results 
against exact diagonalization for small, well-defined quantum systems. These tests
validate the time-averaged evolution variant of TE-QSCI algorithms.
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
    from src.qsci_algorithms import TimeAverageTE_QSCI, SingleTimeTE_QSCI
    from quri_parts.core.operator import Operator
    from quri_parts.core.operator.pauli import PauliLabel, pauli_label
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Time-average TE-QSCI implementation not available: {e}")
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


@pytest.mark.te_qsci_time_average
@pytest.mark.exact_diagonalization
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="Time-average TE-QSCI implementation not available")
class TestTimeAverageTE_QSCI_ExactVerification:
    """
    Verification tests comparing time-average TE-QSCI results against exact diagonalization.
    
    These tests provide mathematical validation by comparing time-average TE-QSCI ground state
    energies and eigenvalues against exact solutions computed using
    scipy.sparse.linalg.eigsh.
    """
    
    def test_simple_2x2_hamiltonian_time_average_te_qsci_exact(self, tolerance_settings, evolution_time_parameters):
        """Test time-average TE-QSCI against exact diagonalization for simple 2x2 case."""
        # Create simple test Hamiltonian
        hamiltonian_matrix, exact_eigenvalues, exact_eigenvectors = create_simple_2x2_hamiltonian()
        
        # Convert to QURI operator
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=1)
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=1, theta=np.pi/4)
        
        # Run time-average TE-QSCI with multiple evolution times
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = TimeAverageTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_times=evolution_time_parameters[:3],  # Use first 3 evolution times
            num_states_pick_out=2**1  # Complete subspace coverage
        )
        
        result = te_qsci.run([diverse_state], total_shots=2000)
        
        # Validate ground state energy
        exact_ground_energy = np.min(exact_eigenvalues)
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        assert energy_error < 1e-8, (
            f"Time-average TE-QSCI ground state energy error "
            f"{energy_error:.6e} exceeds high precision tolerance 1e-8. "
            f"TE-QSCI: {result.ground_state_energy:.6f}, "
            f"Exact: {exact_ground_energy:.6f}"
        )
    
    def test_pauli_hamiltonian_time_average_te_qsci_exact(self, tolerance_settings, evolution_time_parameters):
        """Test time-average TE-QSCI against exact diagonalization for 2-qubit Pauli Hamiltonian."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        
        # Get exact solution
        sparse_ham = csc_matrix(hamiltonian_matrix)
        exact_eigenvalues, exact_eigenvectors = get_exact_solution(sparse_ham, num_states=3)
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/6)
        
        # Run time-average TE-QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = TimeAverageTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_times=evolution_time_parameters,  # Use all evolution times
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        
        result = te_qsci.run([diverse_state], total_shots=2500)
        
        # Validate ground state energy
        exact_ground_energy = exact_eigenvalues[0]
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        assert energy_error < 1e-8, (
            f"2-qubit Pauli Hamiltonian time-average TE-QSCI "
            f"ground state energy error {energy_error:.6e} "
            f"exceeds high precision tolerance 1e-8"
        )
        
        # Validate multiple eigenvalues if available
        if len(result.eigenvalues) > 1:
            for i, qsci_eigenval in enumerate(result.eigenvalues[:len(exact_eigenvalues)]):
                exact_eigenval = exact_eigenvalues[i]
                eigenval_error = abs(qsci_eigenval - exact_eigenval)
                
                assert eigenval_error < 1e-8, (
                    f"Eigenvalue {i} error {eigenval_error:.6e} exceeds high precision tolerance 1e-8. "
                    f"TE-QSCI: {qsci_eigenval:.6f}, Exact: {exact_eigenval:.6f}"
                )
    
    @pytest.mark.parametrize("n_qubits,h_field", [
        (2, 0.5),  # Ferromagnetic phase
        (2, 1.0),  # Critical point
        (2, 2.0),  # Paramagnetic phase
    ])
    def test_tfim_hamiltonian_time_average_te_qsci_exact(self, n_qubits, h_field, tolerance_settings, evolution_time_parameters):
        """Test time-average TE-QSCI against exact diagonalization for TFIM."""
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
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/8)
        
        # Run time-average TE-QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = TimeAverageTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_times=evolution_time_parameters[:3],  # Use first 3 evolution times
            num_states_pick_out=2**n_qubits  # Complete subspace coverage
        )
        
        result = te_qsci.run([diverse_state], total_shots=3000)
        
        # Validate ground state energy
        exact_ground_energy = exact_eigenvalues[0]
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        # Use high precision tolerance for uniform superposition + complete subspace
        assert energy_error < 1e-8, (
            f"TFIM ({n_qubits} qubits, h={h_field}) time-average TE-QSCI "
            f"ground state energy error {energy_error:.6e} exceeds high precision tolerance 1e-8. "
            f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
        )
    
    @pytest.mark.parametrize("n_qubits,jz,jxy", [
        (2, 1.0, 1.0),  # Isotropic Heisenberg
        (2, 2.0, 1.0),  # Easy-axis anisotropic
    ])
    def test_heisenberg_hamiltonian_time_average_te_qsci_exact(self, n_qubits, jz, jxy, tolerance_settings, evolution_time_parameters):
        """Test time-average TE-QSCI against exact diagonalization for Heisenberg model."""
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
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/5)
        
        # Run time-average TE-QSCI
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = TimeAverageTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_times=evolution_time_parameters[:3],  # Use first 3 evolution times
            num_states_pick_out=2**n_qubits  # Complete subspace coverage
        )
        
        result = te_qsci.run([diverse_state], total_shots=3000)
        
        # Validate ground state energy
        exact_ground_energy = exact_eigenvalues[0]
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        # Use high precision tolerance for uniform superposition + complete subspace
        assert energy_error < 1e-8, (
            f"Heisenberg ({n_qubits} qubits, Jz={jz}, Jxy={jxy}) time-average TE-QSCI "
            f"ground state energy error {energy_error:.6e} exceeds high precision tolerance 1e-8. "
            f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
        )
    
    def test_time_averaging_convergence(self, tolerance_settings, evolution_time_parameters):
        """Test convergence properties of time-averaging vs number of time points."""
        # Use simple 2-qubit Pauli Hamiltonian for convergence study
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact ground state energy
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/4)
        
        # Test different numbers of time points for averaging
        time_point_counts = [2, 3, 4]  # Different numbers of evolution times
        energy_errors = []
        
        for num_time_points in time_point_counts:
            evolution_times = evolution_time_parameters[:num_time_points]
            
            # Run time-average TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = TimeAverageTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_times=evolution_times,
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([diverse_state], total_shots=2500)
            
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            energy_errors.append(energy_error)
            
            # Individual validation with high precision
            assert energy_error < 1e-8, (
                f"Time averaging convergence: {num_time_points} time points, "
                f"energy error {energy_error:.6e} exceeds high precision tolerance 1e-8"
            )
        
        # Validate convergence trend (more time points should generally give better or similar results)
        # We're not enforcing strict monotonicity due to statistical noise, but check reasonableness
        max_error = max(energy_errors)
        min_error = min(energy_errors)
        # Handle case where min_error might be 0 (perfect convergence)
        min_error_threshold = max(min_error, 1e-15)  # Use machine precision as minimum
        assert max_error < 2 * min_error_threshold, (
            f"Time averaging convergence: error variation too large across time point counts. "
            f"Errors: {energy_errors}, max: {max_error:.2e}, min: {min_error:.2e}"
        )
    
    def test_statistical_variance_analysis(self, tolerance_settings, evolution_time_parameters):
        """Test statistical properties of time-average TE-QSCI measurements."""
        # Use simple 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/6)
        
        # Run multiple instances of time-average TE-QSCI for variance analysis
        results = []
        for trial in range(3):  # Multiple trials for statistical analysis
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = TimeAverageTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_times=evolution_time_parameters[:3],
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([diverse_state], total_shots=2500)
            results.append(result.ground_state_energy)
        
        # Analyze statistical properties
        mean_energy = np.mean(results)
        std_energy = np.std(results)
        
        # Validate that standard deviation is reasonable (not too large)
        # For a robust algorithm, multiple runs should give consistent results
        assert std_energy < 0.2, (
            f"Statistical variance too large: std={std_energy:.6e} across {len(results)} trials. "
            f"Results: {results}"
        )
        
        # Validate that all individual results are within reasonable bounds
        for i, energy in enumerate(results):
            deviation = abs(energy - mean_energy)
            assert deviation < 3 * std_energy, (
                f"Trial {i} energy {energy:.6f} deviates too much from mean {mean_energy:.6f}"
            )
    
    def test_comparison_with_single_time_variants(self, tolerance_settings, evolution_time_parameters):
        """Compare time-average TE-QSCI with single-time TE-QSCI variants."""
        # Use simple 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact ground state energy
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Create diverse state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        diverse_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/4)
        
        # Run time-average TE-QSCI
        sampler_avg = create_qulacs_vector_concurrent_sampler()
        te_qsci_avg = TimeAverageTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler_avg,
            evolution_times=evolution_time_parameters[:3],
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        result_avg = te_qsci_avg.run([diverse_state], total_shots=2500)
        
        # Run single-time TE-QSCI for comparison (use middle evolution time)
        evolution_time = evolution_time_parameters[1]
        sampler_single = create_qulacs_vector_concurrent_sampler()
        te_qsci_single = SingleTimeTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler_single,
            evolution_time=evolution_time,
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        result_single = te_qsci_single.run([diverse_state], total_shots=2500)
        
        # Validate both methods
        energy_error_avg = abs(result_avg.ground_state_energy - exact_ground_energy)
        energy_error_single = abs(result_single.ground_state_energy - exact_ground_energy)
        
        # Use high precision tolerance for uniform superposition + complete subspace
        assert energy_error_avg < 1e-8, (
            f"Time-average TE-QSCI energy error {energy_error_avg:.6e} exceeds high precision tolerance 1e-8"
        )
        assert energy_error_single < 1e-8, (
            f"Single-time TE-QSCI energy error {energy_error_single:.6e} exceeds high precision tolerance 1e-8"
        )
        
        # Compare performance (both should achieve machine precision)
        avg_vs_single_diff = abs(result_avg.ground_state_energy - result_single.ground_state_energy)
        assert avg_vs_single_diff < 1e-8, (
            f"Time-average vs single-time TE-QSCI difference too large: "
            f"{avg_vs_single_diff:.6e}. Avg: {result_avg.ground_state_energy:.6f}, "
            f"Single: {result_single.ground_state_energy:.6f}"
        )
    
    def test_random_sparse_hamiltonian_time_average_te_qsci_stress(self, tolerance_settings, evolution_time_parameters):
        """Stress test time-average TE-QSCI with random sparse Hamiltonians."""
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
            
            # Create diverse state
            from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
            diverse_state = create_diverse_superposition_state(n_qubits, theta=np.pi/7)
            
            # Run time-average TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = TimeAverageTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_times=evolution_time_parameters[:3],  # Use first 3 evolution times
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([diverse_state], total_shots=2500)
            
            # Validate ground state energy
            exact_ground_energy = exact_eigenvalues[0]
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use high precision tolerance for uniform superposition + complete subspace
            assert energy_error < 1e-8, (
                f"Random sparse Hamiltonian trial {trial} time-average TE-QSCI "
                f"ground state energy error {energy_error:.6e} "
                f"exceeds high precision tolerance 1e-8. "
                f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
            )


if __name__ == "__main__":
    # Run time-average TE-QSCI verification tests
    pytest.main([__file__, "-v", "-m", "te_qsci_time_average"])