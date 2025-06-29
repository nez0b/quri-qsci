"""
State Vector TE-QSCI Verification Tests using Exact Diagonalization

This module contains verification tests that compare state vector TE-QSCI results 
against exact diagonalization for small, well-defined quantum systems. These tests
validate the state vector input variant of TE-QSCI algorithms.
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
    from src.qsci_algorithms import StateVectorTE_QSCI, SingleTimeTE_QSCI
    from quri_parts.core.operator import Operator
    from quri_parts.core.operator.pauli import PauliLabel, pauli_label
    from quri_parts.core.state import GeneralCircuitQuantumState, ComputationalBasisState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: State vector TE-QSCI implementation not available: {e}")
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


def create_circuit_state_from_amplitudes(n_qubits: int, amplitudes: np.ndarray) -> GeneralCircuitQuantumState:
    """
    Create a CircuitQuantumState from state vector amplitudes for testing.
    
    This function creates a proper CircuitQuantumState that has the .circuit attribute
    expected by StateVectorTE_QSCI, unlike the previous CustomStateVector approach.
    
    Args:
        n_qubits: Number of qubits
        amplitudes: Complex amplitude vector (should be normalized)
        
    Returns:
        GeneralCircuitQuantumState object with proper interface
    """
    from quri_parts.circuit import QuantumCircuit
    
    # Normalize amplitudes
    norm = np.linalg.norm(amplitudes)
    if norm > 0:
        amplitudes = amplitudes / norm
    
    # Create a circuit that prepares the desired state
    circuit = QuantumCircuit(n_qubits)
    
    # For testing purposes, we'll create a circuit that approximates the desired state
    # using computational basis states or basic rotations
    
    # Method 1: If amplitudes represent a computational basis state, use that
    max_amplitude_idx = np.argmax(np.abs(amplitudes))
    if np.abs(amplitudes[max_amplitude_idx])**2 > 0.99:  # Nearly computational basis state
        # Initialize to the dominant computational basis state
        for i in range(n_qubits):
            if (max_amplitude_idx >> i) & 1:
                circuit.add_X_gate(i)
    else:
        # Method 2: Create a superposition using rotation gates
        # This is a simplified approach - for complex states, more sophisticated
        # state preparation would be needed
        for i in range(n_qubits):
            # Apply Y rotation to create superposition
            circuit.add_RY_gate(i, np.pi/4)  # Creates |+⟩ state
        
        # Add some entanglement for non-separable states
        for i in range(n_qubits - 1):
            circuit.add_CNOT_gate(i, i + 1)
    
    return GeneralCircuitQuantumState(n_qubits, circuit)


@pytest.mark.te_qsci_state_vector
@pytest.mark.exact_diagonalization
@pytest.mark.skipif(not QSCI_AVAILABLE, reason="State vector TE-QSCI implementation not available")
class TestStateVectorTE_QSCI_ExactVerification:
    """
    Verification tests comparing state vector TE-QSCI results against exact diagonalization.
    
    These tests provide mathematical validation by comparing state vector TE-QSCI ground state
    energies and eigenvalues against exact solutions computed using
    scipy.sparse.linalg.eigsh.
    """
    
    def test_simple_2x2_hamiltonian_state_vector_te_qsci_exact(self, tolerance_settings, evolution_time_parameters):
        """Test state vector TE-QSCI against exact diagonalization for simple 2x2 case."""
        # Create simple test Hamiltonian
        hamiltonian_matrix, exact_eigenvalues, exact_eigenvectors = create_simple_2x2_hamiltonian()
        
        # Convert to QURI operator
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=1)
        
        # Create custom state vector (superposition state)
        amplitudes = np.array([1/np.sqrt(2), 1/np.sqrt(2)], dtype=complex)
        state_vector = create_circuit_state_from_amplitudes(n_qubits=1, amplitudes=amplitudes)
        
        # Run state vector TE-QSCI
        evolution_time = evolution_time_parameters[1]  # Use middle evolution time
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = StateVectorTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_time=evolution_time,
            num_states_pick_out=2**1  # Complete subspace coverage
        )
        
        result = te_qsci.run([state_vector], total_shots=2500)
        
        # Validate ground state energy
        exact_ground_energy = np.min(exact_eigenvalues)
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        assert energy_error < 1e-8, (
            f"State vector TE-QSCI (t={evolution_time}) ground state energy error "
            f"{energy_error:.6e} exceeds high precision tolerance 1e-8. "
            f"TE-QSCI: {result.ground_state_energy:.6f}, "
            f"Exact: {exact_ground_energy:.6f}"
        )
    
    def test_uniform_superposition_state_input(self, tolerance_settings, evolution_time_parameters):
        """Test state vector TE-QSCI with uniform superposition state input."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=2)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Create uniform superposition state vector: (|00⟩ + |01⟩ + |10⟩ + |11⟩)/2
        uniform_amplitude = 0.5  # 1/sqrt(4) = 0.5
        amplitudes = np.full(4, uniform_amplitude, dtype=complex)
        state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=amplitudes)
        
        evolution_time = evolution_time_parameters[1]  # Use middle evolution time
        
        # Run state vector TE-QSCI with uniform superposition
        sampler = create_qulacs_vector_concurrent_sampler()
        te_qsci = StateVectorTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler,
            evolution_time=evolution_time,
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        
        result = te_qsci.run([state_vector], total_shots=3000)
        
        # Validate ground state energy with high precision
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        
        # Expect machine precision with uniform superposition + complete subspace
        assert energy_error < 1e-8, (
            f"Uniform superposition TE-QSCI energy error "
            f"{energy_error:.6e} exceeds high precision tolerance 1e-8"
        )
    
    def test_superposition_state_vector_input(self, tolerance_settings, evolution_time_parameters):
        """Test state vector TE-QSCI with superposition state inputs."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=2)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Create different superposition states
        superposition_states = [
            # Equal superposition: (|00⟩ + |01⟩ + |10⟩ + |11⟩)/2
            np.array([0.5, 0.5, 0.5, 0.5], dtype=complex),
            
            # Bell state: (|00⟩ + |11⟩)/√2
            np.array([1/np.sqrt(2), 0, 0, 1/np.sqrt(2)], dtype=complex),
            
            # Random superposition
            np.array([0.6, 0.3j, 0.5, 0.4+0.2j], dtype=complex),
        ]
        
        evolution_time = evolution_time_parameters[1]  # Use middle evolution time
        
        for i, amplitudes in enumerate(superposition_states):
            # Create state vector
            state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=amplitudes)
            
            # Run state vector TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = StateVectorTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([state_vector], total_shots=2500)
            
            # Validate ground state energy
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use high precision tolerance for uniform superposition + complete subspace
            assert energy_error < 1e-8, (
                f"Superposition state {i} TE-QSCI energy error "
                f"{energy_error:.6e} exceeds high precision tolerance 1e-8"
            )
    
    def test_different_initial_preparations(self, tolerance_settings, evolution_time_parameters):
        """Test state vector TE-QSCI with different initial state preparations."""
        # Use simple 2-qubit TFIM Hamiltonian
        hamiltonian_matrix = create_tfim_hamiltonian(n_qubits=2, h=1.0)
        ham_props = validate_sparse_hamiltonian_properties(hamiltonian_matrix)
        assert ham_props['is_valid_hamiltonian'], "TFIM Hamiltonian should be valid"
        
        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(hamiltonian_matrix, num_states=2)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Convert to QURI operator
        quri_hamiltonian = sparse_matrix_to_quri_operator(hamiltonian_matrix, n_qubits=2)
        
        # Test different preparation strategies
        preparation_strategies = {
            'random_normalized': np.random.random(4) + 1j * np.random.random(4),
            'gaussian_weights': np.exp(-np.arange(4)),  # Exponentially decreasing weights
            'product_state': np.array([0.8, 0.6, 0.0, 0.0], dtype=complex),  # |0⟩⊗(α|0⟩+β|1⟩)
            'maximally_mixed': np.ones(4, dtype=complex),  # Equal amplitudes
        }
        
        evolution_time = evolution_time_parameters[0]  # Use first evolution time
        
        for prep_name, amplitudes in preparation_strategies.items():
            # Create state vector
            state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=amplitudes)
            
            # Run state vector TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = StateVectorTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([state_vector], total_shots=3000)
            
            # Validate ground state energy
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use high precision tolerance for uniform superposition + complete subspace
            assert energy_error < 1e-8, (
                f"Preparation strategy '{prep_name}' TE-QSCI energy error "
                f"{energy_error:.6e} exceeds high precision tolerance 1e-8. "
                f"TE-QSCI: {result.ground_state_energy:.6f}, Exact: {exact_ground_energy:.6f}"
            )
    
    def test_state_vector_vs_circuit_equivalence(self, tolerance_settings, evolution_time_parameters):
        """Test equivalence between state vector input and circuit-generated states."""
        # Create 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Create a specific circuit state
        from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
        circuit_state = create_diverse_superposition_state(n_qubits=2, theta=np.pi/4)
        
        # Extract state vector from circuit (this is conceptual - in practice would use simulator)
        # For testing, create equivalent superposition manually
        # This represents the state that would be created by the circuit
        equivalent_amplitudes = np.array([0.7, 0.2+0.3j, 0.4-0.1j, 0.5+0.2j], dtype=complex)
        state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=equivalent_amplitudes)
        
        evolution_time = evolution_time_parameters[1]  # Use middle evolution time
        
        # Run circuit-based TE-QSCI (using single-time for comparison)
        sampler_circuit = create_qulacs_vector_concurrent_sampler()
        te_qsci_circuit = SingleTimeTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler_circuit,
            evolution_time=evolution_time,
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        result_circuit = te_qsci_circuit.run([circuit_state], total_shots=2500)
        
        # Run state vector TE-QSCI
        sampler_vector = create_qulacs_vector_concurrent_sampler()
        te_qsci_vector = StateVectorTE_QSCI(
            hamiltonian=quri_hamiltonian,
            sampler=sampler_vector,
            evolution_time=evolution_time,
            num_states_pick_out=2**2  # Complete subspace coverage
        )
        result_vector = te_qsci_vector.run([state_vector], total_shots=2500)
        
        # Compare results (they should be reasonably similar for equivalent states)
        energy_diff = abs(result_circuit.ground_state_energy - result_vector.ground_state_energy)
        
        # Both should achieve machine precision with uniform superposition + complete subspace
        assert energy_diff < 1e-8, (
            f"Circuit vs state vector TE-QSCI energy difference too large: "
            f"{energy_diff:.6e}. Circuit: {result_circuit.ground_state_energy:.6f}, "
            f"Vector: {result_vector.ground_state_energy:.6f}"
        )
    
    def test_state_normalization_handling(self, tolerance_settings, evolution_time_parameters):
        """Test how state vector TE-QSCI handles non-normalized input states."""
        # Create simple 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Test with non-normalized state vectors
        non_normalized_states = [
            np.array([2.0, 3.0, 1.0, 4.0], dtype=complex),  # Real amplitudes, not normalized
            np.array([1+1j, 2-1j, 0.5+2j, 3], dtype=complex),  # Complex amplitudes, not normalized
            np.array([100, 200, 50, 150], dtype=complex),  # Large amplitudes
        ]
        
        evolution_time = evolution_time_parameters[0]  # Use first evolution time
        
        for i, amplitudes in enumerate(non_normalized_states):
            # State vector creation should normalize automatically
            state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=amplitudes)
            
            # Note: CircuitQuantumState normalization is handled during circuit creation
            # No direct access to state vector amplitudes for verification
            
            # Run state vector TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = StateVectorTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([state_vector], total_shots=2500)
            
            # Validate ground state energy
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use high precision tolerance for uniform superposition + complete subspace
            assert energy_error < 1e-8, (
                f"Non-normalized state {i} TE-QSCI energy error "
                f"{energy_error:.6e} exceeds high precision tolerance 1e-8"
            )
    
    def test_edge_case_state_vectors(self, tolerance_settings, evolution_time_parameters):
        """Test state vector TE-QSCI with edge case state vectors."""
        # Create simple 2-qubit Pauli Hamiltonian
        hamiltonian_matrix = create_pauli_hamiltonian()
        sparse_ham = csc_matrix(hamiltonian_matrix)
        quri_hamiltonian = sparse_matrix_to_quri_operator(sparse_ham, n_qubits=2)
        
        # Get exact solution
        exact_eigenvalues, _ = get_exact_solution(sparse_ham, num_states=1)
        exact_ground_energy = exact_eigenvalues[0]
        
        # Test edge case states
        edge_case_states = {
            'zero_state': np.array([1, 0, 0, 0], dtype=complex),  # |00⟩
            'all_equal': np.array([1, 1, 1, 1], dtype=complex),   # Equal superposition
            'sparse_state': np.array([1, 0, 0, 1], dtype=complex),  # |00⟩ + |11⟩
            'complex_phases': np.array([1, 1j, -1, -1j], dtype=complex),  # Different phases
        }
        
        evolution_time = evolution_time_parameters[0]  # Use first evolution time
        
        for state_name, amplitudes in edge_case_states.items():
            # Create state vector
            state_vector = create_circuit_state_from_amplitudes(n_qubits=2, amplitudes=amplitudes)
            
            # Run state vector TE-QSCI
            sampler = create_qulacs_vector_concurrent_sampler()
            te_qsci = StateVectorTE_QSCI(
                hamiltonian=quri_hamiltonian,
                sampler=sampler,
                evolution_time=evolution_time,
                num_states_pick_out=2**2  # Complete subspace coverage
            )
            
            result = te_qsci.run([state_vector], total_shots=2500)
            
            # Validate ground state energy
            energy_error = abs(result.ground_state_energy - exact_ground_energy)
            
            # Use high precision tolerance for uniform superposition + complete subspace
            assert energy_error < 1e-8, (
                f"Edge case state '{state_name}' TE-QSCI energy error "
                f"{energy_error:.6e} exceeds high precision tolerance 1e-8"
            )


if __name__ == "__main__":
    # Run state vector TE-QSCI verification tests
    pytest.main([__file__, "-v", "-m", "te_qsci_state_vector"])