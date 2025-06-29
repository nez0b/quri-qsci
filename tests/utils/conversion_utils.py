"""
Conversion utilities for QSCI testing framework.

This module provides functions to convert between different quantum chemistry
frameworks (QURI-Parts, Qiskit Nature, OpenFermion) for integration testing.
"""

import numpy as np
from typing import Dict, Tuple, Optional, Union
import warnings


def extract_operator_terms(operator) -> Dict[str, complex]:
    """
    Extract terms and coefficients from a quantum operator.
    
    This is a generic function that attempts to extract terms from
    various operator formats (QURI-Parts, OpenFermion, etc.)
    
    Args:
        operator: Quantum operator in various formats
        
    Returns:
        Dictionary mapping term strings to coefficients
    """
    terms = {}
    
    # Handle different operator types
    if hasattr(operator, 'terms'):
        # OpenFermion-style operator
        for term, coeff in operator.terms.items():
            terms[str(term)] = complex(coeff)
    elif hasattr(operator, '_terms'):
        # Alternative terms attribute
        for term, coeff in operator._terms.items():
            terms[str(term)] = complex(coeff)
    elif hasattr(operator, 'to_dict'):
        # Dictionary conversion method
        operator_dict = operator.to_dict()
        for term, coeff in operator_dict.items():
            terms[str(term)] = complex(coeff)
    else:
        # If we can't extract terms, return empty dict
        warnings.warn(f"Unknown operator type: {type(operator)}. Cannot extract terms.")
    
    return terms


def create_mock_qiskit_fermionic_op(terms: Dict[str, complex], num_spin_orbitals: int):
    """
    Create a mock Qiskit Nature FermionicOp for testing purposes.
    
    Since we may not have qiskit-nature available, this creates a mock object
    that has the essential properties for testing.
    
    Args:
        terms: Dictionary of fermionic terms and coefficients
        num_spin_orbitals: Number of spin orbitals
        
    Returns:
        Mock FermionicOp object
    """
    
    class MockFermionicOp:
        def __init__(self, terms_dict, num_spin_orbitals):
            self.terms = terms_dict
            self.num_spin_orbitals = num_spin_orbitals
        
        def eigenvalues(self):
            """Mock eigenvalue computation."""
            # This would normally interface with qiskit-addon-sqd
            # For testing, return mock eigenvalues
            return np.array([-1.5, -0.5, 0.5, 1.5])
        
        def to_matrix(self):
            """Convert to matrix representation (simplified)."""
            size = 2**self.num_spin_orbitals
            matrix = np.zeros((size, size), dtype=complex)
            
            # This is a simplified conversion - real implementation would be more complex
            for term_str, coeff in self.terms.items():
                # Add coefficient to diagonal (simplified)
                if term_str == "()":  # Identity term
                    matrix += coeff * np.eye(size)
            
            return matrix
    
    return MockFermionicOp(terms, num_spin_orbitals)


def convert_quri_to_qiskit_format(quri_operator, num_spin_orbitals: int):
    """
    Convert QURI-Parts operator to Qiskit Nature FermionicOp format.
    
    Args:
        quri_operator: QURI-Parts fermionic operator
        num_spin_orbitals: Number of spin orbitals
        
    Returns:
        Mock Qiskit Nature FermionicOp
    """
    # Extract terms from QURI-Parts operator
    terms = extract_operator_terms(quri_operator)
    
    # Convert term format if needed
    qiskit_terms = {}
    for term_str, coeff in terms.items():
        # Convert QURI-Parts term format to Qiskit format
        # This is a simplified conversion - real implementation would handle
        # the full fermionic operator term syntax
        qiskit_terms[term_str] = coeff
    
    return create_mock_qiskit_fermionic_op(qiskit_terms, num_spin_orbitals)


def create_identical_measurement_counts(probabilities: Dict[int, float], 
                                      total_shots: int) -> Dict[int, int]:
    """
    Create deterministic measurement counts from exact probabilities.
    
    This ensures both QSCI implementations start with identical data.
    
    Args:
        probabilities: Dictionary mapping bit strings to probabilities
        total_shots: Total number of measurement shots
        
    Returns:
        Dictionary mapping bit strings to counts
    """
    counts = {}
    remaining_shots = total_shots
    
    # Sort by probability (highest first) for deterministic allocation
    sorted_states = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    
    for i, (state, prob) in enumerate(sorted_states[:-1]):
        # Allocate shots proportionally, ensuring integer counts
        shot_count = int(prob * total_shots)
        counts[state] = shot_count
        remaining_shots -= shot_count
    
    # Allocate remaining shots to the last state
    if sorted_states:
        last_state = sorted_states[-1][0]
        counts[last_state] = remaining_shots
    
    return counts


def create_synthetic_diverse_counts(num_states: int = 10, 
                                  total_shots: int = 1000,
                                  decay_factor: float = 0.7) -> Dict[int, int]:
    """
    Create synthetic measurement counts with decreasing probabilities.
    
    This creates diverse measurement outcomes for testing QSCI state selection.
    
    Args:
        num_states: Number of different measured states
        total_shots: Total number of measurement shots
        decay_factor: How quickly probabilities decay (0 < decay_factor < 1)
        
    Returns:
        Dictionary mapping bit strings to counts
    """
    # Create exponentially decaying probabilities
    probabilities = {}
    for i in range(num_states):
        state = i  # Use simple integer bit patterns
        prob = decay_factor ** i
        probabilities[state] = prob
    
    # Normalize probabilities
    total_prob = sum(probabilities.values())
    for state in probabilities:
        probabilities[state] /= total_prob
    
    return create_identical_measurement_counts(probabilities, total_shots)


def validate_hamiltonian_equivalence(ham1: np.ndarray, 
                                   ham2: np.ndarray, 
                                   tolerance: float = 1e-10) -> bool:
    """
    Validate that two Hamiltonian matrices are equivalent.
    
    Args:
        ham1: First Hamiltonian matrix
        ham2: Second Hamiltonian matrix
        tolerance: Numerical tolerance for comparison
        
    Returns:
        True if matrices are equivalent within tolerance
    """
    if ham1.shape != ham2.shape:
        return False
    
    return np.allclose(ham1, ham2, atol=tolerance)


def compare_eigenvalue_lists(eigvals1: np.ndarray, 
                           eigvals2: np.ndarray,
                           tolerance: float = 1e-9) -> Tuple[bool, float]:
    """
    Compare two lists of eigenvalues, handling potential ordering differences.
    
    Args:
        eigvals1: First set of eigenvalues
        eigvals2: Second set of eigenvalues
        tolerance: Numerical tolerance for comparison
        
    Returns:
        Tuple of (are_equivalent, max_difference)
    """
    if len(eigvals1) != len(eigvals2):
        return False, np.inf
    
    # Sort both sets for comparison
    sorted1 = np.sort(eigvals1)
    sorted2 = np.sort(eigvals2)
    
    differences = np.abs(sorted1 - sorted2)
    max_diff = np.max(differences)
    
    return max_diff < tolerance, max_diff


def create_test_hamiltonian_matrix(system_type: str = "H2") -> Tuple[np.ndarray, dict]:
    """
    Create test Hamiltonian matrices for various molecular systems.
    
    Args:
        system_type: Type of molecular system ("H2", "H4", etc.)
        
    Returns:
        Tuple of (hamiltonian_matrix, metadata)
    """
    if system_type == "H2":
        # Simplified H2 Hamiltonian matrix (4x4 for minimal basis)
        h_matrix = np.array([
            [-1.0966,  0.0000,  0.0000,  0.1809],
            [ 0.0000, -0.4544,  0.1809,  0.0000],
            [ 0.0000,  0.1809, -0.4544,  0.0000],
            [ 0.1809,  0.0000,  0.0000, -1.0966]
        ], dtype=complex)
        
        metadata = {
            'molecule': 'H2',
            'basis': 'sto-3g',
            'bond_length': 0.74,
            'n_qubits': 4,
            'n_electrons': 2,
            'hf_energy': -1.83,
            'exact_ground_state': -1.86
        }
        
    elif system_type == "H4":
        # Simplified H4 Hamiltonian (8x8)
        size = 8
        h_matrix = np.random.rand(size, size) + 1j * np.random.rand(size, size)
        h_matrix = (h_matrix + h_matrix.conj().T) / 2  # Make Hermitian
        
        metadata = {
            'molecule': 'H4',
            'basis': 'sto-3g',
            'n_qubits': size,
            'n_electrons': 4,
            'method': 'mock_data'
        }
        
    else:
        raise ValueError(f"Unknown system type: {system_type}")
    
    return h_matrix, metadata


def mock_qiskit_sqd_diagonalization(hamiltonian: np.ndarray, 
                                   measurement_counts: Dict[int, int],
                                   num_states_select: int) -> Dict[str, Union[np.ndarray, list]]:
    """
    Mock implementation of qiskit-addon-sqd diagonalization for testing.
    
    Args:
        hamiltonian: Hamiltonian matrix
        measurement_counts: Measurement counts dictionary
        num_states_select: Number of states to select for subspace
        
    Returns:
        Dictionary with 'eigenvalues', 'selected_states', 'subspace_dim'
    """
    # Select most frequent states (mock state selection)
    sorted_states = sorted(measurement_counts.items(), key=lambda x: x[1], reverse=True)
    selected_states = [state for state, count in sorted_states[:num_states_select]]
    
    # Create truncated Hamiltonian for selected subspace
    # This is a simplified version - real implementation would be more complex
    subspace_size = len(selected_states)
    
    if subspace_size == 1:
        # Single state - return diagonal element
        state = selected_states[0]
        if state < hamiltonian.shape[0]:
            eigenvalue = hamiltonian[state, state].real
        else:
            eigenvalue = 0.0
        eigenvalues = np.array([eigenvalue])
    else:
        # Multiple states - create and diagonalize subspace Hamiltonian
        subspace_ham = np.zeros((subspace_size, subspace_size), dtype=complex)
        for i, state_i in enumerate(selected_states):
            for j, state_j in enumerate(selected_states):
                if state_i < hamiltonian.shape[0] and state_j < hamiltonian.shape[1]:
                    subspace_ham[i, j] = hamiltonian[state_i, state_j]
        
        eigenvalues, _ = np.linalg.eigh(subspace_ham)
    
    return {
        'eigenvalues': eigenvalues,
        'selected_states': selected_states,
        'subspace_dimension': subspace_size,
        'ground_state_energy': np.min(eigenvalues)
    }