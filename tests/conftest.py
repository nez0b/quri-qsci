"""
Shared pytest fixtures for QSCI testing framework.

This module contains fixtures that can be used across multiple test files
without explicit importing. These fixtures provide common test data and
setup for the QSCI tests.
"""

import pytest
import numpy as np
from .utils.hamiltonian_factory import (
    create_simple_2x2_hamiltonian,
    create_degenerate_hamiltonian,
    create_identity_hamiltonian,
    create_diagonal_hamiltonian,
    create_pauli_hamiltonian,
    create_random_hamiltonian
)


@pytest.fixture
def simple_2x2_hamiltonian():
    """A simple 2x2 Hermitian matrix with known eigenvalues and eigenvectors."""
    hamiltonian, eigenvalues, eigenvectors = create_simple_2x2_hamiltonian()
    return {
        'hamiltonian': hamiltonian,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'size': 2,
        'description': 'Simple 2x2 test case with known solution'
    }


@pytest.fixture
def degenerate_4x4_hamiltonian():
    """A 4x4 matrix with degenerate eigenvalues for testing edge cases."""
    hamiltonian, eigenvalues = create_degenerate_hamiltonian(size=4)
    return {
        'hamiltonian': hamiltonian,
        'known_eigenvalues': eigenvalues,
        'size': 4,
        'degeneracies': [(1, 1), (2, 2)],  # (eigenvalue, multiplicity)
        'description': '4x4 matrix with eigenvalues [1, 2, 2, 4]'
    }


@pytest.fixture
def degenerate_6x6_hamiltonian():
    """A 6x6 matrix with triple degenerate eigenvalues."""
    hamiltonian, eigenvalues = create_degenerate_hamiltonian(size=6)
    return {
        'hamiltonian': hamiltonian,
        'known_eigenvalues': eigenvalues,
        'size': 6,
        'degeneracies': [(1, 1), (2, 3), (5, 1), (6, 1)],  # (eigenvalue, multiplicity)
        'description': '6x6 matrix with eigenvalues [1, 2, 2, 2, 5, 6]'
    }


@pytest.fixture
def identity_hamiltonian():
    """Identity matrix for testing trivial cases."""
    size = 4
    hamiltonian = create_identity_hamiltonian(size)
    return {
        'hamiltonian': hamiltonian,
        'known_eigenvalues': np.ones(size),
        'size': size,
        'description': 'Identity matrix (all eigenvalues = 1)'
    }


@pytest.fixture
def diagonal_hamiltonian():
    """Diagonal matrix with known eigenvalues."""
    eigenvalues = np.array([0.5, 1.0, 2.0, 3.5])
    hamiltonian = create_diagonal_hamiltonian(eigenvalues)
    return {
        'hamiltonian': hamiltonian,
        'known_eigenvalues': eigenvalues,
        'size': len(eigenvalues),
        'description': 'Diagonal matrix with eigenvalues [0.5, 1.0, 2.0, 3.5]'
    }


@pytest.fixture
def pauli_hamiltonian():
    """2-qubit Hamiltonian based on Pauli matrices."""
    hamiltonian = create_pauli_hamiltonian()
    return {
        'hamiltonian': hamiltonian,
        'size': 4,
        'n_qubits': 2,
        'description': 'H = 0.5 * (Z⊗I + I⊗Z + X⊗X)'
    }


@pytest.fixture
def well_conditioned_hamiltonian():
    """Random Hamiltonian with good condition number."""
    size = 8
    hamiltonian = create_random_hamiltonian(size, condition_number=10.0)
    return {
        'hamiltonian': hamiltonian,
        'size': size,
        'condition_number': 10.0,
        'description': 'Random 8x8 Hermitian matrix with condition number ~10'
    }


@pytest.fixture
def ill_conditioned_hamiltonian():
    """Random Hamiltonian with poor condition number."""
    size = 6
    hamiltonian = create_random_hamiltonian(size, condition_number=1e6)
    return {
        'hamiltonian': hamiltonian,
        'size': size,
        'condition_number': 1e6,
        'description': 'Random 6x6 Hermitian matrix with condition number ~1e6'
    }


@pytest.fixture
def zero_eigenvalue_hamiltonian():
    """Hamiltonian with zero eigenvalues for testing numerical edge cases."""
    eigenvalues = np.array([0.0, 0.0, 1.0, 2.0])
    hamiltonian = create_diagonal_hamiltonian(eigenvalues)
    return {
        'hamiltonian': hamiltonian,
        'known_eigenvalues': eigenvalues,
        'size': len(eigenvalues),
        'description': 'Matrix with zero eigenvalues [0, 0, 1, 2]'
    }


@pytest.fixture(params=[1, 2, 3, 4, 5])
def eigenvalue_count_k(request):
    """Parametrized fixture for testing different values of k (number of eigenvalues)."""
    return request.param


@pytest.fixture(params=[True, False], ids=['sparse', 'dense'])
def diagonalization_method(request):
    """Parametrized fixture for testing both sparse and dense diagonalization."""
    return request.param


@pytest.fixture
def synthetic_measurement_counts():
    """Synthetic measurement counts for testing state selection."""
    return {
        'diverse_counts': {
            63: 500,   # Most frequent state
            31: 300,   # Second most frequent
            15: 150,   # Third most frequent
            7: 50,     # Fourth most frequent
            3: 30,     # Fifth most frequent
            1: 20,     # Sixth most frequent
            0: 10      # Least frequent
        },
        'single_state': {
            63: 1000   # Only one state (problematic case we discovered)
        },
        'uniform_counts': {
            i: 100 for i in range(10)  # Ten states with equal counts
        }
    }


@pytest.fixture
def tolerance_settings():
    """Standard numerical tolerances for testing."""
    return {
        'eigenvalue_atol': 1e-9,
        'eigenvector_atol': 1e-7,
        'hermitian_atol': 1e-12,
        'overlap_atol': 1e-7
    }