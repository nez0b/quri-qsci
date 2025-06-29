"""
Shared fixtures and configurations for TE-QSCI exact diagonalization tests.

This module provides common test configurations, fixtures, and utilities
for comprehensive TE-QSCI algorithm verification against exact solutions.
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from src.qsci_algorithms import VanillaQSCI, TimeEvolvedQSCI, SingleTimeTE_QSCI, TimeAverageTE_QSCI, StateVectorTE_QSCI
    from quri_parts.core.operator import Operator
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: QSCI implementation not available: {e}")
    QSCI_AVAILABLE = False


@pytest.fixture(scope="session")
def tolerance_settings():
    """Test tolerance settings for TE-QSCI verification."""
    return {
        'eigenvalue_atol': 1e-6,        # Absolute tolerance for exact methods
        'qsci_verification_atol': 1e-2,  # Realistic tolerance for QSCI approximation
        'qsci_loose_atol': 0.5,         # Loose tolerance for difficult cases
        'relative_tol': 1e-3            # Relative tolerance
    }


@pytest.fixture(scope="session") 
def evolution_time_parameters():
    """Standard evolution time parameters for TE-QSCI tests."""
    return [0.1, 0.3, 1.0, 3.0]  # Geometric progression


@pytest.fixture(scope="session")
def shot_count_parameters():
    """Shot count parameters for convergence analysis."""
    return [500, 1000, 2000]


@pytest.fixture(scope="session")
def test_hamiltonian_parameters():
    """Parameters for test Hamiltonians."""
    return {
        'tfim_fields': [0.5, 1.0, 2.0],  # Ferromagnetic, critical, paramagnetic
        'heisenberg_couplings': [
            (1.0, 1.0),  # Isotropic
            (2.0, 1.0),  # Easy-axis anisotropic
            (0.0, 1.0),  # XX model
        ],
        'system_sizes': [2, 3],  # Keep small for exact diagonalization
        'random_seed': 42
    }


@pytest.fixture
def qulacs_sampler():
    """Create QULACS sampler for QSCI algorithms."""
    if not QSCI_AVAILABLE:
        pytest.skip("QSCI not available")
    return create_qulacs_vector_concurrent_sampler()


def create_hartree_fock_state(n_qubits: int, n_electrons: int) -> GeneralCircuitQuantumState:
    """
    Create a Hartree-Fock initial state.
    
    For testing purposes, creates a simple HF-like state by placing
    electrons in the lowest energy orbitals (first n_electrons/2 spatial orbitals).
    
    Args:
        n_qubits: Number of qubits (should be even for spin orbitals)
        n_electrons: Number of electrons
        
    Returns:
        Quantum state representing HF configuration
    """
    circuit = QuantumCircuit(n_qubits)
    
    # Fill lowest energy orbitals first (spin-up then spin-down)
    for i in range(min(n_electrons, n_qubits)):
        circuit.add_X_gate(i)
    
    return GeneralCircuitQuantumState(n_qubits, circuit)


def create_diverse_superposition_state(n_qubits: int, theta: float = np.pi/6) -> GeneralCircuitQuantumState:
    """
    Create uniform superposition state for perfect QSCI sampling.
    
    This creates a uniform superposition state (H⊗H⊗...⊗H) that provides
    maximum measurement diversity and ensures complete subspace coverage
    for QSCI algorithms.
    
    Args:
        n_qubits: Number of qubits
        theta: Unused parameter (kept for compatibility)
        
    Returns:
        Uniform superposition quantum state
    """
    circuit = QuantumCircuit(n_qubits)
    
    # Apply Hadamard gates to all qubits for uniform superposition
    # This creates (|00...0⟩ + |00...1⟩ + ... + |11...1⟩)/√(2^n)
    for i in range(n_qubits):
        circuit.add_H_gate(i)
    
    return GeneralCircuitQuantumState(n_qubits, circuit)


def create_arbitrary_state(n_qubits: int, seed: int = 42) -> GeneralCircuitQuantumState:
    """
    Create an arbitrary quantum state for testing.
    
    Args:
        n_qubits: Number of qubits
        seed: Random seed for reproducibility
        
    Returns:
        Arbitrary quantum state
    """
    np.random.seed(seed)
    circuit = QuantumCircuit(n_qubits)
    
    # Random rotations
    for i in range(n_qubits):
        circuit.add_RY_gate(i, np.random.uniform(0, 2*np.pi))
        circuit.add_RZ_gate(i, np.random.uniform(0, 2*np.pi))
    
    # Random entanglement
    for i in range(n_qubits - 1):
        if np.random.random() > 0.5:
            circuit.add_CNOT_gate(i, (i + 1) % n_qubits)
    
    return GeneralCircuitQuantumState(n_qubits, circuit)


@pytest.fixture
def initial_state_factory():
    """Factory function for creating different types of initial states."""
    return {
        'hartree_fock': create_hartree_fock_state,
        'diverse_superposition': create_diverse_superposition_state,
        'arbitrary': create_arbitrary_state
    }


# Pytest markers for organizing tests
def pytest_configure(config):
    """Configure pytest markers for TE-QSCI tests."""
    config.addinivalue_line(
        "markers", "te_qsci_single_time: Single-time TE-QSCI algorithm tests"
    )
    config.addinivalue_line(
        "markers", "te_qsci_time_average: Time-average TE-QSCI algorithm tests"
    )
    config.addinivalue_line(
        "markers", "te_qsci_state_vector: State vector TE-QSCI algorithm tests"
    )
    config.addinivalue_line(
        "markers", "initial_state_comparison: Initial state comparison tests"
    )
    config.addinivalue_line(
        "markers", "exact_diagonalization: Tests using exact diagonalization verification"
    )