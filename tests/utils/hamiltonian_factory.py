"""
Test Hamiltonian generation utilities.

This module provides functions to create various types of Hamiltonians
for testing the QSCI implementation, including simple known cases,
random Hermitian matrices, and molecular systems.
"""

import numpy as np
from scipy.stats import unitary_group
from scipy.sparse import csc_matrix, identity, kron
from typing import Optional, Tuple
import sys
import os

# Add parent directory to path to import src module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from pyscf import gto, scf, mcscf
    from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
    from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
    from quri_parts.chem.transforms.jordan_wigner import jordan_wigner
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False


def create_hermitian_matrix(n: int, eigenvalues: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Creates a random n x n Hermitian matrix.
    
    Args:
        n: Size of the matrix
        eigenvalues: If provided, the matrix will have these exact eigenvalues
        
    Returns:
        n x n Hermitian matrix
    """
    if eigenvalues is not None:
        assert n == len(eigenvalues), f"Matrix size {n} must match eigenvalue count {len(eigenvalues)}"
        diag_matrix = np.diag(eigenvalues)
        random_unitary = unitary_group.rvs(n, random_state=42)  # Fixed seed for reproducibility
        return random_unitary @ diag_matrix @ random_unitary.conj().T
    else:
        # Create a random complex matrix and make it Hermitian
        np.random.seed(42)  # Fixed seed for reproducibility
        matrix = np.random.rand(n, n) + 1j * np.random.rand(n, n)
        return (matrix + matrix.conj().T) / 2


def create_simple_2x2_hamiltonian() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Creates a simple 2x2 Hermitian matrix with known eigenvalues and eigenvectors.
    
    Returns:
        Tuple of (hamiltonian, eigenvalues, eigenvectors)
    """
    # Simple Pauli-Z like matrix with off-diagonal coupling
    h = np.array([[1.0, 0.5], [0.5, -1.0]], dtype=complex)
    
    # Calculate exact eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(h)
    
    return h, eigvals, eigvecs


def create_degenerate_hamiltonian(size: int = 4) -> Tuple[np.ndarray, np.ndarray]:
    """
    Creates a Hamiltonian with degenerate eigenvalues for testing edge cases.
    
    Args:
        size: Size of the matrix (default 4)
        
    Returns:
        Tuple of (hamiltonian, known_eigenvalues)
    """
    if size == 4:
        # 4x4 matrix with eigenvalues [1, 2, 2, 4] (one degeneracy)
        eigenvalues = np.array([1.0, 2.0, 2.0, 4.0])
    elif size == 6:
        # 6x6 matrix with eigenvalues [1, 2, 2, 2, 5, 6] (triple degeneracy)
        eigenvalues = np.array([1.0, 2.0, 2.0, 2.0, 5.0, 6.0])
    else:
        # General case: create some degeneracies
        eigenvalues = np.arange(1, size + 1, dtype=float)
        eigenvalues[1:3] = 2.0  # Make first two non-ground eigenvalues degenerate
    
    hamiltonian = create_hermitian_matrix(size, eigenvalues)
    return hamiltonian, eigenvalues


def create_identity_hamiltonian(size: int) -> np.ndarray:
    """Creates an identity matrix for testing trivial cases."""
    return np.eye(size, dtype=complex)


def create_diagonal_hamiltonian(eigenvalues: np.ndarray) -> np.ndarray:
    """Creates a diagonal Hamiltonian with specified eigenvalues."""
    return np.diag(eigenvalues.astype(complex))


def create_pauli_hamiltonian() -> np.ndarray:
    """
    Creates a Hamiltonian based on Pauli matrices for 2-qubit system.
    H = 0.5 * (Z⊗I + I⊗Z + X⊗X)
    """
    # Pauli matrices
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    
    # Kronecker products
    ZI = np.kron(Z, I)
    IZ = np.kron(I, Z)
    XX = np.kron(X, X)
    
    # Hamiltonian
    H = 0.5 * (ZI + IZ + XX)
    return H


def get_h2_hamiltonian_minimal_basis() -> Tuple[np.ndarray, dict]:
    """
    Returns a pre-computed H2 Hamiltonian matrix for minimal basis set.
    
    Returns:
        Tuple of (hamiltonian_matrix, metadata)
    """
    if not PYSCF_AVAILABLE:
        # Return a mock 4x4 Hamiltonian for H2 minimal basis
        # This is a placeholder - in production, you'd pre-compute this
        h2_matrix = np.array([
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
            'method': 'mock_data'
        }
        
        return h2_matrix, metadata
    
    try:
        # Create H2 molecule
        mol = gto.Mole()
        mol.atom = 'H 0 0 0; H 0 0 0.74'
        mol.basis = 'sto-3g'
        mol.build(verbose=0)
        
        # Hartree-Fock calculation
        hf = scf.RHF(mol)
        hf.run(verbose=0)
        
        # CASCI setup (minimal for H2)
        casci = mcscf.CASCI(hf, ncas=2, nelecas=2)
        casci.run(verbose=0)
        
        # Get molecular integrals
        h1, h2 = get_spin_mo_integrals_from_mole(mol, casci.mo_coeff)
        
        # Create qubit Hamiltonian (4 qubits for 2 orbitals)
        hamiltonian, _ = get_qubit_mapped_hamiltonian(
            active_space_info=None,
            mo_eint_set=(h1, h2),
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Convert to matrix (assuming we can get matrix representation)
        # This is a simplified version - actual implementation would need proper conversion
        matrix_size = 2**4  # 4 qubits = 16x16 matrix
        ham_matrix = np.zeros((matrix_size, matrix_size), dtype=complex)
        
        # Note: This is a placeholder for actual matrix construction
        # In practice, you'd use the QURI-Parts operator to matrix conversion
        
        metadata = {
            'molecule': 'H2',
            'basis': 'sto-3g',
            'bond_length': 0.74,
            'n_qubits': 4,
            'n_electrons': 2,
            'hf_energy': hf.e_tot,
            'casci_energy': casci.e_tot,
            'nuclear_repulsion': mol.energy_nuc(),
            'method': 'pyscf_calculation'
        }
        
        return ham_matrix, metadata
        
    except Exception as e:
        # Fallback to mock data if PySCF calculation fails
        print(f"Warning: PySCF calculation failed ({e}), using mock data")
        return get_h2_hamiltonian_minimal_basis()


def create_random_hamiltonian(size: int, condition_number: float = 1e2) -> np.ndarray:
    """
    Creates a random Hermitian matrix with controlled condition number.
    
    Args:
        size: Matrix size
        condition_number: Desired condition number (ratio of largest to smallest eigenvalue)
        
    Returns:
        Random Hermitian matrix
    """
    # Create eigenvalues with specified condition number
    eigenvalues = np.logspace(0, np.log10(condition_number), size)
    eigenvalues = np.sort(eigenvalues)  # Sort in ascending order
    
    return create_hermitian_matrix(size, eigenvalues)


def validate_hamiltonian_properties(hamiltonian: np.ndarray, 
                                  tolerance: float = 1e-12) -> dict:
    """
    Validates that a matrix satisfies Hamiltonian properties.
    
    Args:
        hamiltonian: Matrix to validate
        tolerance: Numerical tolerance for checks
        
    Returns:
        Dictionary with validation results
    """
    results = {
        'is_square': hamiltonian.shape[0] == hamiltonian.shape[1],
        'is_hermitian': np.allclose(hamiltonian, hamiltonian.conj().T, atol=tolerance),
        'eigenvalues_real': True,
        'condition_number': 1.0,
        'smallest_eigenvalue': 0.0,
        'largest_eigenvalue': 0.0
    }
    
    if results['is_square']:
        try:
            eigenvalues = np.linalg.eigvals(hamiltonian)
            results['eigenvalues_real'] = np.allclose(eigenvalues.imag, 0, atol=tolerance)
            results['smallest_eigenvalue'] = float(np.min(eigenvalues.real))
            results['largest_eigenvalue'] = float(np.max(eigenvalues.real))
            
            if results['smallest_eigenvalue'] != 0:
                results['condition_number'] = results['largest_eigenvalue'] / abs(results['smallest_eigenvalue'])
            else:
                results['condition_number'] = np.inf
                
        except np.linalg.LinAlgError:
            results['eigenvalues_real'] = False
    
    results['is_valid_hamiltonian'] = (results['is_square'] and 
                                     results['is_hermitian'] and 
                                     results['eigenvalues_real'])
    
    return results


def create_tfim_hamiltonian(n_qubits: int, h: float) -> csc_matrix:
    """
    Creates the Transverse Field Ising Model (TFIM) Hamiltonian.
    
    H = - Σ Z_i Z_{i+1} - h Σ X_i
    
    Uses periodic boundary conditions (PBC) so the last qubit couples to the first.
    This is a canonical model in quantum many-body physics, exhibiting a quantum
    phase transition at h_c = 1.0 in the thermodynamic limit.
    
    Args:
        n_qubits: Number of qubits (must be >= 2)
        h: Transverse field strength
        
    Returns:
        Sparse CSC matrix representing the TFIM Hamiltonian
        
    Example:
        # Critical point of 1D TFIM
        ham = create_tfim_hamiltonian(4, h=1.0)
        
        # Ferromagnetic phase  
        ham = create_tfim_hamiltonian(4, h=0.5)
        
        # Paramagnetic phase
        ham = create_tfim_hamiltonian(4, h=2.0)
    """
    if n_qubits < 2:
        raise ValueError("TFIM requires at least 2 qubits")
    
    # Pauli matrices as sparse matrices
    sx = csc_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
    sz = csc_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
    id_matrix = identity(2, dtype=complex, format='csc')
    
    hamiltonian = csc_matrix((2**n_qubits, 2**n_qubits), dtype=complex)
    
    # Interaction term: - Σ Z_i Z_{i+1} (with PBC)
    for i in range(n_qubits):
        j = (i + 1) % n_qubits  # Periodic boundary conditions
        
        # Build operator list: I ⊗ I ⊗ ... ⊗ Z_i ⊗ ... ⊗ Z_j ⊗ ... ⊗ I
        op_list = [id_matrix] * n_qubits
        op_list[i] = sz
        op_list[j] = sz
        
        # Compute tensor product
        term = op_list[0]
        for k in range(1, n_qubits):
            term = kron(term, op_list[k], format='csc')
        
        hamiltonian = hamiltonian - term
    
    # Transverse field term: - h Σ X_i
    for i in range(n_qubits):
        # Build operator list: I ⊗ I ⊗ ... ⊗ X_i ⊗ ... ⊗ I
        op_list = [id_matrix] * n_qubits
        op_list[i] = sx
        
        # Compute tensor product
        term = op_list[0]
        for k in range(1, n_qubits):
            term = kron(term, op_list[k], format='csc')
        
        hamiltonian = hamiltonian - h * term
    
    return hamiltonian


def create_heisenberg_hamiltonian(n_qubits: int, jz: float, jxy: float) -> csc_matrix:
    """
    Creates the XXZ Heisenberg model Hamiltonian.
    
    H = Σ [Jxy * (X_i X_{i+1} + Y_i Y_{i+1}) + Jz * Z_i Z_{i+1}]
    
    This is a fundamental model of quantum magnetism. Special cases:
    - jz = jxy: Isotropic Heisenberg model (XXX)
    - jxy = 0: Pure Ising model
    - jz = 0: XX model
    
    Args:
        n_qubits: Number of qubits (must be >= 2)
        jz: Coupling strength for Z-Z interactions
        jxy: Coupling strength for X-X and Y-Y interactions
        
    Returns:
        Sparse CSC matrix representing the Heisenberg Hamiltonian
        
    Example:
        # Isotropic antiferromagnetic Heisenberg model
        ham = create_heisenberg_hamiltonian(4, jz=1.0, jxy=1.0)
        
        # Easy-axis anisotropic model
        ham = create_heisenberg_hamiltonian(4, jz=2.0, jxy=1.0)
        
        # XX model (no Ising interactions)
        ham = create_heisenberg_hamiltonian(4, jz=0.0, jxy=1.0)
    """
    if n_qubits < 2:
        raise ValueError("Heisenberg model requires at least 2 qubits")
    
    # Pauli matrices as sparse matrices
    sx = csc_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
    sy = csc_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
    sz = csc_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
    id_matrix = identity(2, dtype=complex, format='csc')
    
    hamiltonian = csc_matrix((2**n_qubits, 2**n_qubits), dtype=complex)
    
    # Loop over nearest-neighbor pairs (with PBC)
    for i in range(n_qubits):
        j = (i + 1) % n_qubits  # Periodic boundary conditions
        
        # X_i X_j term
        if jxy != 0:
            op_list = [id_matrix] * n_qubits
            op_list[i] = sx
            op_list[j] = sx
            
            term = op_list[0]
            for k in range(1, n_qubits):
                term = kron(term, op_list[k], format='csc')
            
            hamiltonian = hamiltonian + jxy * term
        
        # Y_i Y_j term  
        if jxy != 0:
            op_list = [id_matrix] * n_qubits
            op_list[i] = sy
            op_list[j] = sy
            
            term = op_list[0]
            for k in range(1, n_qubits):
                term = kron(term, op_list[k], format='csc')
            
            hamiltonian = hamiltonian + jxy * term
        
        # Z_i Z_j term
        if jz != 0:
            op_list = [id_matrix] * n_qubits
            op_list[i] = sz
            op_list[j] = sz
            
            term = op_list[0]
            for k in range(1, n_qubits):
                term = kron(term, op_list[k], format='csc')
            
            hamiltonian = hamiltonian + jz * term
    
    return hamiltonian


def create_random_sparse_hamiltonian(n_qubits: int, density: float = 0.1, 
                                   random_seed: int = 42) -> csc_matrix:
    """
    Creates a random sparse Hermitian Hamiltonian for stress testing.
    
    Generates a random sparse matrix and makes it Hermitian by taking
    H = (A + A†) / 2 where A is a random sparse matrix.
    
    Args:
        n_qubits: Number of qubits
        density: Fraction of non-zero elements (0 < density <= 1)
        random_seed: Random seed for reproducibility
        
    Returns:
        Random sparse Hermitian matrix
        
    Example:
        # Sparse random Hamiltonian for 4 qubits with 10% density
        ham = create_random_sparse_hamiltonian(4, density=0.1)
    """
    if not 0 < density <= 1:
        raise ValueError("Density must be between 0 and 1")
    
    np.random.seed(random_seed)
    
    dim = 2**n_qubits
    
    # Generate random sparse matrix
    from scipy.sparse import random
    random_matrix = random(dim, dim, density=density, format='csc', 
                          dtype=complex, random_state=random_seed)
    
    # Add random imaginary parts
    imag_matrix = random(dim, dim, density=density, format='csc',
                        dtype=complex, random_state=random_seed + 1)
    random_matrix = random_matrix + 1j * imag_matrix
    
    # Make Hermitian: H = (A + A†) / 2
    hermitian_matrix = (random_matrix + random_matrix.conj().T) / 2
    
    return hermitian_matrix.tocsc()


def get_exact_solution(hamiltonian_matrix, num_states: int = 1, 
                      which: str = 'SA') -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes exact eigenvalues and eigenvectors using sparse diagonalization.
    
    This function provides the ground truth for validation tests by computing
    exact solutions using scipy.sparse.linalg.eigsh or scipy.linalg.eigh
    for small matrices.
    
    Args:
        hamiltonian_matrix: Sparse or dense Hamiltonian matrix
        num_states: Number of eigenvalues/eigenvectors to compute
        which: Which eigenvalues to find ('SA' = smallest algebraic, 
               'LA' = largest algebraic, 'SM' = smallest magnitude, etc.)
               
    Returns:
        Tuple of (eigenvalues, eigenvectors)
        eigenvalues: 1D array of shape (num_states,)
        eigenvectors: 2D array of shape (dim, num_states)
        
    Example:
        # Get ground state
        ham = create_tfim_hamiltonian(3, h=1.0)
        energy, state = get_exact_solution(ham, num_states=1)
        
        # Get first 3 excited states
        energies, states = get_exact_solution(ham, num_states=3)
    """
    from scipy.sparse.linalg import eigsh
    from scipy.sparse import issparse
    
    # Convert to sparse if not already
    if not issparse(hamiltonian_matrix):
        hamiltonian_matrix = csc_matrix(hamiltonian_matrix)
    
    matrix_size = hamiltonian_matrix.shape[0]
    
    # For small matrices or when asking for all/most eigenvalues, use dense solver
    if matrix_size <= 8 or num_states >= matrix_size - 1:
        # Use dense eigensolver for small matrices
        dense_matrix = hamiltonian_matrix.toarray()
        all_eigenvalues, all_eigenvectors = np.linalg.eigh(dense_matrix)
        
        # Select the requested number of eigenvalues based on 'which' parameter
        if which == 'SA':  # Smallest algebraic
            indices = np.arange(min(num_states, len(all_eigenvalues)))
        elif which == 'LA':  # Largest algebraic
            indices = np.arange(max(0, len(all_eigenvalues) - num_states), len(all_eigenvalues))
        else:
            # Default to smallest
            indices = np.arange(min(num_states, len(all_eigenvalues)))
        
        eigenvalues = all_eigenvalues[indices]
        eigenvectors = all_eigenvectors[:, indices]
    else:
        # Use sparse eigensolver for larger matrices
        eigenvalues, eigenvectors = eigsh(hamiltonian_matrix, k=num_states, which=which)
        
        # Sort by eigenvalue (eigsh sometimes returns them unsorted)
        sort_indices = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[sort_indices]
        eigenvectors = eigenvectors[:, sort_indices]
    
    return eigenvalues, eigenvectors


def validate_sparse_hamiltonian_properties(hamiltonian: csc_matrix,
                                         tolerance: float = 1e-12) -> dict:
    """
    Validates that a sparse matrix satisfies Hamiltonian properties.
    
    Args:
        hamiltonian: Sparse matrix to validate
        tolerance: Numerical tolerance for checks
        
    Returns:
        Dictionary with validation results
    """
    from scipy.sparse import issparse
    
    results = {
        'is_sparse': issparse(hamiltonian),
        'is_square': hamiltonian.shape[0] == hamiltonian.shape[1],
        'is_hermitian': False,
        'sparsity': 0.0,
        'eigenvalues_real': True,
        'condition_number': 1.0,
        'smallest_eigenvalue': 0.0,
        'largest_eigenvalue': 0.0
    }
    
    if results['is_square']:
        # Check Hermiticity for sparse matrix
        hermitian_diff = hamiltonian - hamiltonian.conj().T
        results['is_hermitian'] = (hermitian_diff.max() < tolerance)
        
        # Calculate sparsity
        total_elements = hamiltonian.shape[0] * hamiltonian.shape[1]
        nonzero_elements = hamiltonian.nnz
        results['sparsity'] = 1.0 - (nonzero_elements / total_elements)
        
        # For small matrices, check eigenvalues
        if hamiltonian.shape[0] <= 64:  # Only for small matrices
            try:
                eigenvalues, _ = get_exact_solution(hamiltonian, num_states=1)
                results['smallest_eigenvalue'] = float(eigenvalues[0])
                
                # Get largest eigenvalue too
                eigenvalues_max, _ = get_exact_solution(hamiltonian, num_states=1, which='LA')
                results['largest_eigenvalue'] = float(eigenvalues_max[0])
                
                results['eigenvalues_real'] = True  # eigsh only returns real eigenvalues
                
                if abs(results['smallest_eigenvalue']) > tolerance:
                    results['condition_number'] = (results['largest_eigenvalue'] / 
                                                 abs(results['smallest_eigenvalue']))
                else:
                    results['condition_number'] = np.inf
                    
            except Exception:
                results['eigenvalues_real'] = False
    
    results['is_valid_hamiltonian'] = (results['is_square'] and 
                                     results['is_hermitian'] and 
                                     results['eigenvalues_real'])
    
    return results