"""
Molecular system setup utilities for ffsim/UCJ ansatz integration.

This module provides utilities for creating molecular systems (H2, N2) using PySCF
and converting them to formats compatible with both ffsim and QURI Parts.

This module requires ffsim to be installed: pip install quri-qsci[ffsim]
"""

from typing import Dict, Optional, Tuple, Any
import numpy as np
from dataclasses import dataclass

# Quantum chemistry imports
from pyscf import gto, scf, fci
from quri_parts.core.operator import Operator
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian
from quri_parts.openfermion.transforms import jordan_wigner
from quri_parts.chem.mol import ActiveSpace, MolecularOrbitals

# ffsim imports - with proper error handling
try:
    import ffsim
    FFSIM_AVAILABLE = True
except ImportError:
    FFSIM_AVAILABLE = False
    
def _require_ffsim():
    """Check if ffsim is available and raise helpful error if not."""
    if not FFSIM_AVAILABLE:
        raise ImportError(
            "ffsim is required for this functionality. "
            "Install with: pip install quri-qsci[ffsim]"
        )


@dataclass
class MolecularSystem:
    """Container for molecular system data."""
    name: str
    geometry: str
    basis: str
    charge: int
    spin: int
    bond_length: float
    
    # PySCF objects
    mole: Any
    scf_result: Any
    
    # Quantum chemistry data
    active_space: Any
    mo_integrals: Any
    
    # Energy references
    hartree_fock_energy: float
    fci_energy: float
    
    # Hamiltonian representations
    quri_hamiltonian: Operator
    ffsim_mol_data: Any
    jw_mapping: Any


def create_h2_molecule(
    basis: str = "6-31g", 
    bond_length: float = 0.74,
    charge: int = 0,
    spin: int = 0
) -> MolecularSystem:
    """Create H2 molecule system for testing and validation.
    
    Args:
        basis: Basis set (default: 6-31g)
        bond_length: H-H bond length in Angstrom (default: 0.74)
        charge: Molecular charge (default: 0)
        spin: Spin multiplicity (default: 0 for singlet)
        
    Returns:
        MolecularSystem object with all necessary data
        
    Raises:
        ImportError: If ffsim is not installed
    """
    _require_ffsim()
    print(f"Creating H2 molecule with {basis} basis, bond length {bond_length} Å")
    
    # Create geometry string
    geometry = f"H 0 0 0; H 0 0 {bond_length}"
    
    try:
        # Create PySCF molecule
        geometry_list = [
            ["H", (0.0, 0.0, 0.0)],
            ["H", (0.0, 0.0, bond_length)]
        ]
        
        mole = gto.M(
            atom=geometry_list,
            basis=basis,
            charge=charge,
            spin=spin,
            verbose=0
        )
        
        # Perform Hartree-Fock calculation
        mf = scf.RHF(mole)
        mf.verbose = 0
        mf.run()
        
        if not mf.converged:
            raise RuntimeError("Hartree-Fock calculation did not converge")
        
        # Generate molecular orbital integrals for QURI Parts
        active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
        
        # Create QURI Parts Hamiltonian
        quri_hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
            active_space,
            mo_eint_set,
            sz=None,
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Calculate FCI reference energy
        fci_energy = _calculate_fci_energy(mole, mf)
        
        # Create ffsim molecular data
        ffsim_mol_data = ffsim.MolecularData.from_scf(mf)
        
        print(f"✓ H2 system created successfully")
        print(f"  - Active space: {active_space.n_active_ele} electrons, {active_space.n_active_orb} orbitals")
        print(f"  - Qubit count: {2 * active_space.n_active_orb}")
        print(f"  - HF energy: {mf.e_tot:.6f} Ha")
        print(f"  - FCI energy: {fci_energy:.6f} Ha")
        print(f"  - Correlation energy: {mf.e_tot - fci_energy:.6f} Ha")
        
        return MolecularSystem(
            name="H2",
            geometry=geometry,
            basis=basis,
            charge=charge,
            spin=spin,
            bond_length=bond_length,
            mole=mole,
            scf_result=mf,
            active_space=active_space,
            mo_integrals=mo_eint_set,
            hartree_fock_energy=mf.e_tot,
            fci_energy=fci_energy,
            quri_hamiltonian=quri_hamiltonian,
            ffsim_mol_data=ffsim_mol_data,
            jw_mapping=jw_mapping
        )
        
    except Exception as e:
        print(f"Error creating H2 molecule: {e}")
        raise


def create_n2_molecule(
    basis: str = "6-31g",
    bond_length: float = 1.0,
    charge: int = 0,
    spin: int = 0,
    active_space: Optional[Tuple[int, int]] = None
) -> MolecularSystem:
    """Create N2 molecule system for benchmark studies.
    
    Args:
        basis: Basis set (default: 6-31g) 
        bond_length: N-N bond length in Angstrom (default: 1.0)
        charge: Molecular charge (default: 0)
        spin: Spin multiplicity (default: 0 for singlet)
        active_space: (n_electrons, n_orbitals) for active space (None for full space)
        
    Returns:
        MolecularSystem object with all necessary data
        
    Raises:
        ImportError: If ffsim is not installed
    """
    _require_ffsim()
    print(f"Creating N2 molecule with {basis} basis, bond length {bond_length} Å")
    
    # Create geometry string
    geometry = f"N 0 0 0; N 0 0 {bond_length}"
    
    try:
        # Create PySCF molecule
        geometry_list = [
            ["N", (0.0, 0.0, 0.0)],
            ["N", (0.0, 0.0, bond_length)]
        ]
        
        mole = gto.M(
            atom=geometry_list,
            basis=basis,
            charge=charge,
            spin=spin,
            verbose=0
        )
        
        # Perform Hartree-Fock calculation
        mf = scf.RHF(mole)
        mf.verbose = 0
        mf.run()
        
        if not mf.converged:
            raise RuntimeError("Hartree-Fock calculation did not converge")
        
        # Generate molecular orbital integrals for QURI Parts
        active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mf.mo_coeff)
        
        # Create QURI Parts Hamiltonian
        quri_hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
            active_space,
            mo_eint_set,
            sz=None,
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Calculate FCI reference energy
        fci_energy = _calculate_fci_energy(mole, mf)
        
        # Create ffsim molecular data
        ffsim_mol_data = ffsim.MolecularData.from_scf(mf)
        
        print(f"✓ N2 system created successfully")
        print(f"  - Active space: {active_space.n_active_ele} electrons, {active_space.n_active_orb} orbitals")
        print(f"  - Qubit count: {2 * active_space.n_active_orb}")
        print(f"  - HF energy: {mf.e_tot:.6f} Ha")
        print(f"  - FCI energy: {fci_energy:.6f} Ha")
        print(f"  - Correlation energy: {mf.e_tot - fci_energy:.6f} Ha")
        print(f"  - ffsim molecule details:")
        print(f"    - norb: {ffsim_mol_data.norb}")
        print(f"    - nelec: {ffsim_mol_data.nelec}")
        print(f"    - Expected ffsim Fock space: C({ffsim_mol_data.norb},{ffsim_mol_data.nelec[0]}) * C({ffsim_mol_data.norb},{ffsim_mol_data.nelec[1]})")
        
        from scipy.special import comb
        expected_ffsim_dim = int(comb(ffsim_mol_data.norb, ffsim_mol_data.nelec[0]) * comb(ffsim_mol_data.norb, ffsim_mol_data.nelec[1]))
        print(f"    - Expected ffsim dimension: {expected_ffsim_dim}")
        
        return MolecularSystem(
            name="N2",
            geometry=geometry,
            basis=basis,
            charge=charge,
            spin=spin,
            bond_length=bond_length,
            mole=mole,
            scf_result=mf,
            active_space=active_space,
            mo_integrals=mo_eint_set,
            hartree_fock_energy=mf.e_tot,
            fci_energy=fci_energy,
            quri_hamiltonian=quri_hamiltonian,
            ffsim_mol_data=ffsim_mol_data,
            jw_mapping=jw_mapping
        )
        
    except Exception as e:
        print(f"Error creating N2 molecule: {e}")
        raise


def create_n2_homo_lumo_molecule(
    basis: str = "sto-3g",
    bond_length: float = 1.0,
    charge: int = 0,
    spin: int = 0,
    active_space_size: int = 4
) -> MolecularSystem:
    """Create N2 molecule with HOMO-LUMO focused active space for efficient optimization.
    
    Args:
        basis: Basis set (default: sto-3g)
        bond_length: N-N bond length in Angstrom (default: 1.0)
        charge: Molecular charge (default: 0)
        spin: Spin multiplicity (default: 0 for singlet)
        active_space_size: Number of orbitals around HOMO-LUMO (default: 4)
        
    Returns:
        MolecularSystem object with small active space for efficient optimization
    """
    print(f"Creating N2 molecule with HOMO-LUMO active space ({active_space_size} orbitals) in {basis} basis")
    
    try:
        # Build N2 molecule
        mol = gto.Mole()
        mol.build(
            atom=[["N", (0, 0, 0)], ["N", (bond_length, 0, 0)]],
            basis=basis,
            symmetry="Dooh",
            charge=charge,
            spin=spin,
            verbose=0
        )
        
        # Perform Hartree-Fock calculation
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.run()
        
        if not mf.converged:
            raise RuntimeError("Hartree-Fock calculation did not converge")
        
        # Calculate HOMO-LUMO indices for active space selection
        n_electrons = mol.nelectron
        homo_idx = n_electrons // 2 - 1  # Highest occupied orbital (0-indexed)
        
        # Define active space around HOMO-LUMO region
        half_active = active_space_size // 2
        start_orbital = max(0, homo_idx - half_active + 1)
        end_orbital = min(mol.nao_nr(), start_orbital + active_space_size)
        active_orbitals = list(range(start_orbital, end_orbital))
        
        # Calculate electrons in active space
        # Count occupied orbitals in active space
        n_occ_in_active = 0
        for orb_idx in active_orbitals:
            if orb_idx <= homo_idx:
                n_occ_in_active += 1
        
        n_active_electrons = n_occ_in_active * 2  # Each occupied orbital has 2 electrons
        n_active_orbitals = len(active_orbitals)
        
        print(f"HOMO-LUMO active space configuration:")
        print(f"  - Total orbitals: {mol.nao_nr()}")
        print(f"  - HOMO index: {homo_idx}")
        print(f"  - Active orbital indices: {active_orbitals}")
        print(f"  - Active orbitals: {n_active_orbitals}")
        print(f"  - Active electrons: {n_active_electrons}")
        
        # Create a modified molecule object for the active space
        from copy import deepcopy
        
        # Extract active space MO coefficients
        active_mo_coeff = mf.mo_coeff[:, active_orbitals]
        
        # Create a temporary SCF object with modified MO coefficients for active space
        temp_mf = deepcopy(mf)
        temp_mf.mo_coeff = active_mo_coeff
        
        # Create a temporary molecule with adjusted dimensions
        temp_mol = deepcopy(mol)
        temp_mol.nao_nr = lambda *args: n_active_orbitals
        temp_mol.nelectron = n_active_electrons
        
        # Generate molecular orbital integrals using QURI Parts for active space  
        _, mo_integrals = get_spin_mo_integrals_from_mole(temp_mol, active_mo_coeff)
        
        # Create active space object
        active_space_obj = ActiveSpace(
            n_active_ele=n_active_electrons,
            n_active_orb=n_active_orbitals
        )
        
        # Calculate frozen orbital energy (everything outside active space)
        frozen_energy = 0.0
        
        # Calculate energy of doubly occupied orbitals below active space
        if start_orbital > 0:
            frozen_occ_coeff = mf.mo_coeff[:, :start_orbital]
            h1_ao = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
            frozen_h1 = np.einsum('pi,pq,qj->ij', frozen_occ_coeff, h1_ao, frozen_occ_coeff)
            frozen_energy += 2 * np.trace(frozen_h1)
        
        # Add nuclear repulsion and any other constant terms
        total_constant = mo_integrals.const + frozen_energy
        mo_integrals_updated = mo_integrals._replace(const=total_constant)
        
        # Create QURI Parts Hamiltonian
        quri_hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
            active_space_obj,
            mo_integrals_updated,
            sz=None,
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Calculate FCI reference energy for full system
        fci_energy = _calculate_fci_energy(mol, mf)
        
        # Create ffsim molecular data for active space
        active_ffsim_data = ffsim.MolecularData(
            one_body_tensor=mo_integrals_updated.mo_1e_int,
            two_body_tensor=mo_integrals_updated.mo_2e_int,
            constant=total_constant,
            norb=n_active_orbitals,
            nelec=(n_active_electrons//2, n_active_electrons//2)  # Assume closed shell
        )
        
        print(f"✓ N2 HOMO-LUMO active space system created successfully")
        print(f"  - Active space: {n_active_electrons} electrons, {n_active_orbitals} orbitals")
        print(f"  - Qubit count: {2 * n_active_orbitals}")
        print(f"  - HF energy: {mf.e_tot:.6f} Ha")
        print(f"  - FCI energy: {fci_energy:.6f} Ha")
        print(f"  - Active space constant: {total_constant:.6f} Ha")
        
        # Estimate expected parameter count
        if n_active_orbitals == 4:
            print(f"  - Estimated LUCJ parameters: ~30-40 (much more feasible!)")
            print(f"  - Estimated UCJ parameters: ~50-70 (much more feasible!)")
        
        return MolecularSystem(
            name="N2_HOMO_LUMO",
            geometry=f"N 0 0 0; N 0 0 {bond_length}",
            basis=basis,
            charge=charge,
            spin=spin,
            bond_length=bond_length,
            mole=mol,
            scf_result=mf,
            active_space=active_space_obj,
            mo_integrals=mo_integrals_updated,
            hartree_fock_energy=mf.e_tot,
            fci_energy=fci_energy,
            quri_hamiltonian=quri_hamiltonian,
            ffsim_mol_data=active_ffsim_data,
            jw_mapping=jw_mapping
        )
        
    except Exception as e:
        print(f"Error creating N2 HOMO-LUMO molecule: {e}")
        raise


def create_n2_active_space_molecule(
    n_frozen: int = 2,
    basis: str = "6-31g", 
    bond_length: float = 1.0,
    charge: int = 0,
    spin: int = 0,
    use_homo_lumo: bool = False,
    homo_lumo_size: int = 4
) -> MolecularSystem:
    """Create N2 molecule with active space using frozen core orbitals or HOMO-LUMO.
    
    Args:
        n_frozen: Number of frozen core orbitals (default: 2)
        basis: Basis set (default: 6-31g)
        bond_length: N-N bond length in Angstrom (default: 1.0)
        charge: Molecular charge (default: 0)
        spin: Spin multiplicity (default: 0 for singlet)
        use_homo_lumo: If True, use HOMO-LUMO focused active space (default: False)
        homo_lumo_size: Number of orbitals around HOMO-LUMO (default: 4)
        
    Returns:
        MolecularSystem object with active space
    """
    print(f"Creating N2 molecule with frozen core ({n_frozen} frozen orbitals) in {basis} basis")
    
    try:
        # Specify molecule properties
        open_shell = False
        spin_sq = 0
        
        # Build N2 molecule using exact config specified
        mol = gto.Mole()
        mol.build(
            atom=[["N", (0, 0, 0)], ["N", (bond_length, 0, 0)]],
            basis=basis,
            symmetry="Dooh",
            charge=charge,
            spin=spin,
            verbose=0
        )
        
        # Perform Hartree-Fock calculation
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.run()
        
        if not mf.converged:
            raise RuntimeError("Hartree-Fock calculation did not converge")
        
        # Define active space
        if use_homo_lumo:
            # Use HOMO-LUMO focused active space
            n_electrons = mol.nelectron
            homo_idx = n_electrons // 2 - 1  # Highest occupied orbital (0-indexed)
            
            # Define active space around HOMO-LUMO region
            half_active = homo_lumo_size // 2
            start_orbital = max(0, homo_idx - half_active + 1)
            end_orbital = min(mol.nao_nr(), start_orbital + homo_lumo_size)
            active_space_orbitals = list(range(start_orbital, end_orbital))
            
            # Calculate electrons in active space
            n_occ_in_active = 0
            for orb_idx in active_space_orbitals:
                if orb_idx <= homo_idx:
                    n_occ_in_active += 1
            
            n_active_electrons = n_occ_in_active * 2  # Each occupied orbital has 2 electrons
            n_active_orbitals = len(active_space_orbitals)
            
            print(f"HOMO-LUMO active space configuration:")
            print(f"  - HOMO index: {homo_idx}")
            print(f"  - Active orbital indices: {active_space_orbitals}")
            
            # Extract HOMO-LUMO active space MO coefficients
            active_mo_coeff = mf.mo_coeff[:, active_space_orbitals]
        else:
            # Use frozen core approach
            active_space_orbitals = list(range(n_frozen, mol.nao_nr()))
            n_active_orbitals = len(active_space_orbitals)
            n_active_electrons = mol.nelectron - 2 * n_frozen  # Remove electrons in frozen orbitals
            
            # Extract active space MO coefficients
            active_mo_coeff = mf.mo_coeff[:, n_frozen:]
        
        print(f"Active space configuration:")
        print(f"  - Total orbitals: {mol.nao_nr()}")
        print(f"  - Active orbitals: {n_active_orbitals}")
        print(f"  - Active electrons: {n_active_electrons}")
        
        # Create a temporary SCF object with modified MO coefficients for active space
        from copy import deepcopy
        temp_mf = deepcopy(mf)
        temp_mf.mo_coeff = active_mo_coeff
        
        # Generate molecular orbital integrals using QURI Parts for active space  
        _, mo_integrals = get_spin_mo_integrals_from_mole(mol, active_mo_coeff)
        
        # Create a new active space object with correct dimensions
        active_space_obj = ActiveSpace(
            n_active_ele=n_active_electrons,
            n_active_orb=n_active_orbitals
        )
        
        # Calculate frozen core energy contribution
        frozen_core_energy = 0.0
        if use_homo_lumo:
            # For HOMO-LUMO, calculate energy of orbitals outside active space
            if start_orbital > 0:
                frozen_mo_coeff = mf.mo_coeff[:, :start_orbital]
                h1_ao = mol.intor('int1e_kin') + mol.intor('int1e_nuc') 
                frozen_h1 = np.einsum('pi,pq,qj->ij', frozen_mo_coeff, h1_ao, frozen_mo_coeff)
                frozen_core_energy = 2 * np.trace(frozen_h1)  # Factor of 2 for doubly occupied
        elif n_frozen > 0:
            # Calculate energy of frozen orbitals
            frozen_mo_coeff = mf.mo_coeff[:, :n_frozen]
            h1_ao = mol.intor('int1e_kin') + mol.intor('int1e_nuc') 
            frozen_h1 = np.einsum('pi,pq,qj->ij', frozen_mo_coeff, h1_ao, frozen_mo_coeff)
            frozen_core_energy = 2 * np.trace(frozen_h1)  # Factor of 2 for doubly occupied
        
        # Add frozen core energy to the constant term
        total_constant = mo_integrals.const + frozen_core_energy
        
        # Update the original mo_integrals with new constant
        mo_integrals_updated = mo_integrals._replace(const=total_constant)
        
        # Create QURI Parts Hamiltonian
        quri_hamiltonian, jw_mapping = get_qubit_mapped_hamiltonian(
            active_space_obj,
            mo_integrals_updated,
            sz=None,
            fermion_qubit_mapping=jordan_wigner
        )
        
        # Calculate FCI reference energy
        fci_energy = _calculate_fci_energy(mol, mf)
        
        # Create ffsim molecular data for active space
        active_ffsim_data = ffsim.MolecularData(
            one_body_tensor=mo_integrals_updated.mo_1e_int,
            two_body_tensor=mo_integrals_updated.mo_2e_int,
            constant=total_constant,
            norb=n_active_orbitals,
            nelec=(n_active_electrons//2, n_active_electrons//2)  # Assume closed shell
        )
        
        print(f"✓ N2 active space system created successfully")
        print(f"  - Active space: {n_active_electrons} electrons, {n_active_orbitals} orbitals")
        print(f"  - Qubit count: {2 * n_active_orbitals}")
        print(f"  - HF energy: {mf.e_tot:.6f} Ha")
        print(f"  - FCI energy: {fci_energy:.6f} Ha")
        print(f"  - Active space constant: {total_constant:.6f} Ha")
        
        return MolecularSystem(
            name="N2_ACTIVE_SPACE",
            geometry=f"N 0 0 0; N 0 0 {bond_length}",
            basis=basis,
            charge=charge,
            spin=spin,
            bond_length=bond_length,
            mole=mol,
            scf_result=mf,
            active_space=active_space_obj,
            mo_integrals=mo_integrals_updated,
            hartree_fock_energy=mf.e_tot,
            fci_energy=fci_energy,
            quri_hamiltonian=quri_hamiltonian,
            ffsim_mol_data=active_ffsim_data,
            jw_mapping=jw_mapping
        )
        
    except Exception as e:
        print(f"Error creating N2 active space molecule: {e}")
        raise


def _calculate_fci_energy(mole: Any, scf_result: Any) -> float:
    """Calculate FCI ground state energy for exact reference.
    
    Args:
        mole: PySCF molecule object
        scf_result: SCF calculation result
        
    Returns:
        FCI ground state energy in Hartree
    """
    try:
        # Create FCI solver
        cisolver = fci.FCI(scf_result)
        
        # Solve for ground state
        fci_energy, fci_civec = cisolver.run()
        
        return fci_energy
        
    except Exception as e:
        print(f"Warning: FCI calculation failed ({e}), using HF energy as fallback")
        return scf_result.e_tot


def get_reference_energies(molecule: str, basis: str = "6-31g") -> Dict[str, float]:
    """Get reference energies for common molecules and basis sets.
    
    Args:
        molecule: Molecule name ("H2" or "N2")
        basis: Basis set name
        
    Returns:
        Dictionary with reference energies in Hartree
    """
    # Reference energies from literature/high-level calculations
    references = {
        "H2": {
            "6-31g": {
                "hartree_fock": -1.123,  # Approximate HF energy
                "fci": -1.135,           # Approximate FCI energy
                "experimental": -1.139   # Approximate experimental
            },
            "sto-3g": {
                "hartree_fock": -1.117,
                "fci": -1.127,
                "experimental": -1.139
            }
        },
        "N2": {
            "6-31g": {
                "hartree_fock": -108.9,  # Approximate HF energy
                "fci": -109.1,           # Approximate FCI energy  
                "experimental": -109.3   # Target benchmark energy
            },
            "sto-3g": {
                "hartree_fock": -107.5,
                "fci": -107.7,
                "experimental": -109.3
            }
        }
    }
    
    molecule_upper = molecule.upper()
    basis_lower = basis.lower()
    
    if molecule_upper in references and basis_lower in references[molecule_upper]:
        return references[molecule_upper][basis_lower]
    else:
        print(f"Warning: No reference energies available for {molecule} with {basis} basis")
        return {"hartree_fock": 0.0, "fci": 0.0, "experimental": 0.0}


def validate_molecular_system(system: MolecularSystem) -> bool:
    """Validate that a molecular system is properly constructed.
    
    Args:
        system: MolecularSystem to validate
        
    Returns:
        True if system is valid, False otherwise
    """
    try:
        # Check basic attributes
        assert system.mole is not None
        assert system.scf_result is not None
        assert system.scf_result.converged
        
        # Check active space
        assert system.active_space is not None
        assert system.active_space.n_active_ele > 0
        assert system.active_space.n_active_orb > 0
        
        # Check Hamiltonians
        assert system.quri_hamiltonian is not None
        assert system.ffsim_mol_data is not None
        
        # Check energy ordering (FCI should be lower than HF)
        assert system.fci_energy <= system.hartree_fock_energy
        
        print(f"✓ {system.name} molecular system validation passed")
        return True
        
    except AssertionError as e:
        print(f"✗ {system.name} molecular system validation failed: {e}")
        return False
    except Exception as e:
        print(f"✗ {system.name} molecular system validation error: {e}")
        return False