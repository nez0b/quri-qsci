#!/usr/bin/env python3
"""
CASCI calculation for N2 with (10e, 18o) active space for stochastic QSCI testing.

This script generates the reference data needed for stochastic QSCI calculations,
using a smaller active space (10e, 18o) for initial testing before scaling to (10e, 22o).

System specifications:
- Molecule: N2 at 1.09 Å bond length
- Basis set: cc-pVDZ 
- Active space: 10 electrons in 18 orbitals
- Method: CASCI (Complete Active Space Configuration Interaction)

Output files:
- n2_10e18o_1.09A_ci_vector.npy: CI coefficients
- n2_10e18o_1.09A_fci_energy.npy: FCI reference energy
- n2_10e18o_1.09A_h1.npy: One-electron integrals
- n2_10e18o_1.09A_h2.npy: Two-electron integrals  
- n2_10e18o_1.09A_ecore.npy: Core energy
- n2_10e18o_1.09A_config.npy: System configuration
"""

import numpy as np
import time
from pathlib import Path

# PySCF imports
try:
    from pyscf import gto, scf, mcscf
    from pyscf.fci import cistring
except ImportError as e:
    raise ImportError(
        f"PySCF is required for CASCI calculations: {e}\n"
        f"Install with: pip install pyscf"
    )

def create_n2_molecule(bond_length=1.09, basis="cc-pVDZ", verbose=0):
    """Create N2 molecule with specified parameters.
    
    Args:
        bond_length: N-N bond length in Angstroms
        basis: Basis set name
        verbose: PySCF verbosity level
        
    Returns:
        PySCF Mole object
    """
    mol = gto.M(
        atom=f'N 0 0 0; N 0 0 {bond_length}',
        basis=basis,
        charge=0,
        spin=0,  # Singlet ground state
        verbose=verbose
    )
    
    print(f"✓ Created N2 molecule:")
    print(f"  Bond length: {bond_length} Å")
    print(f"  Basis: {basis}")
    print(f"  Total electrons: {mol.nelectron}")
    print(f"  Total orbitals: {mol.nao}")
    
    return mol

def run_casci_calculation(mol, n_electrons=10, n_orbitals=18):
    """Run CASCI calculation for given molecule.
    
    Args:
        mol: PySCF molecule object
        n_electrons: Number of active electrons
        n_orbitals: Number of active orbitals
        
    Returns:
        CASCI solver object with results
    """
    print(f"\nRunning CASCI calculation...")
    print(f"Active space: ({n_electrons}e, {n_orbitals}o)")
    start_time = time.time()
    
    # Run HF calculation first
    print("  Running HF calculation...")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.run()
    
    print(f"  HF energy: {mf.e_tot:.6f} Ha")
    
    # Run CASCI calculation
    print("  Running CASCI calculation...")
    cas = mcscf.CASCI(mf, n_orbitals, n_electrons)
    cas.verbose = 0
    cas.run()
    
    execution_time = time.time() - start_time
    print(f"✓ CASCI calculation completed in {execution_time:.2f} seconds")
    print(f"  CASCI energy: {cas.e_tot:.6f} Ha")
    print(f"  CI vector shape: {cas.ci.shape}")
    print(f"  Total CI coefficients: {cas.ci.size}")
    
    # Calculate expected number of determinants
    try:
        na, nb = cas.nelecas if isinstance(cas.nelecas, tuple) else (cas.nelecas[0], cas.nelecas[1])
        n_alpha_strings = cistring.num_strings(n_orbitals, na)
        n_beta_strings = cistring.num_strings(n_orbitals, nb)
        expected_determinants = n_alpha_strings * n_beta_strings
        print(f"  Expected determinants: {expected_determinants}")
        print(f"  Alpha strings: {n_alpha_strings}")
        print(f"  Beta strings: {n_beta_strings}")
    except Exception as e:
        print(f"  Could not calculate expected determinants: {e}")
    
    return cas

def save_casci_data(cas, bond_length=1.09, output_dir="stochastic_data"):
    """Save CASCI results for stochastic QSCI calculations.
    
    Args:
        cas: CASCI solver object
        bond_length: Bond length for file naming
        output_dir: Output directory for data files
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    print(f"\nSaving CASCI data to {output_dir}/...")
    
    # Generate base filename
    bond_str = f"{bond_length:.2f}A".replace(".", "_")
    base_name = f"n2_10e18o_{bond_str}"
    
    # Get Hamiltonian components (memory-efficient approach)
    print("  Extracting Hamiltonian components...")
    h1, ecore = cas.get_h1eff()
    h2 = cas.get_h2eff()
    
    # Save all required data
    data_files = [
        (f"{base_name}_ci_vector.npy", cas.ci, "CI coefficients"),
        (f"{base_name}_fci_energy.npy", cas.e_tot, "FCI energy"),
        (f"{base_name}_h1.npy", h1, "One-electron integrals"),
        (f"{base_name}_h2.npy", h2, "Two-electron integrals"),
        (f"{base_name}_ecore.npy", ecore, "Core energy"),
        (f"{base_name}_config.npy", {
            'norb': cas.ncas, 
            'nelec': cas.nelecas,
            'bond_length': bond_length,
            'basis': 'cc-pVDZ',
            'n_electrons': 10,
            'n_orbitals': 18
        }, "System configuration")
    ]
    
    for filename, data, description in data_files:
        filepath = output_path / filename
        np.save(filepath, data)
        print(f"  ✓ Saved {description}: {filepath}")
    
    print(f"✓ All data saved with base name: {base_name}")
    
    # Print CI coefficient analysis
    ci_flat = cas.ci.flatten() if hasattr(cas.ci, 'flatten') else cas.ci
    ci_probs = np.abs(ci_flat)**2
    ci_probs_sorted = np.sort(ci_probs)[::-1]
    
    print(f"\nCI coefficient analysis:")
    print(f"  Largest probability: {ci_probs_sorted[0]:.6f} ({ci_probs_sorted[0]*100:.2f}%)")
    print(f"  Second largest: {ci_probs_sorted[1]:.6f} ({ci_probs_sorted[1]*100:.2f}%)")
    print(f"  Third largest: {ci_probs_sorted[2]:.6f} ({ci_probs_sorted[2]*100:.2f}%)")
    print(f"  Sum of top 10: {np.sum(ci_probs_sorted[:10]):.6f} ({np.sum(ci_probs_sorted[:10])*100:.2f}%)")
    print(f"  Sum of top 100: {np.sum(ci_probs_sorted[:100]):.6f} ({np.sum(ci_probs_sorted[:100])*100:.2f}%)")
    
    return base_name

def main():
    """Main function to run CASCI calculation and save data."""
    print("="*80)
    print("N2 CASCI CALCULATION FOR STOCHASTIC QSCI TESTING")
    print("="*80)
    print("System: N2 (10e, 18o) active space at 1.09 Å")
    print("Purpose: Generate reference data for stochastic QSCI convergence analysis")
    
    try:
        # Create molecule
        mol = create_n2_molecule(bond_length=1.09, basis="cc-pVDZ")
        
        # Run CASCI calculation
        cas = run_casci_calculation(mol, n_electrons=10, n_orbitals=18)
        
        # Save data for stochastic QSCI
        base_name = save_casci_data(cas, bond_length=1.09)
        
        print(f"\n✓ CASCI calculation completed successfully!")
        print(f"  Base name: {base_name}")
        print(f"  Files saved in: stochastic_data/")
        print(f"  Ready for stochastic QSCI testing")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in CASCI calculation: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)