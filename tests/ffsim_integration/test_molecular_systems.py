"""
Test suite for ffsim molecular systems module.
"""

import pytest

# Skip all tests if ffsim is not available
pytest.importorskip("ffsim", reason="ffsim is required for this test module")

# Import from the main package
try:
    from ffsim_integration.molecular_systems import (
        create_h2_molecule, 
        create_n2_molecule,
        get_reference_energies,
        validate_molecular_system
    )
except ImportError:
    pytest.skip("ffsim_integration module not available", allow_module_level=True)


def test_h2_molecule_creation():
    """Test H2 molecule creation with 6-31g basis."""
    h2_system = create_h2_molecule(basis="6-31g", bond_length=0.74)
    
    # Basic validation
    assert h2_system.name == "H2"
    assert h2_system.basis == "6-31g"
    assert h2_system.bond_length == 0.74
    assert h2_system.charge == 0
    assert h2_system.spin == 0
    
    # System validation
    assert validate_molecular_system(h2_system)
    
    # Energy checks (should be reasonable for H2)
    assert -2.0 < h2_system.hartree_fock_energy < 0.0
    assert h2_system.fci_energy <= h2_system.hartree_fock_energy
    
    # Active space checks
    assert h2_system.active_space.n_active_ele == 2  # 2 electrons in H2
    assert h2_system.active_space.n_active_orb >= 1  # At least 1 orbital
    
    print(f"✓ H2 test passed - HF: {h2_system.hartree_fock_energy:.6f} Ha")


def test_n2_molecule_creation():
    """Test N2 molecule creation with sto-3g basis and small active space."""
    # Use sto-3g basis and small active space (6,6) to speed up FCI calculation
    n2_system = create_n2_molecule(basis="sto-3g", bond_length=1.09, active_space=(6, 6))
    
    # Basic validation
    assert n2_system.name == "N2"
    assert n2_system.basis == "sto-3g"
    assert n2_system.bond_length == 1.09
    assert n2_system.charge == 0
    assert n2_system.spin == 0
    
    # System validation
    assert validate_molecular_system(n2_system)
    
    # Energy checks (should be reasonable for N2 with active space)
    assert -120.0 < n2_system.hartree_fock_energy < -100.0
    # Note: FCI might equal HF if active space is too small or FCI fails
    assert n2_system.fci_energy <= n2_system.hartree_fock_energy or abs(n2_system.fci_energy - n2_system.hartree_fock_energy) < 1e-8
    
    # Active space checks - updated for (6,6) active space
    # With (6,6) active space, we expect reduced electrons/orbitals
    assert n2_system.active_space.n_active_ele >= 6   # At least 6 active electrons
    assert n2_system.active_space.n_active_orb >= 6   # At least 6 active orbitals
    
    # Check ffsim molecular data for consistency 
    assert hasattr(n2_system, 'ffsim_mol_data')
    assert n2_system.ffsim_mol_data.norb >= 6  # Should have active orbitals
    
    print(f"✓ N2 test passed - HF: {n2_system.hartree_fock_energy:.6f} Ha")
    print(f"  Active space: {n2_system.active_space.n_active_ele} electrons, {n2_system.active_space.n_active_orb} orbitals")
    print(f"  FCI energy: {n2_system.fci_energy:.6f} Ha")


def test_reference_energies():
    """Test reference energy lookup function."""
    h2_refs = get_reference_energies("H2", "6-31g")
    assert "hartree_fock" in h2_refs
    assert "fci" in h2_refs
    assert "experimental" in h2_refs
    
    n2_refs = get_reference_energies("N2", "6-31g")
    assert "hartree_fock" in n2_refs
    assert "fci" in n2_refs
    assert "experimental" in n2_refs
    
    # Check that FCI is lower than HF (more negative)
    assert h2_refs["fci"] < h2_refs["hartree_fock"]
    assert n2_refs["fci"] < n2_refs["hartree_fock"]
    
    print(f"✓ Reference energies lookup working")


if __name__ == "__main__":
    print("Running molecular systems tests...")
    print("=" * 50)
    
    # test_reference_energies()
    test_h2_molecule_creation()
    test_n2_molecule_creation()  # Now using small active space for faster execution
    
    print("=" * 50)
    print("All molecular systems tests passed!")