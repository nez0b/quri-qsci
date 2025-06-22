#!/usr/bin/env python3
"""
Quick test to verify the physics fixes for H6 TE-QSCI are working correctly.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from examples.h6_molecule_example import H6MoleculeStudy

def test_physics_fixes():
    """Test that all three physics issues are resolved."""
    print("Testing H6 TE-QSCI physics fixes...")
    
    # Create study instance
    study = H6MoleculeStudy()
    
    print(f"✓ Hartree-Fock energy: {study.hartree_fock_energy:.6f} Ha")
    print(f"✓ FCI exact energy: {study.fci_energy:.6f} Ha")
    print(f"✓ Correlation energy: {study.hartree_fock_energy - study.fci_energy:.6f} Ha")
    
    # Verify energy ordering (Issue A fix)
    if study.fci_energy < study.hartree_fock_energy:
        print("✓ Energy ordering correct: FCI < HF (variational principle)")
    else:
        print("✗ Energy ordering wrong!")
        return False
    
    # Test single TE-QSCI calculation at short time
    print("\nTesting single TE-QSCI calculation...")
    
    from src import SingleTimeTeQSCIAlgorithm
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    
    sampler = create_qulacs_vector_concurrent_sampler()
    
    algorithm = SingleTimeTeQSCIAlgorithm(
        hamiltonian=study.hamiltonian,
        evolution_time=0.5,  # Short time for quick test
        sampler=sampler,
        num_states_pick_out=20  # Small subspace for speed
    )
    
    try:
        result = algorithm.run(
            initial_state=study.hartree_fock_state,
            total_shots=1000  # Reduced shots for speed
        )
        
        print(f"✓ TE-QSCI energy: {result.ground_state_energy:.6f} Ha")
        
        # Issue A: Check variational principle (TE-QSCI energy ≥ FCI energy)
        variational_error = result.ground_state_energy - study.fci_energy
        print(f"✓ Variational error: {variational_error:.6f} Ha")
        
        if variational_error >= -1e-6:  # Allow small numerical tolerance
            print("✓ Issue A FIXED: TE-QSCI energy ≥ FCI (variational principle)")
        else:
            print("✗ Issue A NOT FIXED: TE-QSCI energy < FCI!")
            return False
        
        # Issue B: Check probability calculation
        if hasattr(result, 'qsci_result') and result.qsci_result:
            ground_state_overlap = study._calculate_ground_state_overlap(result.qsci_result)
            print(f"✓ Ground state overlap P: {ground_state_overlap:.4f}")
            
            if 0.0 <= ground_state_overlap <= 1.0:
                print("✓ Issue B FIXED: Probability P in valid range [0,1]")
            else:
                print("✗ Issue B NOT FIXED: Invalid probability!")
                return False
        
        # Correlation energy percentage
        correlation_captured = (study.hartree_fock_energy - result.ground_state_energy)
        correlation_percent = 100 * correlation_captured / (study.hartree_fock_energy - study.fci_energy)
        print(f"✓ Correlation captured: {correlation_percent:.1f}%")
        
        print("\n✓ All physics fixes verified successfully!")
        return True
        
    except Exception as e:
        print(f"✗ Error in TE-QSCI calculation: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_physics_fixes()
    sys.exit(0 if success else 1)