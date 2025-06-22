#!/usr/bin/env python3
"""
Test script to verify exact vs Trotter evolution implementation.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src import SingleTimeTeQSCIAlgorithm
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
from quri_parts.core.operator import Operator, pauli_label
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.circuit import QuantumCircuit

def create_simple_hamiltonian():
    """Create a simple 4-qubit Hamiltonian for testing."""
    hamiltonian = Operator()
    
    # Add some Z terms
    for i in range(4):
        hamiltonian += Operator({pauli_label(f"Z{i}"): -0.5})
    
    # Add some XX interactions
    for i in range(3):
        hamiltonian += Operator({pauli_label(f"X{i} X{i+1}"): 0.25})
    
    hamiltonian.constant = -1.0
    return hamiltonian

def create_simple_initial_state():
    """Create a simple initial state."""
    circuit = QuantumCircuit(4)
    # Prepare |0101⟩ state
    circuit.add_X_gate(1)
    circuit.add_X_gate(3)
    return GeneralCircuitQuantumState(4, circuit)

def test_exact_vs_trotter():
    """Test exact vs Trotter evolution methods."""
    print("Testing exact vs Trotter evolution...")
    
    # Create test system
    hamiltonian = create_simple_hamiltonian()
    initial_state = create_simple_initial_state()
    sampler = create_qulacs_vector_concurrent_sampler()
    
    # Test parameters
    evolution_time = 0.5
    total_shots = 1000
    num_states = 20
    
    results = {}
    
    for method in ["exact", "trotter"]:
        print(f"\nTesting {method} evolution...")
        
        try:
            # Create algorithm
            algorithm = SingleTimeTeQSCIAlgorithm(
                hamiltonian=hamiltonian,
                evolution_time=evolution_time,
                sampler=sampler,
                num_states_pick_out=num_states,
                time_evolution_method=method
            )
            
            print(f"Algorithm: {algorithm.name}")
            
            # Run algorithm
            result = algorithm.run(
                initial_state=initial_state,
                total_shots=total_shots
            )
            
            energy = result.ground_state_energy
            subspace_dim = result.qsci_result.subspace_dimension
            
            results[method] = {
                "energy": energy,
                "subspace_dimension": subspace_dim,
                "eigenvalues": result.qsci_result.eigenvalues[:3]  # First 3 eigenvalues
            }
            
            print(f"✓ {method.capitalize()} energy: {energy:.6f}")
            print(f"✓ Subspace dimension: {subspace_dim}")
            print(f"✓ First 3 eigenvalues: {result.qsci_result.eigenvalues[:3]}")
            
        except Exception as e:
            print(f"✗ Error with {method} method: {e}")
            import traceback
            traceback.print_exc()
            results[method] = {"error": str(e)}
    
    # Compare results
    print("\n" + "="*50)
    print("COMPARISON RESULTS")
    print("="*50)
    
    if "exact" in results and "trotter" in results:
        exact_energy = results["exact"].get("energy")
        trotter_energy = results["trotter"].get("energy")
        
        if exact_energy is not None and trotter_energy is not None:
            energy_diff = trotter_energy - exact_energy
            print(f"Exact energy:    {exact_energy:.6f}")
            print(f"Trotter energy:  {trotter_energy:.6f}")
            print(f"Difference:      {energy_diff:.6f}")
            print(f"Relative error:  {abs(energy_diff/exact_energy)*100:.4f}%")
            
            if abs(energy_diff) < 0.01:  # Reasonable tolerance for Trotter approximation
                print("✓ Results are consistent!")
            else:
                print("⚠ Large difference detected - this may be expected for Trotter approximation")
        else:
            print("✗ Could not compare energies due to errors")
    
    return results

if __name__ == "__main__":
    test_exact_vs_trotter()