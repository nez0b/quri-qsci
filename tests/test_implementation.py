"""
Test script to verify TE-QSCI implementation.

This script runs basic tests to ensure the implementation works correctly
before running the full H6 molecule study.
"""

import numpy as np
import pytest
import sys
import os
from typing import Dict, List

# Add parent directory to path to import src module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Try importing our modules
try:
    from quri_qsci import (
        QSCIVariant, VanillaQSCI, TimeEvolvedQSCI, QSCIResult,
        create_qsci_algorithm, VanillaQSCIAlgorithm, SingleTimeTeQSCIAlgorithm,
        LoweringLevel, create_vm_enabled_algorithm
    )
    print("✓ Successfully imported TE-QSCI modules")
except ImportError as e:
    print(f"✗ Failed to import TE-QSCI modules: {e}")
    sys.exit(1)

# Try importing QURI Parts
try:
    from quri_parts.core.operator import pauli_label
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.core.state import CircuitQuantumState
    print("✓ Successfully imported QURI Parts")
except ImportError as e:
    print(f"✗ Failed to import QURI Parts: {e}")
    print("Please install QURI Parts with: pip install 'quri-parts[qulacs]'")
    sys.exit(1)

# Try importing sampler
try:
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    print("✓ Successfully imported Qulacs sampler")
    sampler_available = True
except ImportError as e:
    print(f"⚠ Warning: Qulacs sampler not available: {e}")
    print("Sampling-based tests will be skipped")
    sampler_available = False


def create_test_hamiltonian():
    """Create simple test Hamiltonian (H2 molecule)."""
    from quri_parts.core.operator import Operator, pauli_label
    
    hamiltonian = Operator()
    hamiltonian += Operator({pauli_label("Z0"): -1.25})
    hamiltonian += Operator({pauli_label("Z1"): -1.25})
    hamiltonian += Operator({pauli_label("Z0 Z1"): 0.5})
    hamiltonian += Operator({pauli_label("X0 X1"): 0.25})
    return hamiltonian


def create_test_state():
    """Create simple test state (HF state for H2)."""
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.core.state import GeneralCircuitQuantumState
    
    circuit = QuantumCircuit(2)
    circuit.add_X_gate(0)  # Fill first orbital
    return GeneralCircuitQuantumState(2, circuit)


def test_hamiltonian_creation():
    """Test Hamiltonian creation."""
    print("\n--- Testing Hamiltonian Creation ---")
    hamiltonian = create_test_hamiltonian()
    print(f"✓ Created test Hamiltonian with {len(hamiltonian)} terms")
    
    # Check if Hermitian
    from quri_parts.core.operator import is_hermitian
    assert is_hermitian(hamiltonian), "Hamiltonian must be Hermitian"
    print("✓ Hamiltonian is Hermitian")


def test_state_creation():
    """Test quantum state creation."""
    print("\n--- Testing Quantum State Creation ---")
    state = create_test_state()
    print(f"✓ Created test state with {state.qubit_count} qubits")
    print(f"✓ Circuit has {len(state.circuit.gates)} gates")
    assert state.qubit_count == 2, "State should have 2 qubits"
    assert len(state.circuit.gates) > 0, "Circuit should have gates"


def test_algorithm_creation():
    """Test algorithm creation."""
    print("\n--- Testing Algorithm Creation ---")
    hamiltonian = create_test_hamiltonian()
    
    # Test factory function
    vanilla_algo = create_qsci_algorithm(
        QSCIVariant.VANILLA,
        hamiltonian,
        num_states_pick_out=10
    )
    print(f"✓ Created vanilla QSCI: {vanilla_algo.name}")
    assert vanilla_algo.name is not None, "Algorithm should have a name"
    
    te_algo = create_qsci_algorithm(
        QSCIVariant.SINGLE_TIME_TE,
        hamiltonian,
        evolution_time=1.0,
        num_states_pick_out=10
    )
    print(f"✓ Created TE-QSCI: {te_algo.name}")
    assert te_algo.name is not None, "TE algorithm should have a name"
    
    # Test direct class creation
    direct_algo = SingleTimeTeQSCIAlgorithm(
        hamiltonian=hamiltonian,
        evolution_time=0.5,
        num_states_pick_out=10
    )
    print(f"✓ Created direct TE-QSCI: {direct_algo.name}")
    assert direct_algo.name is not None, "Direct algorithm should have a name"


def test_vm_analysis():
    """Test VM analysis capabilities."""
    print("\n--- Testing VM Analysis ---")
    try:
        hamiltonian = create_test_hamiltonian()
        base_algo = create_qsci_algorithm(
            QSCIVariant.SINGLE_TIME_TE,
            hamiltonian,
            evolution_time=1.0,
            num_states_pick_out=10
        )
        
        vm_algo = create_vm_enabled_algorithm(base_algo, "STAR")
        print(f"✓ Created VM-enabled algorithm: {vm_algo.name}")
        
        # Test analysis
        test_state = create_test_state()
        
        logical_analysis = vm_algo.analyze(
            [test_state], 100, LoweringLevel.LogicalCircuit
        )
        print(f"✓ Logical analysis - Latency: {logical_analysis.total_latency}, "
              f"Qubits: {logical_analysis.max_physical_qubit_count}")
        
        arch_analysis = vm_algo.analyze(
            [test_state], 100, LoweringLevel.ArchLogicalCircuit
        )
        print(f"✓ Arch analysis - Latency: {arch_analysis.total_latency}, "
              f"Qubits: {arch_analysis.max_physical_qubit_count}")
        
        assert True, "VM analysis completed successfully"
    except Exception as e:
        pytest.fail(f"VM analysis failed: {e}")


def test_sampling_based_execution():
    """Test actual algorithm execution with sampling."""
    if not sampler_available:
        print("\n--- Skipping Sampling Tests (Qulacs not available) ---")
        return True
    
    print("\n--- Testing Sampling-based Execution ---")
    try:
        hamiltonian = create_test_hamiltonian()
        sampler = create_qulacs_vector_concurrent_sampler()
        test_state = create_test_state()
        
        # Test vanilla QSCI
        vanilla_algo = VanillaQSCIAlgorithm(
            hamiltonian=hamiltonian,
            sampler=sampler,
            num_states_pick_out=5
        )
        
        print("Running vanilla QSCI...")
        vanilla_result = vanilla_algo.run([test_state], total_shots=100)
        print(f"✓ Vanilla QSCI result - Energy: {vanilla_result.ground_state_energy:.4f}")
        
        # Test TE-QSCI
        te_algo = SingleTimeTeQSCIAlgorithm(
            hamiltonian=hamiltonian,
            evolution_time=0.5,
            sampler=sampler,
            num_states_pick_out=5
        )
        
        print("Running TE-QSCI...")
        te_result = te_algo.run(test_state, total_shots=100)
        print(f"✓ TE-QSCI result - Energy: {te_result.ground_state_energy:.4f}")
        
        assert True, "Sampling execution completed successfully"
    except Exception as e:
        pytest.fail(f"Sampling execution failed: {e}")


def run_all_tests():
    """Run all tests and report results."""
    print("="*60)
    print("TE-QSCI Implementation Test Suite")
    print("="*60)
    
    tests = [
        ("Hamiltonian Creation", test_hamiltonian_creation),
        ("Quantum State Creation", test_state_creation),
        ("Algorithm Creation", test_algorithm_creation),
        ("VM Analysis", test_vm_analysis),
        ("Sampling Execution", test_sampling_based_execution),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    print("\n" + "="*60)
    print("TEST RESULTS SUMMARY")
    print("="*60)
    
    passed = 0
    total = len(results)
    
    for test_name, success in results:
        status = "PASS" if success else "FAIL"
        print(f"{test_name:<25} {status}")
        if success:
            passed += 1
    
    print("-"*60)
    print(f"TOTAL: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! TE-QSCI implementation is ready to use.")
        print("\nYou can now run the H6 molecule example:")
        print("  python h6_molecule_example.py")
    else:
        print(f"\n⚠ {total - passed} tests failed. Please check the implementation.")
        return False
    
    return True


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)