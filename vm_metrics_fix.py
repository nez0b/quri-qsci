#!/usr/bin/env python3
"""
VM Metrics Fix Analysis
=======================

Analyze how to properly fix the overhead and fidelity calculations.
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))

from quri_parts.circuit import QuantumCircuit
from quri_parts.backend.units import TimeValue, TimeUnit
from src.qsci_vm_analysis import VMFactory, QURIVMInterface
from src.qsci_algo_interface import LoweringLevel

def analyze_fix_options():
    """Analyze how to fix the metrics calculation."""
    print("VM Metrics Fix Analysis")
    print("=" * 40)
    
    # Create test circuit
    circuit = QuantumCircuit(4)
    circuit.add_H_gate(0)
    circuit.add_CNOT_gate(0, 1)
    
    # Test what happens with consistent VM types
    print("\n=== OPTION 1: Compare same error rates ===")
    
    error_rates = [0.0, 1e-4, 1e-3, 1e-2]
    
    for i, error_rate in enumerate(error_rates):
        print(f"\nError rate: {error_rate}")
        
        # Create two STAR VMs with same error rate for logical/arch comparison
        logical_vm = VMFactory.create_star_vm(0.0)  # Baseline STAR (no error)
        arch_vm = VMFactory.create_star_vm(error_rate)  # Test STAR with error
        
        logical_result = logical_vm.analyze(circuit)
        arch_result = arch_vm.analyze(circuit)
        
        print(f"  Logical STAR (ε=0.0): latency={logical_result.latency}, fidelity={logical_result.fidelity:.4f}")
        print(f"  Arch STAR (ε={error_rate}): latency={arch_result.latency}, fidelity={arch_result.fidelity:.4f}")
        
        # Calculate proper metrics
        if logical_result.latency and arch_result.latency:
            overhead = arch_result.latency.in_ns() / logical_result.latency.in_ns()
            print(f"  Overhead: {overhead:.2f}x")
        
        fidelity_loss = logical_result.fidelity - arch_result.fidelity
        print(f"  Fidelity loss: {fidelity_loss:.1%}")
    
    print("\n=== OPTION 2: Use relative comparison ===")
    
    # For each error rate, compare against ideal baseline
    baseline_vm = VMFactory.create_star_vm(0.0)
    baseline_result = baseline_vm.analyze(circuit)
    
    print(f"\nBaseline STAR (ε=0.0): latency={baseline_result.latency}, fidelity={baseline_result.fidelity:.4f}")
    
    for error_rate in [1e-4, 1e-3, 1e-2]:
        error_vm = VMFactory.create_star_vm(error_rate)
        error_result = error_vm.analyze(circuit)
        
        if baseline_result.latency and error_result.latency:
            overhead = error_result.latency.in_ns() / baseline_result.latency.in_ns()
        else:
            overhead = 1.0
            
        fidelity_loss = baseline_result.fidelity - error_result.fidelity
        
        print(f"Error ε={error_rate}: overhead={overhead:.2f}x, fidelity_loss={fidelity_loss:.1%}")
    
    print("\n=== RECOMMENDATION ===")
    print("Fix approach: For consistency, always compare STAR VMs:")
    print("- Logical level: STAR VM with ε=0.0 (ideal STAR)")  
    print("- Arch level: STAR VM with specified error rate")
    print("- This gives meaningful overhead and fidelity comparisons")

if __name__ == "__main__":
    analyze_fix_options()