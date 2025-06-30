# Molecular Systems with VM-Aware QSCI

The VM-Aware QSCI module is designed to handle realistic quantum chemistry problems, providing authentic resource estimation for molecular electronic structure calculations on quantum hardware.

## Overview

Molecular quantum chemistry represents one of the most promising applications of quantum computing. The VM-Aware QSCI module enables realistic assessment of quantum resources required for molecular electronic structure problems, bridging theoretical quantum algorithms with practical quantum hardware constraints.

## Supported Molecular Systems

### Linear Hydrogen Chains

The primary focus is on linear hydrogen chain molecules, which provide scalable test systems for quantum algorithms:

#### H2 (Hydrogen Dimer)
- **Geometry**: H 0 0 0; H 0 0 1
- **Basis**: STO-3G  
- **Qubits**: 4
- **Hamiltonian Terms**: ~15
- **Typical Circuit**: 48 gates, depth 24

#### H4 (Hydrogen Tetramer)  
- **Geometry**: H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3
- **Basis**: STO-3G
- **Qubits**: 8  
- **Hamiltonian Terms**: ~185
- **Typical Circuit**: 924 gates, depth 814

#### H6 (Hydrogen Hexamer)
- **Geometry**: H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3; H 0 0 4; H 0 0 5
- **Basis**: STO-3G
- **Qubits**: 12
- **Hamiltonian Terms**: ~630  
- **Typical Circuit**: 3000+ gates, depth 2500+

## Practical Examples

### H2 System Analysis

Complete example for H2 molecular system:

```python
from quri_parts.chem.mol import make_h2_molecule
from quri_parts.algo.utils.operator import create_molecular_hamiltonian
from src.qsci_vm_interface import create_vm_aware_te_qsci
from src.qsci_vm_analysis import create_star_vm_for_h6

def analyze_h2_system(error_rate=1e-3):
    """Complete H2 system analysis with VM-aware QSCI"""
    
    # 1. Create H2 molecule
    h2_molecule = make_h2_molecule(distance=1.0)  # 1 Bohr separation
    
    # 2. Generate molecular Hamiltonian  
    hamiltonian = create_molecular_hamiltonian(h2_molecule)
    print(f"H2 Hamiltonian: {len(hamiltonian.mapping)} terms")
    
    # 3. Configure VM interface
    vm_interface = create_star_vm_for_h6(error_rate)
    
    # 4. Create VM-aware algorithm
    algorithm = create_vm_aware_te_qsci(
        hamiltonian=hamiltonian,
        evolution_time=1.0,
        trotter_steps=4,
        vm_interface=vm_interface,
        use_vm_sampling=True
    )
    
    # 5. Run calculation
    result = algorithm.run()
    
    # 6. Analyze results
    print("=== H2 System Results ===")
    print(f"Ground state energy: {result.qsci_result.ground_state_energy:.6f} Ha")
    print(f"Selected energies: {len(result.qsci_result.selected_energies)} states")
    
    if result.vm_analysis:
        logical = result.vm_analysis.logical_analysis
        arch = result.vm_analysis.arch_analysis
        print(f"\nCircuit Resources:")
        print(f"  Logical: {logical.gate_count} gates, depth {logical.circuit_depth}")
        print(f"  STAR: {arch.gate_count} gates, depth {arch.circuit_depth}")
        print(f"  Overhead: {arch.gate_count / logical.gate_count:.1f}x")
    
    return result

# Run H2 analysis
h2_result = analyze_h2_system()
```

### H4 System for Development

H4 provides a good balance between realism and computational efficiency:

```python
def create_h4_hamiltonian():
    """Create H4 linear chain Hamiltonian"""
    from quri_parts.core.operator import Operator, pauli_label
    
    # Simplified H4 Hamiltonian (in practice, use quantum chemistry package)
    # This is a representative example - actual H4 would have ~185 terms
    hamiltonian = Operator({
        # Single-qubit terms
        pauli_label("Z0"): -0.8105,
        pauli_label("Z1"): -0.8105, 
        pauli_label("Z2"): -0.8105,
        pauli_label("Z3"): -0.8105,
        pauli_label("Z4"): -0.8105,
        pauli_label("Z5"): -0.8105,
        pauli_label("Z6"): -0.8105,
        pauli_label("Z7"): -0.8105,
        
        # Two-qubit interactions (representative subset)
        pauli_label("Z0 Z1"): -0.3123,
        pauli_label("Z2 Z3"): -0.3123,
        pauli_label("Z4 Z5"): -0.3123,
        pauli_label("Z6 Z7"): -0.3123,
        
        # Inter-fragment interactions
        pauli_label("Z0 Z2"): -0.0987,
        pauli_label("Z1 Z3"): -0.0987,
        
        # Exchange terms
        pauli_label("X0 X1"): -0.1812,
        pauli_label("Y0 Y1"): -0.1812,
        pauli_label("X2 X3"): -0.1812,
        pauli_label("Y2 Y3"): -0.1812,
        
        # Additional molecular interactions...
        # (In practice, generated from quantum chemistry calculation)
    })
    
    return hamiltonian

def analyze_h4_scaling():
    """Analyze H4 system resource scaling"""
    hamiltonian = create_h4_hamiltonian()
    error_rates = [0.0, 1e-5, 1e-4, 1e-3]
    
    results = {}
    for error_rate in error_rates:
        vm_interface = create_star_vm_for_h6(error_rate)
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=1.0,
            vm_interface=vm_interface
        )
        
        result = algorithm.run()
        results[error_rate] = result
    
    # Compare results across error rates
    print("=== H4 Error Rate Impact ===")
    for error_rate, result in results.items():
        energy = result.qsci_result.ground_state_energy
        print(f"Error rate {error_rate}: Energy = {energy:.6f} Ha")
        
        if result.vm_analysis:
            arch = result.vm_analysis.arch_analysis
            if hasattr(arch, 'circuit_fidelity') and arch.circuit_fidelity:
                print(f"  Circuit fidelity: {arch.circuit_fidelity:.4f}")
    
    return results
```

### H6 System for Research

H6 represents a computationally challenging molecular system suitable for research applications:

```python
def analyze_h6_performance():
    """Comprehensive H6 system performance analysis"""
    
    # Note: H6 Hamiltonian would typically be generated from quantum chemistry
    # This example shows the analysis structure
    
    # Different evolution times to study time dependence
    evolution_times = [0.5, 1.0, 2.0, 4.0]
    error_rate = 1e-3
    
    vm_interface = create_star_vm_for_h6(error_rate)
    
    performance_data = []
    
    for t in evolution_times:
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=h6_hamiltonian,  # Assume pre-generated
            evolution_time=t,
            trotter_steps=4,
            vm_interface=vm_interface
        )
        
        # Time the calculation
        import time
        start_time = time.time()
        result = algorithm.run()
        end_time = time.time()
        
        data = {
            'evolution_time': t,
            'calculation_time': end_time - start_time,
            'energy': result.qsci_result.ground_state_energy,
        }
        
        if result.vm_analysis:
            arch = result.vm_analysis.arch_analysis
            data.update({
                'gates': arch.gate_count,
                'depth': arch.circuit_depth,
                'fidelity': getattr(arch, 'circuit_fidelity', None)
            })
        
        performance_data.append(data)
    
    # Analysis results
    print("=== H6 Performance Analysis ===")
    for data in performance_data:
        print(f"t={data['evolution_time']}: {data['gates']} gates, "
              f"fidelity={data.get('fidelity', 'N/A')}")
    
    return performance_data
```

## Resource Scaling Studies

### Cross-System Comparison

Compare resource requirements across different molecular systems:

```python
def compare_molecular_systems():
    """Compare H2, H4, H6 resource requirements"""
    
    systems = {
        'H2': {'qubits': 4, 'hamiltonian': h2_hamiltonian},
        'H4': {'qubits': 8, 'hamiltonian': h4_hamiltonian}, 
        'H6': {'qubits': 12, 'hamiltonian': h6_hamiltonian}
    }
    
    error_rate = 1e-3
    vm_interface = create_star_vm_for_h6(error_rate)
    
    comparison_data = {}
    
    for system_name, system_data in systems.items():
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=system_data['hamiltonian'],
            evolution_time=1.0,
            vm_interface=vm_interface
        )
        
        # Generate circuit for analysis
        circuit = algorithm._generate_authentic_time_evolution_circuit(
            system_data['qubits']
        )
        
        # Analyze at both levels
        logical = vm_interface.analyze_circuit_at_level(circuit, "LogicalCircuit")
        arch = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")
        
        comparison_data[system_name] = {
            'qubits': system_data['qubits'],
            'hamiltonian_terms': len(system_data['hamiltonian'].mapping),
            'logical_gates': logical.gate_count,
            'logical_depth': logical.circuit_depth,
            'arch_gates': arch.gate_count,
            'arch_depth': arch.circuit_depth,
            'overhead': arch.gate_count / logical.gate_count
        }
    
    # Display comparison
    print("=== Molecular System Comparison ===")
    print(f"{'System':<6} {'Qubits':<7} {'H-terms':<8} {'Log Gates':<10} {'Arch Gates':<11} {'Overhead':<9}")
    print("-" * 70)
    
    for system, data in comparison_data.items():
        print(f"{system:<6} {data['qubits']:<7} {data['hamiltonian_terms']:<8} "
              f"{data['logical_gates']:<10} {data['arch_gates']:<11} "
              f"{data['overhead']:<9.1f}x")
    
    return comparison_data
```

### Time Evolution Parameter Studies

Study how algorithm parameters affect molecular calculations:

```python
def study_trotter_parameters(hamiltonian, n_qubits):
    """Study Trotter step impact on molecular calculations"""
    
    trotter_steps = [2, 4, 8, 16]
    evolution_time = 1.0
    error_rate = 1e-3
    
    vm_interface = create_star_vm_for_h6(error_rate)
    results = {}
    
    for steps in trotter_steps:
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=evolution_time,
            trotter_steps=steps,
            vm_interface=vm_interface
        )
        
        result = algorithm.run()
        
        # Generate circuit for resource analysis
        circuit = algorithm._generate_authentic_time_evolution_circuit(n_qubits)
        arch_analysis = vm_interface.analyze_circuit_at_level(circuit, "ArchLogicalCircuit")
        
        results[steps] = {
            'energy': result.qsci_result.ground_state_energy,
            'gates': arch_analysis.gate_count,
            'depth': arch_analysis.circuit_depth,
            'fidelity': getattr(arch_analysis, 'circuit_fidelity', None)
        }
    
    print("=== Trotter Steps Impact ===")
    print(f"{'Steps':<6} {'Energy':<12} {'Gates':<8} {'Depth':<8} {'Fidelity':<10}")
    print("-" * 50)
    
    for steps, data in results.items():
        fidelity_str = f"{data['fidelity']:.4f}" if data['fidelity'] else "N/A"
        print(f"{steps:<6} {data['energy']:<12.6f} {data['gates']:<8} "
              f"{data['depth']:<8} {fidelity_str:<10}")
    
    return results
```

## Error Analysis for Molecular Systems

### Quantum Error Impact

Study how quantum errors affect molecular energy calculations:

```python
def analyze_molecular_error_sensitivity(hamiltonian, n_qubits):
    """Analyze molecular calculation sensitivity to quantum errors"""
    
    error_rates = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    baseline_energy = None
    
    error_analysis = {}
    
    for error_rate in error_rates:
        vm_interface = create_star_vm_for_h6(error_rate)
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=1.0,
            vm_interface=vm_interface
        )
        
        result = algorithm.run()
        energy = result.qsci_result.ground_state_energy
        
        if baseline_energy is None:
            baseline_energy = energy
        
        energy_error = abs(energy - baseline_energy)
        
        error_analysis[error_rate] = {
            'energy': energy,
            'energy_error': energy_error,
            'relative_error': energy_error / abs(baseline_energy) if baseline_energy != 0 else 0
        }
    
    print("=== Molecular Error Sensitivity ===")
    print(f"{'Error Rate':<12} {'Energy':<12} {'Energy Error':<13} {'Relative Error':<15}")
    print("-" * 60)
    
    for error_rate, data in error_analysis.items():
        print(f"{error_rate:<12.1e} {data['energy']:<12.6f} "
              f"{data['energy_error']:<13.6f} {data['relative_error']:<15.2e}")
    
    return error_analysis
```

## Best Practices

### System Selection

**For Development and Testing**:
- Start with H2 (4 qubits) for rapid prototyping
- Use H4 (8 qubits) for realistic development
- Reserve H6 (12 qubits) for final validation

**For Research**:
- H6 and larger systems for realistic studies
- Consider computational resources and time constraints
- Use appropriate error rates for target quantum hardware

### Performance Optimization

```python
def optimize_molecular_calculation(hamiltonian, n_qubits, target_accuracy=1e-6):
    """Optimize calculation parameters for target accuracy"""
    
    # Start with conservative parameters
    base_params = {
        'evolution_time': 1.0,
        'trotter_steps': 4,
        'error_rate': 1e-3
    }
    
    # Test parameter variations
    param_variations = [
        {'trotter_steps': 2},   # Faster, less accurate
        {'trotter_steps': 8},   # Slower, more accurate
        {'evolution_time': 0.5}, # Shorter evolution
        {'error_rate': 1e-4},   # Lower error rate
    ]
    
    results = {}
    
    # Baseline calculation
    vm_interface = create_star_vm_for_h6(base_params['error_rate'])
    baseline_algorithm = create_vm_aware_te_qsci(
        hamiltonian=hamiltonian,
        evolution_time=base_params['evolution_time'],
        trotter_steps=base_params['trotter_steps'],
        vm_interface=vm_interface
    )
    baseline_result = baseline_algorithm.run()
    baseline_energy = baseline_result.qsci_result.ground_state_energy
    
    results['baseline'] = {
        'params': base_params,
        'energy': baseline_energy,
        'gates': 'TBD'  # Would extract from circuit analysis
    }
    
    # Test variations
    for i, variation in enumerate(param_variations):
        params = {**base_params, **variation}
        
        vm_interface = create_star_vm_for_h6(params['error_rate'])
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=params['evolution_time'],
            trotter_steps=params['trotter_steps'],
            vm_interface=vm_interface
        )
        
        result = algorithm.run()
        energy = result.qsci_result.ground_state_energy
        energy_diff = abs(energy - baseline_energy)
        
        results[f'variant_{i}'] = {
            'params': params,
            'energy': energy,
            'energy_diff': energy_diff,
            'meets_accuracy': energy_diff <= target_accuracy
        }
    
    return results
```

### Integration with Quantum Chemistry

For realistic molecular calculations, integrate with quantum chemistry packages:

```python
# Example integration structure (implementation depends on specific packages)
def create_molecular_hamiltonian_from_chemistry(molecule_spec, basis='sto-3g'):
    """Create molecular Hamiltonian from quantum chemistry calculation"""
    
    # This would typically use packages like:
    # - PySCF for quantum chemistry
    # - quri-parts molecular modules
    # - OpenFermion for operator manipulation
    
    # Placeholder for actual implementation
    pass

def run_molecular_qsci_study(molecule_specs):
    """Run comprehensive QSCI study across multiple molecules"""
    
    results = {}
    
    for mol_name, mol_spec in molecule_specs.items():
        print(f"Analyzing {mol_name}...")
        
        # Generate Hamiltonian from quantum chemistry
        hamiltonian = create_molecular_hamiltonian_from_chemistry(mol_spec)
        
        # Run VM-aware QSCI
        vm_interface = create_star_vm_for_h6(1e-3)
        algorithm = create_vm_aware_te_qsci(
            hamiltonian=hamiltonian,
            evolution_time=1.0,
            vm_interface=vm_interface
        )
        
        result = algorithm.run()
        results[mol_name] = result
    
    return results
```

## Troubleshooting

### Common Issues

**Large Circuit Generation**: 
For H6 and larger systems, circuit generation can be slow. Consider:
- Using H4 for development
- Reducing Trotter steps initially
- Implementing circuit caching

**Memory Usage**:
Large molecular Hamiltonians can consume significant memory:
- Monitor Hamiltonian term count
- Use appropriate VM configurations
- Consider system resource limitations

**Accuracy vs Performance**:
Balance calculation accuracy with computational resources:
- Start with lower Trotter steps
- Increase gradually until convergence
- Monitor energy convergence patterns

## Next Steps

- **[VM Analysis](vm_analysis.md)** - Understand resource estimation details
- **[Authentic Time Evolution](authentic_time_evolution.md)** - Deep dive into circuit generation
- **[API Reference](../../api/vm_qsci.md)** - Complete technical documentation