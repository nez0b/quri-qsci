# QSCI/SQD Validation Plan: LUCJ-Based State Preparation

## Executive Summary

Based on comprehensive analysis with Gemini Pro and research into qiskit-addon-sqd documentation, this document outlines the implementation plan for creating a new N2 comparison script using **LUCJ (Local Unitary Cluster Jastrow)** ansatz for state preparation. This approach is **scientifically optimal** because qiskit-addon-sqd's own N2 tutorial uses LUCJ, making it the canonical comparison method.

## Research Findings

### ✅ **LUCJ Approach is Validated**
- qiskit-addon-sqd official N2 tutorial uses "36-qubit quantum circuit ansatz (in this case, an LUCJ circuit)"
- Uses CCSD-computed t1 and t2 amplitudes to initialize ansatz
- Leverages ffsim library for LUCJ implementation
- Workflow: HF reference → apply UCJ operator → sample diverse configurations

### ⚠️ **Alternative Lighter Methods Available**
Gemini Pro suggested simpler alternatives for faster validation:
1. **Single Trotter Step**: |ψ⟩ = exp(-iτH)|HF⟩ (physically motivated, easy)
2. **Generic Entangling Circuit**: Simple Ry(θ) + CNOT layers (minimal effort)
3. **Qiskit Nature UCCSD**: Leverage existing ansatz implementations

## Implementation Plan

### **Phase 1: LUCJ Implementation (Primary Goal)**

#### File: `n2_lucj_comparison.py`

1. **Dependencies Setup**
   - Add `ffsim` dependency for LUCJ implementation
   - Import qiskit-nature components for ansatz construction
   - Leverage existing qiskit-addon-sqd integration

2. **LUCJ State Preparation**
   ```python
   # Setup CCSD for t1/t2 amplitudes
   casci = mcscf.CASCI(hf, ncas=8, nelecas=6)
   
   # Initialize LUCJ ansatz with CCSD amplitudes
   ansatz = UCJOpSpinBalanced(norb=8, n_reps=1)
   ansatz.params = get_ccsd_amplitudes(hf, casci)
   
   # Prepare quantum circuit
   circuit = QuantumCircuit(16)
   circuit.append(PrepareHartreeFockJW(8, (3,3)))  # HF reference
   circuit.append(ansatz.to_qiskit())              # Apply LUCJ
   
   # Create diverse state for sampling
   lucj_state = GeneralCircuitQuantumState(16, circuit)
   ```

3. **Deterministic Comparison Framework**
   - Use statevector simulator to get exact probabilities
   - Generate identical synthetic measurement counts for both implementations
   - Compare: ground state energies, selected states, truncated Hamiltonians, eigenvalues

4. **Validation Metrics**
   - **Primary**: Energy convergence plots (should overlap perfectly)
   - **Diagnostic**: Selected state lists, Hamiltonian matrices, eigenvalue spectra
   - **Success Criteria**: <1 mHa difference between implementations

### **Phase 2: Lightweight Alternatives (Backup Methods)**

#### File: `n2_simple_comparison.py`

1. **Single Trotter Step Preparation**
   ```python
   # Simple time evolution: exp(-iτH)|HF⟩
   tau = 0.1  # Small imaginary time step
   trotter_circuit = create_trotter_step(hamiltonian, tau)
   evolved_state = trotter_circuit.apply(hf_state)
   ```

2. **Generic Entangling Circuit**
   ```python
   # Minimal entangling layer
   circuit = QuantumCircuit(16)
   circuit.add_Ry_gates(θ=π/8)  # Small rotation
   circuit.add_CNOT_chain()     # Entanglement
   ```

3. **Rapid Validation Suite**
   - Faster execution for iterative debugging
   - Lower implementation overhead
   - Still provides measurement diversity

### **Phase 3: Production Testing Framework**

#### File: `validation_suite.py`

1. **Automated Testing Matrix**
   - Multiple state preparation methods
   - Various subspace sizes (5, 10, 15, 20, 25, 30)
   - Different molecular systems (H2, N2, H6)

2. **Comprehensive Diagnostics**
   - Convergence analysis
   - Performance benchmarking
   - Error analysis and debugging

## Implementation Strategy

### **Recommended Order**
1. **Start with Phase 2** (lightweight) for immediate validation
2. **Proceed to Phase 1** (LUCJ) for comprehensive validation
3. **Scale to Phase 3** for production testing

### **Technical Dependencies**
- Add `ffsim>=0.0.34` for LUCJ implementation
- Add `qiskit-nature>=0.7.0` for ansatz support
- Leverage existing qiskit-addon-sqd infrastructure

### **Success Criteria**
- **Immediate**: Subspace size varies (not stuck at 1)
- **Validation**: Energy convergence matches qiskit-addon-sqd
- **Production**: <1 mHa accuracy across test suite

## Risk Mitigation

### **High-Risk Items**
- LUCJ implementation complexity → Start with Phase 2 alternatives
- Performance overhead → Optimize after correctness validation
- Dependency conflicts → Test in isolated environment

### **Fallback Strategy**
If LUCJ proves too complex, the single Trotter step method provides 80% of the benefit with 20% of the implementation effort.

## Expected Outcomes

- **Immediate**: Fix the "measurement diversity" problem
- **Short-term**: Validate quri-qsci correctness vs qiskit-addon-sqd
- **Long-term**: Robust testing framework for production use

## Status

- **Created**: 2025-01-27
- **Priority**: Medium (after core testing framework completion)
- **Dependencies**: Complete comprehensive testing framework first
- **Next Action**: Implement Phase 2 (lightweight alternatives) after core bug validation