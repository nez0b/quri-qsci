# Comprehensive TE-QSCI Testing Suite

This directory contains extensive testing for TE-QSCI (Time-Evolved Quantum Selected Configuration Interaction) algorithms with exact diagonalization verification.

## Answer to User's Question

**Question**: "In 'tests/test_qsci_verification.py', you tested various Hamiltonian on VanillaQSCI. Did you input HF state as ψ_init and sample from ⟨ψ_init | H_test | ψ_init⟩?"

**Answer**: **NO**, the existing VanillaQSCI tests do **NOT** use Hartree-Fock (HF) states as initial states. Instead, they use `create_diverse_superposition_state()` which creates engineered superposition states with:

- RY rotations with varying angles per qubit: `θ * (i + 1) / n_qubits`
- CNOT gates for entanglement between adjacent qubits
- RZ rotations for additional phase complexity: `θ / 4`

This choice ensures **diverse measurement outcomes** for robust QSCI testing, rather than using potentially localized HF states.

## Test Organization

```
/tests/exact_diagonalizations/
├── conftest.py                        # Shared fixtures and configurations
├── test_te_qsci_single_time.py       # Single-time TE-QSCI verification tests
├── test_te_qsci_time_average.py      # Time-average TE-QSCI verification tests
├── test_te_qsci_state_vector.py      # State vector input TE-QSCI tests
├── test_initial_state_comparison.py   # HF vs superposition vs arbitrary comparison
└── README.md                          # This documentation
```

## Test Methodology

### 1. Algorithm Coverage

- **Single-Time TE-QSCI**: Tests evolution at specific time points
- **Time-Average TE-QSCI**: Tests time-averaged measurements over multiple evolution periods
- **State Vector TE-QSCI**: Tests direct state vector input variants
- **VanillaQSCI**: Baseline comparison using existing methodology

### 2. Hamiltonian Test Suite

All tests use the same comprehensive Hamiltonian suite for consistency:

- **TFIM (Transverse Field Ising Model)**: Various field strengths (0.5, 1.0, 2.0) across quantum phase transitions
- **Heisenberg Model**: Different coupling ratios (isotropic, anisotropic, XX model)
- **Random Sparse**: Controlled condition numbers and sparsity patterns
- **Simple Test Cases**: 2x2 matrices and Pauli Hamiltonians for exact validation

### 3. Verification Framework

- **Exact Solutions**: Using `scipy.sparse.linalg.eigsh` from `hamiltonian_factory` utilities
- **Realistic Tolerances**: 1e-2 to 0.5 for QSCI approximation methods (not exact method tolerances)
- **System Sizes**: 2-3 qubits optimal for exact diagonalization verification
- **Statistical Analysis**: Multiple trials and convergence studies

## Key Findings

### Initial State Performance Comparison

Based on comprehensive testing across multiple Hamiltonian types:

| Initial State Type | Average Error | Measurement Diversity | Recommendation |
|-------------------|---------------|----------------------|----------------|
| **Diverse Superposition** | Best-Good | **High** | ✅ **Recommended** |
| Hartree-Fock | Good | Low-Medium | ⚠️ Limited diversity |
| Arbitrary Random | Variable | Medium-High | 🔄 Inconsistent |
| Computational Basis | Poor | **Low** | ❌ Not recommended |

### Why Superposition States Are Preferred

1. **Measurement Diversity**: Higher Shannon entropy and effective dimension
2. **Robust Testing**: Consistent performance across different Hamiltonian types
3. **QSCI Convergence**: Comparable or better convergence rates
4. **Statistical Sampling**: More diverse computational basis outcomes for robust QSCI subspace selection

## Running the Tests

### Individual Test Modules

```bash
# Test specific TE-QSCI variants
pytest tests/exact_diagonalizations/test_te_qsci_single_time.py -v
pytest tests/exact_diagonalizations/test_te_qsci_time_average.py -v
pytest tests/exact_diagonalizations/test_te_qsci_state_vector.py -v

# Test initial state comparison (answers user's question)
pytest tests/exact_diagonalizations/test_initial_state_comparison.py -v
```

### Test by Markers

```bash
# Test all TE-QSCI variants
pytest -m "te_qsci_single_time" -v
pytest -m "te_qsci_time_average" -v
pytest -m "te_qsci_state_vector" -v

# Test initial state comparisons
pytest -m "initial_state_comparison" -v

# Test all exact diagonalization verification
pytest -m "exact_diagonalization" -v
```

### Complete Test Suite

```bash
# Run all TE-QSCI tests
pytest tests/exact_diagonalizations/ -v
```

## Test Parameters

### Evolution Time Parameters
- **Geometric Progression**: [0.1, 0.3, 1.0, 3.0] atomic units
- **Coverage**: Short-time, intermediate, and long-time evolution regimes

### Shot Count Analysis
- **Convergence Study**: [500, 1000, 2000] shots
- **Statistical Validation**: Multiple trials for variance analysis

### Tolerance Settings
- **Exact Methods**: 1e-6 (for reference solutions)
- **QSCI Verification**: 1e-2 (realistic for approximation methods)
- **Loose Tolerance**: 0.5 (for difficult cases and random Hamiltonians)

## Implementation Notes

### Shared Utilities

- **Hamiltonian Factory**: Reuses existing utilities from `../utils/hamiltonian_factory.py`
- **State Creation**: Factories for HF, superposition, and arbitrary states in `conftest.py`
- **Pauli Decomposition**: Matrix to QURI operator conversion for small systems
- **Exact Solutions**: Consistent exact diagonalization using sparse methods

### Test Architecture

- **Parametrized Tests**: Systematic coverage of parameter spaces
- **Fixture-Based**: Shared configurations and tolerance settings
- **Modular Design**: Each TE-QSCI variant in separate test files
- **Cross-Validation**: Comparison tests between algorithm variants

## Validation Results

The comprehensive testing validates that:

1. **Existing Methodology is Sound**: The choice of diverse superposition states in `test_qsci_verification.py` is empirically justified
2. **TE-QSCI Variants Work**: All TE-QSCI algorithms pass verification against exact diagonalization
3. **Initial State Strategy**: Superposition states provide optimal balance of accuracy and measurement diversity
4. **Method Robustness**: Results are consistent across different QSCI algorithm variants

## Future Extensions

This testing framework can be extended for:

- **Larger System Sizes**: Using approximate validation methods
- **Noise Models**: Adding quantum error effects
- **Performance Benchmarking**: Systematic timing and scaling studies
- **Parameter Optimization**: Automated hyperparameter tuning

---

**Conclusion**: The existing VanillaQSCI tests use diverse superposition states (NOT Hartree-Fock states) as initial states, and this choice is empirically validated by our comprehensive testing to provide optimal measurement diversity and robust QSCI convergence across multiple Hamiltonian types.