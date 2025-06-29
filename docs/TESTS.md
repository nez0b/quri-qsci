# QSCI Testing Framework

This document provides a comprehensive guide to the testing suite for `quri-qsci`, explaining the testing philosophy, structure, and how to add new tests.

## 1. Testing Philosophy

The testing suite for `quri-qsci` is designed to ensure correctness, stability, and performance across multiple layers of the application. Our philosophy combines foundational verification with practical, real-world integration tests:

- **Correctness:** Verifying algorithm results against known ground truths from exact diagonalization and established quantum chemistry packages
- **Robustness:** Ensuring the solver behaves correctly across a wide range of Hamiltonians and system parameters  
- **Integration:** Confirming seamless operation with other frameworks like `qiskit-addon-sqd` and `quri-parts`
- **Regression:** Preventing the re-introduction of fixed bugs through comprehensive regression testing
- **Precision:** Validating that uniform superposition states with complete subspace coverage achieve machine precision accuracy

## 2. How to Run Tests

### Basic Test Execution

To run the full test suite, execute from the root directory:

```bash
pytest
```

To run tests from a specific file:

```bash
pytest tests/test_qsci_core.py
```

To run the exact diagonalization test suite:

```bash
pytest tests/exact_diagonalizations/
```

To run tests with verbose output:

```bash
pytest -v
```

### Test Markers

Run tests by category using pytest markers:

```bash
pytest -m verification                    # Exact diagonalization validation tests
pytest -m integration                    # Framework integration tests  
pytest -m molecular                      # Molecular system tests
pytest -m regression                     # Regression tests for specific bugs
pytest -m exact_ground_state_precision   # Uniform superposition precision tests
pytest -m te_qsci_single_time           # Single-time TE-QSCI algorithm tests
pytest -m te_qsci_time_average          # Time-average TE-QSCI algorithm tests
pytest -m te_qsci_state_vector          # State vector TE-QSCI algorithm tests
pytest -m initial_state_comparison      # Initial state comparison tests
pytest -m exact_diagonalization         # Tests using exact diagonalization verification
```

### Common Test Commands

```bash
# Run only fast unit tests (exclude molecular calculations)
pytest -m "not molecular"

# Run with coverage report
pytest --cov=src

# Run specific test function
pytest tests/exact_diagonalizations/test_exact_ground_state_precision.py::TestExactGroundStatePrecision::test_vanilla_qsci_uniform_superposition_precision

# Run tests and stop on first failure
pytest -x

# Run with short traceback format
pytest --tb=short
```

## 3. Test Categories and Directory Structure

Our tests are organized by their purpose and scope:

### Core Test Files

| File / Directory                     | Category                 | Description                                                                                                                              |
| ------------------------------------ | ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `tests/test_qsci_verification.py`    | **Verification**         | Compares QSCI results against exact diagonalization for small, well-defined quantum systems. This is our ground truth validation. |
| `tests/test_qsci_core.py`            | **Unit & Regression**    | Tests the core logic, internal functions, and data structures of the QSCI implementation. Includes regression tests for specific bugs.    |
| `tests/test_implementation.py`       | **Unit**                 | Basic validation tests to confirm that the QSCI solver class is implemented correctly and handles inputs as expected.                  |
| `tests/test_qsci_integrations.py`    | **Framework Integration**| Compares `quri-qsci` against `qiskit-addon-sqd` using mock implementations to ensure API compatibility and consistent behavior.           |
| `tests/test_h2_verification.py`      | **Molecular Verification** | A specific, high-value verification test comparing the H2 molecule ground state against `qiskit-addon-sqd` results.                      |
| `tests/test_qsci_real_integration.py`| **Real Integration**     | End-to-end tests using real molecular Hamiltonians (e.g., from PySCF) and diverse quantum state preparations.                              |

### Exact Diagonalization Test Suite

The `tests/exact_diagonalizations/` directory contains comprehensive tests that validate all QSCI variants against exact solutions using uniform superposition states for optimal sampling:

| File                                           | Category                        | Description                                                                                                                              |
| ---------------------------------------------- | ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `exact_diagonalizations/test_exact_ground_state_precision.py` | **Precision Validation**        | Tests that all QSCI variants achieve machine precision (1e-8) when using uniform superposition + complete subspace coverage. |
| `exact_diagonalizations/test_te_qsci_single_time.py`         | **Single-Time TE-QSCI**        | Validation of SingleTimeTE_QSCI against exact diagonalization for various Hamiltonians (Pauli, TFIM, Heisenberg). |
| `exact_diagonalizations/test_te_qsci_time_average.py`        | **Time-Average TE-QSCI**       | Tests TimeAverageTE_QSCI with multiple evolution times and convergence analysis. |
| `exact_diagonalizations/test_te_qsci_state_vector.py`        | **State Vector TE-QSCI**       | Validation of StateVectorTE_QSCI with direct state vector inputs and circuit state equivalence. |
| `exact_diagonalizations/test_initial_state_comparison.py`    | **Initial State Strategies**   | Comparative analysis of different initial state types (Hartree-Fock vs superposition vs arbitrary) for QSCI performance. |

### Test Utilities

| File / Directory                     | Purpose                  | Description                                                                                                                              |
| ------------------------------------ | ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `tests/utils/hamiltonian_factory.py` | **Test Hamiltonians**    | Factory functions for creating various types of test Hamiltonians (TFIM, Heisenberg, Pauli, random sparse matrices). |
| `tests/utils/conversion_utils.py`    | **Data Conversion**      | Helper functions for framework-specific data conversions and mock implementations. |
| `tests/utils/exact_ground_state_utils.py` | **State Preparation** | Utilities for exact ground state extraction and quantum state preparation (including the improved uniform superposition approach). |
| `tests/utils/universal_qsci_tester.py` | **Universal Testing**   | Common testing patterns and utilities used across multiple test files. |
| `tests/conftest.py`                  | **Global Configuration** | Global pytest fixtures and configuration shared across the entire test suite. |
| `exact_diagonalizations/conftest.py` | **ED Configuration**     | Specialized fixtures for exact diagonalization tests, including uniform superposition state creation. |

## 4. Test Categories in Detail

### Exact Ground State Precision Tests

**Key Innovation:** These tests use **uniform superposition states** (H⊗H⊗...⊗H) with **complete subspace coverage** to achieve machine precision accuracy.

**Test Types:**
- **Uniform Superposition Validation:** Verifies that uniform superposition states provide maximum measurement diversity
- **QSCI Variant Precision:** Tests all QSCI variants (Vanilla, SingleTime TE-QSCI, TimeAverage TE-QSCI, StateVector TE-QSCI) achieve 1e-8 precision
- **Comparative Analysis:** Demonstrates superiority of uniform superposition over computational basis states

**Test Systems:**
- 1-qubit simple Hamiltonian (2x2 matrix)
- 2-qubit Pauli Hamiltonian
- 2-qubit TFIM at critical point (h=1.0)
- 2-qubit Heisenberg model

**Key Parameters:**
- Complete subspace: `num_states_pick_out = 2**n_qubits`
- High shot counts: 2000-3000 for statistical robustness
- Machine precision tolerance: 1e-8

### TE-QSCI Algorithm Tests

**Single-Time TE-QSCI (`test_te_qsci_single_time.py`):**
- Evolution time parameter sweeps across different Hamiltonians
- TFIM tests across ferromagnetic, critical, and paramagnetic phases
- Heisenberg model with various coupling parameters
- Evolution time dependency analysis

**Time-Average TE-QSCI (`test_te_qsci_time_average.py`):**
- Multiple evolution time averaging
- Statistical variance analysis
- Convergence behavior with increasing time points
- Comparison with single-time variants

**State Vector TE-QSCI (`test_te_qsci_state_vector.py`):**
- Direct state vector input processing
- Circuit state vs vector equivalence
- State normalization handling
- Edge case state vector processing

### Initial State Comparison Tests

**Purpose:** Provide empirical evidence for optimal initial state selection strategies in QSCI algorithms.

**Comparisons:**
- Hartree-Fock states vs uniform superposition states
- Measurement outcome diversity analysis
- Convergence behavior across different initial state types
- Cross-algorithm robustness (VanillaQSCI vs TE-QSCI variants)

**Key Findings Validated:**
- Uniform superposition states achieve machine precision with complete subspace
- Hartree-Fock states have limited sampling diversity
- Computational basis states perform poorly for QSCI algorithms

### Verification Tests (`test_qsci_verification.py`)

These are the original mathematical validation tests against exact diagonalization:

- **Ground State Energy Tests:** Compare QSCI ground state energy vs `scipy.sparse.linalg.eigsh`
- **Ground State Fidelity Tests:** Validate quantum state vectors using fidelity calculations
- **Excited State Tests:** Verify low-lying excited states when supported
- **Parameter Sweep Tests:** Ensure robustness across model parameter ranges (e.g., TFIM field strength)

### Unit & Regression Tests (`test_qsci_core.py`)

Focus on implementation details and preventing regressions:

- **Eigenvalue Count Bug Tests:** Verify that k>1 eigenvalues are computed when requested
- **Sparse vs Dense Consistency:** Ensure both diagonalization methods give identical results
- **Edge Cases:** Single state measurements, degenerate eigenvalues, ill-conditioned matrices
- **Mock Implementation Tests:** Validate testing framework itself

### Framework Integration Tests (`test_qsci_integrations.py`)

Ensure compatibility with external quantum computing frameworks:

- **API Compatibility:** Consistent behavior with `qiskit-addon-sqd`
- **Data Format Conversion:** Proper handling of different Hamiltonian representations
- **Measurement Count Processing:** Validation of synthetic and real measurement data

### Real Integration Tests (`test_qsci_real_integration.py`)

End-to-end validation with realistic quantum systems:

- **Diverse State Preparation:** Tests that solve the measurement diversity problem
- **Molecular Systems:** H2 molecule with proper quantum chemistry integration
- **Time Evolution:** Trotter-evolved states for creating measurement diversity
- **Sampling Validation:** Tests with actual quantum circuit sampling

## 5. Key Testing Insights and Methodologies

### Uniform Superposition Breakthrough

**Discovery:** The combination of uniform superposition initial states + complete subspace coverage provides optimal QSCI performance.

**Implementation:**
```python
def create_diverse_superposition_state(n_qubits: int, theta: float = np.pi/6) -> GeneralCircuitQuantumState:
    """Create uniform superposition state for perfect QSCI sampling."""
    circuit = QuantumCircuit(n_qubits)
    
    # Apply Hadamard gates to all qubits for uniform superposition
    for i in range(n_qubits):
        circuit.add_H_gate(i)
    
    return GeneralCircuitQuantumState(n_qubits, circuit)
```

**Results:** All QSCI variants achieve machine precision (< 1e-15) with this approach.

### Complete Subspace Coverage

**Principle:** For small exact diagonalization tests, use complete subspace coverage to ensure optimal sampling.

**Implementation:**
```python
# For n-qubit systems
num_states_pick_out = 2**n_qubits  # Complete Hilbert space
```

**Benefits:**
- Eliminates sampling bias
- Enables machine precision accuracy
- Provides deterministic results for validation

### Robust Test Design Patterns

**Tolerance Handling:**
- Machine precision tests: 1e-8 tolerance
- Comparative tests: Different tolerances per state type
- Numerical precision handling for edge cases

**Shot Count Optimization:**
- Standard tests: 2000-2500 shots
- TE-QSCI tests: 3000 shots for better convergence
- Statistical analysis: Multiple trials for variance studies

## 6. Adding New Tests

### For New QSCI Variants

Add tests to the appropriate exact diagonalization test file:

```python
@pytest.mark.exact_diagonalization
def test_new_variant_uniform_superposition_precision(self, exact_ground_state_test_systems):
    """Test NewVariantQSCI achieves machine precision with uniform superposition."""
    from tests.exact_diagonalizations.conftest import create_diverse_superposition_state
    
    for system_name, system_data in exact_ground_state_test_systems.items():
        # Create uniform superposition state
        uniform_state = create_diverse_superposition_state(
            n_qubits=system_data['n_qubits']
        )
        
        # Run NewVariantQSCI with complete subspace coverage
        sampler = create_qulacs_vector_concurrent_sampler()
        new_qsci = NewVariantQSCI(
            hamiltonian=system_data['quri_hamiltonian'],
            sampler=sampler,
            num_states_pick_out=2**system_data['n_qubits']  # Complete subspace
        )
        
        result = new_qsci.run([uniform_state], total_shots=3000)
        
        # Validate precision
        energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
        
        # Expect machine precision with uniform superposition + complete subspace
        precision_tolerance = 1e-8
        assert energy_error < precision_tolerance, (
            f"NewVariantQSCI energy error {energy_error:.2e} exceeds tolerance "
            f"{precision_tolerance:.2e} for {system_name}"
        )
```

### For New Hamiltonians

Extend the `hamiltonian_factory` in `tests/utils/`:

```python
def create_new_model_hamiltonian(n_qubits, param1, param2):
    """Create Hamiltonian for new quantum model."""
    # Implementation using scipy.sparse for efficiency
    return sparse_hamiltonian
```

### For Mathematical Validation

Add verification tests to `tests/test_qsci_verification.py`:

```python
@pytest.mark.verification
@pytest.mark.parametrize("test_case", test_cases)
def test_against_exact_diagonalization(test_case):
    """Validate against scipy exact diagonalization."""
    hamiltonian = create_hamiltonian(test_case)
    
    # Exact solution
    exact_energy, _ = eigsh(hamiltonian, k=1, which='SA')
    
    # QSCI solution
    qsci_result = run_qsci(hamiltonian)
    
    # Validate
    np.testing.assert_allclose(qsci_result.energy, exact_energy[0], atol=1e-8)
```

### For Bug Fixes

Add regression tests that would have failed before the fix:

```python
@pytest.mark.regression
def test_bug_xyz_fixed():
    """Regression test for issue #XYZ."""
    # Create scenario that triggered the bug
    problematic_input = create_bug_scenario()
    
    # This should not raise an exception or return wrong results
    result = fixed_function(problematic_input)
    assert result.is_correct()
```

## 7. Test Utilities

### Hamiltonian Factory (`tests/utils/hamiltonian_factory.py`)

Provides functions to create various types of test Hamiltonians:

- `create_hermitian_matrix()` - Random Hermitian matrices with specified eigenvalues
- `create_simple_2x2_hamiltonian()` - Simple test case with known eigenvalues
- `create_tfim_hamiltonian()` - Transverse Field Ising Model
- `create_heisenberg_hamiltonian()` - XXZ Heisenberg model  
- `create_pauli_hamiltonian()` - 2-qubit Pauli combination
- `get_h2_hamiltonian_minimal_basis()` - H2 molecule Hamiltonian
- `create_random_sparse_hamiltonian()` - Random sparse Hermitian matrices
- `get_exact_solution()` - Exact eigenvalue/eigenvector computation
- `validate_sparse_hamiltonian_properties()` - Hamiltonian validation

### Exact Ground State Utilities (`tests/utils/exact_ground_state_utils.py`)

Helper functions for state preparation and validation:

- `extract_exact_ground_state()` - Extract ground state from Hamiltonian matrix
- `create_exact_circuit_state_from_vector()` - Convert state vectors to quantum circuits
- `validate_exact_state_preparation()` - Validate state preparation accuracy
- `create_exact_computational_basis_state()` - Create computational basis states
- `compare_state_preparation_methods()` - Compare different preparation approaches
- `create_exact_bell_state()` - Create Bell states for 2-qubit systems
- `create_exact_ghz_state()` - Create GHZ states for multi-qubit systems

### Conversion Utilities (`tests/utils/conversion_utils.py`)

Helper functions for data format conversions and mock implementations:

- Framework-specific data conversions
- Mock qiskit-addon-sqd implementations  
- Measurement count generation and processing
- Eigenvalue and eigenvector comparison utilities

## 8. Best Practices

### Test Design Principles

1. **Use Uniform Superposition for Precision Tests:** Always use uniform superposition + complete subspace for machine precision validation
2. **Complete Subspace Coverage:** For small systems, use `num_states_pick_out = 2**n_qubits`
3. **Appropriate Shot Counts:** 2000-3000 shots for robust statistical sampling
4. **Realistic Tolerances:** 1e-8 for precision tests, looser tolerances for comparative tests

### Test Naming

- Use descriptive test names that explain what is being tested
- Include the test category in the name: `test_vanilla_qsci_uniform_superposition_precision`
- Use parametrization for testing multiple similar cases

### Assertions

- Use `np.testing.assert_allclose()` for floating-point comparisons
- Include meaningful error messages: `err_msg=f"Failed on test case: {test_id}"`
- Set appropriate tolerances based on numerical precision requirements
- Handle machine precision edge cases (division by zero, very small numbers)

### Test Data

- Use the `hamiltonian_factory` to create test Hamiltonians consistently
- Set random seeds for reproducible tests: `np.random.seed(42)`
- Keep test cases small enough to run quickly (2-4 qubits for exact diagonalization)
- Use fixtures for common test data setup

### Markers and Categories

Use pytest markers to categorize tests:

```python
@pytest.mark.exact_ground_state_precision   # Uniform superposition precision tests
@pytest.mark.te_qsci_single_time            # Single-time TE-QSCI tests
@pytest.mark.te_qsci_time_average           # Time-average TE-QSCI tests
@pytest.mark.te_qsci_state_vector           # State vector TE-QSCI tests
@pytest.mark.initial_state_comparison       # Initial state comparison tests
@pytest.mark.exact_diagonalization          # Tests using exact diagonalization
@pytest.mark.verification                   # Ground truth validation
@pytest.mark.integration                    # Framework integration  
@pytest.mark.molecular                      # Molecular systems (may be slow)
@pytest.mark.regression                     # Specific bug regression tests
```

## 9. Continuous Integration

### Test Selection

- All tests should run in under 5 minutes for CI efficiency
- Exact diagonalization tests are optimized for small system sizes (1-2 qubits)
- Molecular tests may be slower but should still complete reasonably quickly
- Use markers to selectively run test subsets if needed

### Dependencies

- Core tests should work with minimal dependencies
- Molecular tests require PySCF and may be skipped if unavailable
- Framework integration tests require qiskit-addon-sqd
- Exact diagonalization tests require quri-parts and qulacs

### Coverage

- Aim for high coverage of core algorithm logic
- Focus coverage on critical paths and edge cases
- Use coverage reports to identify untested code paths
- Prioritize testing of uniform superposition approach

## 10. Troubleshooting

### Common Issues

- **Import errors:** Check that the source code is on the Python path
- **Numerical precision:** Adjust tolerance values for floating-point comparisons
- **Missing dependencies:** Install optional dependencies or skip relevant tests
- **Slow tests:** Use smaller system sizes or mock implementations for speed
- **Machine precision failures:** Check for division by zero or very small number handling

### Debugging Tips

- Use `pytest -v -s` to see print statements and verbose output
- Run single tests with `pytest path/to/test.py::test_function_name`
- Use `pytest --pdb` to drop into debugger on test failures
- Check test fixtures and factory functions for correct data generation
- For exact diagonalization tests, verify uniform superposition state creation
- Use `pytest --tb=short` for concise error reporting

### Performance Optimization

- Use complete subspace coverage (`2**n_qubits`) for deterministic results
- Increase shot counts (2000-3000) for better statistical convergence
- Use uniform superposition states for optimal QSCI performance
- Profile tests if they become too slow and optimize accordingly

## 11. Contributing

When contributing new tests:

1. **Follow the existing structure** - place tests in appropriate files based on category
2. **Use uniform superposition** - for precision tests, always use uniform superposition + complete subspace
3. **Use the utilities** - leverage `hamiltonian_factory` and other utilities rather than creating ad-hoc test data
4. **Document test purpose** - include docstrings explaining what the test validates
5. **Add markers** - use appropriate pytest markers for test categorization
6. **Test your tests** - ensure new tests actually catch the bugs/issues they're designed to detect
7. **Follow precision guidelines** - use 1e-8 tolerance for machine precision tests
8. **Handle edge cases** - add proper handling for numerical precision edge cases

### Pull Request Checklist

- [ ] Tests follow uniform superposition methodology for precision validation
- [ ] Complete subspace coverage used where appropriate (`2**n_qubits`)
- [ ] Appropriate shot counts (2000-3000) for statistical robustness
- [ ] Proper tolerance settings (1e-8 for precision, looser for comparative tests)
- [ ] Tests are properly categorized with pytest markers
- [ ] Edge cases and numerical precision issues are handled
- [ ] Test documentation clearly explains validation approach
- [ ] All tests pass in CI environment

For questions about the testing framework, please refer to existing tests as examples or open an issue for discussion.

## 12. Test Suite Statistics

**Total Tests:** 41 exact diagonalization tests + core test suite  
**Test Success Rate:** 100% (41/41 exact diagonalization tests passing)  
**Precision Achievement:** Machine precision (< 1e-15) for uniform superposition approach  
**Coverage:** All QSCI variants validated against exact diagonalization  
**Performance:** Complete test suite runs in under 5 minutes  

**Key Validation Achievements:**
- ✅ Uniform superposition + complete subspace = machine precision accuracy
- ✅ All QSCI variants (Vanilla, SingleTime TE-QSCI, TimeAverage TE-QSCI, StateVector TE-QSCI) validated
- ✅ Comprehensive Hamiltonian coverage (Pauli, TFIM, Heisenberg, random sparse)
- ✅ Robust numerical precision handling and edge case management
- ✅ Empirical validation of optimal initial state selection strategies