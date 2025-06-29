# TODO: Fix and Complete TE-QSCI Testing Suite

## Context

We have implemented a comprehensive TE-QSCI testing suite in `/tests/exact_diagonalizations/` but encountered import and class structure issues when attempting to run the tests. The tests were designed to validate TE-QSCI algorithms against exact diagonalization for small quantum systems.

## Current Status

### ✅ Completed
- Created complete test infrastructure in `/tests/exact_diagonalizations/`:
  - `conftest.py` - Shared fixtures and configurations
  - `test_te_qsci_single_time.py` - Single-time TE-QSCI verification tests
  - `test_te_qsci_time_average.py` - Time-average TE-QSCI verification tests
  - `test_te_qsci_state_vector.py` - State vector input TE-QSCI tests
  - `test_initial_state_comparison.py` - HF vs superposition vs arbitrary state comparisons
  - `README.md` - Comprehensive documentation

### ❌ Issues Found
1. **Import Errors**: Test files expect classes `SingleTimeTE_QSCI`, `TimeAverageTE_QSCI`, `StateVectorTE_QSCI` which don't exist
2. **Class Structure Mismatch**: Actual implementation uses `TimeEvolvedQSCI` with different run methods
3. **Missing pytest installation** (now fixed)

## Key Finding from Analysis

**User's Question Answered**: The existing VanillaQSCI tests in `/tests/test_qsci_verification.py` use `create_diverse_superposition_state()` (NOT Hartree-Fock states) as initial states. This creates engineered superposition states with RY rotations, CNOT entanglement, and RZ rotations for measurement diversity.

## Priority Tasks

### 1. Fix Class Structure Issues (HIGH PRIORITY)

**Problem**: Tests expect individual classes but implementation uses one `TimeEvolvedQSCI` class with multiple run methods.

**Solution**: Create wrapper classes or modify tests to use existing API.

**Files to modify**:
- `/src/qsci_algorithms.py` (lines 353-730) - Add wrapper classes OR
- All test files in `/tests/exact_diagonalizations/` - Update to use existing `TimeEvolvedQSCI` API

**Reference Implementation**: 
```python
# Current implementation in src/qsci_algorithms.py:
class TimeEvolvedQSCI(QSCIBase):
    def run_single_time(self, initial_state, evolution_time, total_shots, ...)
    def run_time_average(self, initial_state, evolution_times, shots_per_time, ...)
    def run_state_vector(self, initial_state, evolution_time, num_eigenstates, ...)
```

**Expected by tests**:
```python
# What tests expect:
class SingleTimeTE_QSCI(QSCIBase): ...
class TimeAverageTE_QSCI(QSCIBase): ...
class StateVectorTE_QSCI(QSCIBase): ...
```

### 2. Fix Import Issues (HIGH PRIORITY)

**Files affected**:
- `/tests/exact_diagonalizations/test_te_qsci_single_time.py` (line 25)
- `/tests/exact_diagonalizations/test_te_qsci_time_average.py` (line 25)
- `/tests/exact_diagonalizations/test_te_qsci_state_vector.py` (line 25)
- `/tests/exact_diagonalizations/test_initial_state_comparison.py` (line 29)

**Current error**: `cannot import name 'SingleTimeTE_QSCI' from 'src.qsci_algorithms'`

### 3. Run Tests and Gather Error Analysis (HIGH PRIORITY)

**Command to run**: 
```bash
cd /Users/po-jenwang/Code/QSCI/quri-qsci
source .venv/bin/activate
python -m pytest tests/exact_diagonalizations/ -v
```

**Expected behavior**: Tests should validate that QSCI algorithms find exact eigenstates for small systems with diverse input states.

**Analysis needed**: If tests fail with energy errors > tolerance, analyze why QSCI isn't finding exact solutions.

## Implementation Options

### Option A: Create Wrapper Classes (RECOMMENDED)

Add to `/src/qsci_algorithms.py`:

```python
class SingleTimeTE_QSCI(TimeEvolvedQSCI):
    def run(self, input_states, total_shots, **kwargs):
        evolution_time = kwargs.get('evolution_time', 1.0)
        return self.run_single_time(input_states[0], evolution_time, total_shots, **kwargs)

class TimeAverageTE_QSCI(TimeEvolvedQSCI):
    def __init__(self, hamiltonian, sampler, evolution_times, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_times = evolution_times
    
    def run(self, input_states, total_shots, **kwargs):
        shots_per_time = total_shots // len(self.evolution_times)
        return self.run_time_average(input_states[0], self.evolution_times, shots_per_time, **kwargs)

class StateVectorTE_QSCI(TimeEvolvedQSCI):
    def run(self, input_states, total_shots, **kwargs):
        evolution_time = kwargs.get('evolution_time', 1.0)
        return self.run_state_vector(input_states[0], evolution_time, **kwargs)
```

### Option B: Update All Test Files

Modify all test files to use `TimeEvolvedQSCI` with appropriate method calls instead of expecting separate classes.

## Deep Analysis Required (USE GEMINI PRO)

Once tests run, if energy errors exceed tolerances, analyze:

### Potential Issues:
1. **Subspace Selection**: QSCI may not be selecting the right computational basis states
2. **Hamiltonian Construction**: Truncated Hamiltonian in subspace may have errors
3. **Measurement Sampling**: Insufficient shots or poor initial state diversity
4. **Numerical Precision**: Eigenvalue computation errors
5. **Algorithm Parameters**: Wrong `num_states_pick_out` or other settings

### Analysis Questions for Gemini Pro:
1. Why would QSCI fail to find exact eigenstates on small (2-3 qubit) systems?
2. What are optimal `num_states_pick_out` values for exact recovery?
3. How does measurement diversity affect QSCI subspace quality?
4. Are there systematic biases in the state selection algorithm?

## Reference Files for Context

### Core Implementation:
- `/src/qsci_algorithms.py` - Main QSCI algorithm implementations
- `/tests/test_qsci_verification.py` - Original VanillaQSCI tests (working reference)
- `/tests/utils/hamiltonian_factory.py` - Hamiltonian generation utilities

### Test Infrastructure:
- `/tests/exact_diagonalizations/conftest.py` - Shared test fixtures
- `/tests/exact_diagonalizations/README.md` - Complete documentation of test methodology

### Key Findings:
- VanillaQSCI tests use superposition states (NOT HF states) via `create_diverse_superposition_state()`
- Existing tolerances: 1e-6 (exact), 1e-2 (QSCI verification), 0.5 (loose)
- Test Hamiltonians: TFIM, Heisenberg, Pauli, random sparse (all validated)

## Success Criteria

1. **All tests pass** with realistic tolerances for QSCI approximation
2. **Energy errors < 1e-2** for most test cases (QSCI should find near-exact solutions on small systems)
3. **Comprehensive analysis** of any failures using Gemini Pro
4. **Documentation** of optimal QSCI parameters for exact recovery

## Commands for Quick Reference

```bash
# Activate environment
cd /Users/po-jenwang/Code/QSCI/quri-qsci
source .venv/bin/activate

# Run specific test modules
python -m pytest tests/exact_diagonalizations/test_te_qsci_single_time.py -v
python -m pytest tests/exact_diagonalizations/test_initial_state_comparison.py -v

# Run by test markers
python -m pytest -m "te_qsci_single_time" -v
python -m pytest -m "initial_state_comparison" -v

# Run all exact diagonalization tests
python -m pytest tests/exact_diagonalizations/ -v

# Check what's available in algorithms
python -c "from src.qsci_algorithms import *; print([x for x in dir() if 'QSCI' in x])"
```

This TODO provides complete context for fixing and analyzing the TE-QSCI testing suite.