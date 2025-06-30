# Contributing to QSCI

Thank you for your interest in contributing to QSCI! This document provides guidelines for contributing to the project.

## Code of Conduct

We are committed to providing a welcoming and inspiring community for all. Please read and follow our Code of Conduct.

## Getting Started

### Development Setup

1. **Fork and Clone**
   ```bash
   git clone https://github.com/nez0b/quri-qsci.git
   cd quri-qsci
   ```

2. **Set Up Environment**
   ```bash
   uv venv
   source .venv/bin/activate
   uv sync
   ```

3. **Verify Installation**
   ```bash
   pytest tests/ -v
   ```

### Development Workflow

1. **Create a Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Changes**
   - Write code following our style guidelines
   - Add tests for new functionality
   - Update documentation as needed

3. **Test Your Changes**
   ```bash
   # Run the full test suite
   pytest
   
   # Run specific tests
   pytest tests/exact_diagonalizations/ -v
   
   # Check code coverage
   pytest --cov=src
   ```

4. **Commit and Push**
   ```bash
   git add .
   git commit -m "feat: add new QSCI variant"
   git push origin feature/your-feature-name
   ```

5. **Create Pull Request**
   - Use descriptive title and description
   - Reference any related issues
   - Ensure all checks pass

## Types of Contributions

### 1. Algorithm Development

#### Adding New QSCI Variants

When implementing a new QSCI variant:

```python
class NewQSCIVariant(QSCIBase):
    """New QSCI algorithm variant."""
    
    def __init__(self, hamiltonian, sampler, custom_parameter, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.custom_parameter = custom_parameter
    
    def run(self, input_states, total_shots, **kwargs):
        """Implement the new algorithm logic."""
        # Your implementation here
        pass
```

**Requirements:**
- Inherit from `QSCIBase`
- Implement required abstract methods
- Add comprehensive tests
- Update documentation

#### Extending Time Evolution Methods

For new time evolution approaches:

```python
def _create_custom_time_evolution(
    self,
    initial_state: CircuitQuantumState,
    evolution_time: float,
    custom_params: Dict[str, Any]
) -> CircuitQuantumState:
    """Implement custom time evolution method."""
    # Your implementation here
    pass
```

**Guidelines:**
- Follow the existing time evolution interface
- Support both exact and approximate methods
- Include proper error handling
- Validate against known solutions

### 2. Testing Contributions

#### Adding New Tests

Follow the uniform superposition testing methodology:

```python
@pytest.mark.exact_diagonalization
def test_new_variant_uniform_superposition_precision(self, exact_ground_state_test_systems):
    """Test NewVariant achieves machine precision with uniform superposition."""
    
    for system_name, system_data in exact_ground_state_test_systems.items():
        # Create uniform superposition state
        uniform_state = create_uniform_superposition(system_data['n_qubits'])
        
        # Run with complete subspace coverage
        algorithm = NewQSCIVariant(
            hamiltonian=system_data['quri_hamiltonian'],
            sampler=create_qulacs_vector_concurrent_sampler(),
            num_states_pick_out=2**system_data['n_qubits']
        )
        
        result = algorithm.run([uniform_state], total_shots=2500)
        
        # Validate machine precision
        energy_error = abs(result.ground_state_energy - system_data['exact_ground_energy'])
        assert energy_error < 1e-8
```

**Testing Checklist:**
- [ ] Use uniform superposition initial states
- [ ] Apply complete subspace coverage (`2**n_qubits`)
- [ ] Set appropriate shot counts (2000-3000)
- [ ] Use machine precision tolerance (1e-8)
- [ ] Add proper test markers
- [ ] Include comprehensive docstrings

#### Test Utilities

When creating new test utilities:

```python
def create_new_model_hamiltonian(n_qubits: int, parameter: float) -> np.ndarray:
    """Create Hamiltonian for new quantum model.
    
    Args:
        n_qubits: Number of qubits
        parameter: Model parameter
        
    Returns:
        Sparse Hamiltonian matrix
    """
    # Implementation using scipy.sparse
    hamiltonian = create_sparse_matrix(n_qubits, parameter)
    
    # Validate properties
    assert is_hermitian(hamiltonian), "Hamiltonian must be Hermitian"
    
    return hamiltonian
```

### 3. Documentation

#### API Documentation

Use Google-style docstrings:

```python
def example_function(param1: int, param2: str = "default") -> bool:
    """Brief description of the function.
    
    Longer description with more details about the function's behavior,
    use cases, and important notes.
    
    Args:
        param1: Description of param1
        param2: Description of param2 with default value
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: When param1 is negative
        TypeError: When param2 is not a string
        
    Examples:
        Basic usage:
        
        ```python
        result = example_function(42, "test")
        print(result)  # True
        ```
        
        Advanced usage with error handling:
        
        ```python
        try:
            result = example_function(-1)
        except ValueError as e:
            print(f"Error: {e}")
        ```
    """
    if param1 < 0:
        raise ValueError("param1 must be non-negative")
    
    return isinstance(param2, str) and param1 > 0
```

#### User Guides

When adding user guides:

1. **Start with motivation** - Why is this feature useful?
2. **Provide complete examples** - Working code that users can copy
3. **Explain key concepts** - Background theory when relevant
4. **Include best practices** - How to use the feature effectively
5. **Add troubleshooting** - Common issues and solutions

### 4. Performance Optimization

#### Profiling

Before optimizing, profile your code:

```python
import cProfile
import pstats

def profile_algorithm():
    """Profile QSCI algorithm performance."""
    profiler = cProfile.Profile()
    profiler.enable()
    
    # Run your algorithm
    result = algorithm.run(initial_state, total_shots=1000)
    
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.sort_stats('cumulative').print_stats(10)
```

#### Optimization Guidelines

1. **Measure first** - Always profile before optimizing
2. **Focus on bottlenecks** - Optimize the slowest parts first
3. **Maintain correctness** - Verify results don't change
4. **Document changes** - Explain performance improvements
5. **Add benchmarks** - Include performance tests

## Code Style Guidelines

### Python Style

Follow PEP 8 with these specifics:

- **Line length**: 88 characters (Black default)
- **Imports**: Use absolute imports, group by standard/third-party/local
- **Type hints**: Required for all public functions
- **Docstrings**: Google style for all public classes and functions

### Code Formatting

Use Black for consistent formatting:

```bash
black src/ tests/
```

### Import Organization

```python
# Standard library
import os
import sys
from typing import Dict, List, Optional

# Third party
import numpy as np
import scipy.sparse

# QURI ecosystem
from quri_parts.core.operator import Operator
from quri_parts.core.state import GeneralCircuitQuantumState

# Local imports
from src.qsci_algorithms import VanillaQSCI
```

## Testing Standards

### Test Requirements

- **100% coverage** for new algorithm implementations
- **Machine precision validation** for exact diagonalization tests
- **Cross-platform compatibility** (Windows, macOS, Linux)
- **Performance regression tests** for optimization changes

### Test Categories

Mark your tests appropriately:

```python
@pytest.mark.exact_ground_state_precision   # Machine precision tests
@pytest.mark.te_qsci_single_time            # Algorithm-specific tests
@pytest.mark.verification                   # Mathematical validation
@pytest.mark.integration                    # Framework integration
@pytest.mark.molecular                      # Molecular systems (slow)
@pytest.mark.regression                     # Bug regression tests
```

### Test Data

Use the standardized test utilities:

```python
# Use hamiltonian factory
hamiltonian = create_tfim_hamiltonian(n_qubits=2, h=1.0)

# Use uniform superposition
initial_state = create_uniform_superposition(n_qubits=2)

# Use complete subspace coverage
num_states_pick_out = 2**n_qubits
```

## Pull Request Guidelines

### PR Checklist

- [ ] **Tests pass**: All existing tests continue to pass
- [ ] **New tests added**: For new functionality
- [ ] **Documentation updated**: API docs and user guides
- [ ] **Performance verified**: No unexpected regressions
- [ ] **Type hints included**: For all new public functions
- [ ] **Changelog updated**: For user-facing changes

### PR Description Template

```markdown
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix (non-breaking change fixing an issue)
- [ ] New feature (non-breaking change adding functionality)
- [ ] Breaking change (fix/feature causing existing functionality to change)
- [ ] Documentation update

## Testing
- [ ] All existing tests pass
- [ ] New tests added for new functionality
- [ ] Manual testing performed

## Performance Impact
- [ ] No performance impact
- [ ] Performance improvement (include benchmarks)
- [ ] Performance regression (justified and documented)

## Related Issues
Fixes #123
Related to #456
```

### Review Process

1. **Automated checks** must pass (tests, linting, type checking)
2. **Code review** by at least one maintainer
3. **Documentation review** for user-facing changes
4. **Performance review** for optimization changes
5. **Final approval** and merge

## Release Process

### Version Numbering

We follow [Semantic Versioning](https://semver.org/):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist

- [ ] **All tests pass** on all supported platforms
- [ ] **Documentation updated** including API changes
- [ ] **Changelog updated** with user-facing changes
- [ ] **Version bumped** in appropriate files
- [ ] **Release notes** prepared
- [ ] **Performance benchmarks** run and documented

## Getting Help

### Development Questions

- **GitHub Discussions**: For design questions and ideas
- **GitHub Issues**: For bugs and feature requests
- **Documentation**: Check existing docs first

### Maintainer Contact

- **Response time**: Usually within 2-3 business days
- **Urgency**: Use appropriate labels (bug, enhancement, question)
- **Context**: Include minimal reproduction examples

## Recognition

Contributors are recognized in:

- **CHANGELOG.md**: For significant contributions
- **GitHub contributors**: Automatic recognition
- **Documentation**: For major feature additions
- **Release notes**: For important improvements

Thank you for contributing to QSCI! Your contributions help advance quantum computing research and make quantum algorithms more accessible to the scientific community.