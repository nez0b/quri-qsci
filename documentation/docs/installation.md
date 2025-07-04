# Installation

This guide will help you install QSCI and all its dependencies.

## Prerequisites

- Python 3.9 or higher
- pip or uv package manager
- Git (for development installation)

## Quick Installation

### Basic Installation

Install the core QSCI package:

=== "Using pip"

    ```bash
    pip install quri-qsci
    ```

=== "Using uv (Recommended)"

    ```bash
    uv add quri-qsci
    ```

### With Optional Dependencies

Install with ffsim integration for advanced molecular ansatz support:

=== "Using pip"

    ```bash
    # With ffsim integration
    pip install quri-qsci[ffsim]
    
    # With all optional dependencies
    pip install quri-qsci[all]
    ```

=== "Using uv (Recommended)"

    ```bash
    # With ffsim integration
    uv add quri-qsci[ffsim]
    
    # With all optional dependencies  
    uv add quri-qsci[all]
    ```

## Development Installation

For development or to get the latest features:

### 1. Clone the Repository

```bash
git clone https://github.com/QunaSys/quri-qsci.git
cd quri-qsci
```

### 2. Create Virtual Environment

=== "Using uv"

    ```bash
    uv venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```

=== "Using venv"

    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    ```

### 3. Install Dependencies

=== "Using uv"

    ```bash
    uv sync
    ```

=== "Using pip"

    ```bash
    pip install -e .
    ```

## Dependencies

QSCI depends on the following packages:

### Core Dependencies
- **quri-parts**: Quantum circuit and operator manipulation
- **numpy**: Numerical computations
- **scipy**: Sparse matrix operations and linear algebra

### Optional Dependencies

#### ffsim Integration (`quri-qsci[ffsim]`)
- **ffsim**: Advanced molecular ansatz library for UCJ/LUCJ ansatz
- Provides the `ffsim_integration` submodule with:
  - UCJ (Unitary Coupled Cluster) ansatz support  
  - LUCJ (Low-rank UCJ) ansatz support
  - Molecular system creation utilities (H2, N2, etc.)
  - State conversion between ffsim and QURI Parts
  - High-level QSCI-ffsim integration interface

#### Other Optional Dependencies
- **quri-algo**: Time evolution and algorithm interfaces
- **qulacs**: High-performance quantum circuit simulator
- **pyscf**: Quantum chemistry calculations (for molecular examples)

### Development Dependencies
- **pytest**: Testing framework
- **mkdocs**: Documentation generation

## Verification

### Basic Installation

Verify your core installation:

```bash
# Test core functionality
python -c "import quri_qsci; print('✓ Core QSCI installed successfully')"
```

### ffsim Integration (if installed)

Verify ffsim integration:

```bash
# Test ffsim integration
python -c "import ffsim_integration; print('✓ ffsim integration available')"
```

### Running Tests

Test your installation with the test suite:

```bash
# Run all tests (automatically skips optional dependencies not installed)
pytest

# Run only core tests (skip ffsim tests)  
pytest -m "not ffsim"

# Run only ffsim tests (requires ffsim installation)
pytest -m "ffsim"
```

### Complete Example

Test with a simple QSCI calculation:

```python
from quri_qsci import VanillaQSCI
import numpy as np

# Create a simple 2x2 Hamiltonian
hamiltonian = np.array([[1.0, 0.5], [0.5, 2.0]])

# Run QSCI
qsci = VanillaQSCI(hamiltonian, num_states_pick_out=2)
result = qsci.run(num_eigenstates=1)

print("✓ QSCI installation successful!")
print(f"Ground state energy: {result.ground_state_energy:.6f}")
```

### ffsim Integration Example (requires ffsim)

Test molecular calculations with ffsim:

```python
from ffsim_integration import create_h2_molecule, run_lucj_qsci

# Create H2 molecule system
h2_system = create_h2_molecule(basis="sto-3g", bond_length=0.74)

# Run LUCJ-QSCI calculation
result = run_lucj_qsci(
    molecular_system=h2_system,
    n_reps=1,
    num_states_pick_out=10,
    max_iterations=5
)

print("✓ ffsim integration working!")
print(f"H2 LUCJ-QSCI energy: {result.qsci_energy:.6f} Ha")
```

## Optional: Quantum Chemistry Support

For molecular examples and PySCF integration:

```bash
pip install pyscf
```

This enables:
- Molecular Hamiltonian generation
- STO-3G and other basis sets  
- H2, H6 molecular examples

For Configuration Interaction calcualtions, we recommend [pyCI](https://pyci.qcdevs.org/intro.html)

```bash
git clone https://github.com/theochem/pyci.git && cd pyci
make
pip install pyCI
```


## Getting Help

If you encounter issues:

1. Check the [FAQ](faq.md)
2. Search existing [GitHub issues](https://github.com/QunaSys/quri-qsci/issues)
3. Create a new issue with:
   - Python version
   - Operating system
   - Complete error message
   - Minimal reproduction example

## Next Steps

Once installed, continue to:

- **[Quick Start](quickstart.md)**: Run your first QSCI calculation
- **[Examples](examples/basic_usage.md)**: Explore practical examples
- **[API Reference](api/)**: Dive into the technical details