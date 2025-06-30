# Installation

This guide will help you install QSCI and all its dependencies.

## Prerequisites

- Python 3.9 or higher
- pip or uv package manager
- Git (for development installation)

## Quick Installation

### Using pip

```bash
pip install quri-qsci
```

### Using uv (Recommended)

```bash
uv add quri-qsci
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
- **quri-algo**: Time evolution and algorithm interfaces
- **qulacs**: High-performance quantum circuit simulator
- **pyscf**: Quantum chemistry calculations (for molecular examples)

### Development Dependencies
- **pytest**: Testing framework
- **mkdocs**: Documentation generation

## Verification

Verify your installation by running:

```python
import src.qsci_algorithms as qsci
from quri_parts.core.operator import Operator, pauli_label

# Create a simple Hamiltonian
hamiltonian = Operator()
hamiltonian += Operator({pauli_label("Z0"): -1.0})
hamiltonian += Operator({pauli_label("Z1"): -1.0})

print("QSCI installation successful!")
print(f"Hamiltonian: {hamiltonian}")
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