# TE-QSCI Project Structure

This project implements the Time-Evolved Quantum Selected Configuration Interaction (TE-QSCI) algorithm using the QURI Parts framework.

## Directory Structure

```
te-qsci/
├── src/                    # Core implementation modules
│   ├── __init__.py        # Package initialization with exports
│   ├── qsci_algorithms.py # Core QSCI and TE-QSCI algorithms
│   ├── qsci_algo_interface.py # QURI-algo compatible interfaces
│   └── qsci_vm_analysis.py # QURI VM integration for circuit analysis
├── examples/              # Example scripts and demonstrations
│   ├── h6_molecule_example.py # Complete H6 molecule study
│   ├── simplified_h6_example.py # Simplified demonstration
│   ├── quick_test.py      # Quick functionality test
│   └── main.py           # Entry point example
├── tests/                 # Test files
│   └── test_implementation.py # Implementation tests
├── requirements.txt       # Python dependencies
├── pyproject.toml        # Project configuration
└── README.md             # Main project documentation
```

## Core Modules

### `src/qsci_algorithms.py`
- **VanillaQSCI**: Standard quantum selected configuration interaction
- **TimeEvolvedQSCI**: TE-QSCI implementation with variants:
  - Single-time TE-QSCI
  - Time-average TE-QSCI  
  - State vector TE-QSCI

### `src/qsci_algo_interface.py`
- QURI-algo compatible algorithm interfaces
- Algorithm, AlgorithmResult, and Analysis classes
- Factory functions for easy algorithm creation

### `src/qsci_vm_analysis.py`
- QURI VM integration for circuit analysis
- LogicalCircuit and ArchLogicalCircuit level analysis
- Resource estimation and architecture optimization

## Examples

### Quick Test
```bash
source .venv/bin/activate
python examples/quick_test.py
```

### H6 Molecule Study
```bash
source .venv/bin/activate
python examples/h6_molecule_example.py
```

### Simplified Demonstration
```bash
source .venv/bin/activate
python examples/simplified_h6_example.py
```

## Key Features

✅ **Optimization-free state preparation** - No VQE or circuit optimization needed  
✅ **Multiple algorithm variants** - Single-time, time-average, state vector  
✅ **Architecture awareness** - STAR topology optimization  
✅ **Extensible design** - Easy to add new variants and architectures  
✅ **Full QURI SDK integration** - Compatible with quri-parts, quri-algo, quri-vm  

## Installation

1. Create virtual environment:
```bash
source .venv/bin/activate
```

2. Install dependencies:
```bash
uv add numpy scipy matplotlib 'quri-parts[qulacs]'
```

3. Run tests:
```bash
python examples/quick_test.py
```

## Algorithm Theory

The TE-QSCI algorithm leverages time evolution to generate important electronic configurations:

```
|ψ(t)⟩ = e^(-iĤt)|ψ_I⟩ = |ψ_I⟩ - iĤt|ψ_I⟩ + (−iĤt)²/2!|ψ_I⟩ + ...
```

Where the k-th order term includes up to 2k-th order excitations, systematically exploring the configuration space needed for accurate quantum chemistry calculations.