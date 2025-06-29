# Architecture and System Overview

## System Architecture

The TE-QSCI system is designed as a modular, layered architecture that integrates seamlessly with the QURI ecosystem. Each layer has well-defined responsibilities and interfaces, ensuring maintainability and extensibility.

```
┌─────────────────────────────────────────────────────────────┐
│                        TE-QSCI System                       │
├─────────────────────────────────────────────────────────────┤
│  User Interface & Examples                                  │
│  ├── examples/h6_molecule_example.py (H6 study)            │
│  ├── examples/quick_test.py (basic demo)                   │
│  └── tests/test_implementation.py (validation)             │
├─────────────────────────────────────────────────────────────┤
│  Algorithm Interface Layer (quri-algo compatible)          │
│  ├── QSCIAlgorithmBase (base interface)                    │
│  ├── VanillaQSCIAlgorithm (standard QSCI)                  │
│  ├── SingleTimeTeQSCIAlgorithm (TE-QSCI variant)           │
│  ├── TimeAverageTeQSCIAlgorithm (time-average variant)     │
│  └── StateVectorTeQSCIAlgorithm (exact variant)            │
├─────────────────────────────────────────────────────────────┤
│  Core Algorithm Implementation                              │
│  ├── QSCIBase (abstract base)                              │
│  ├── VanillaQSCI (standard implementation)                 │
│  └── TimeEvolvedQSCI (TE-QSCI implementation)              │
├─────────────────────────────────────────────────────────────┤
│  QURI VM Integration Layer                                  │
│  ├── VMCircuitAnalyzer (circuit analysis)                  │
│  ├── CircuitResourceEstimate (resource estimation)         │
│  └── QSCIAnalysis (algorithm-specific analysis)            │
├─────────────────────────────────────────────────────────────┤
│  QURI Parts Foundation                                      │
│  ├── Operators (Hamiltonians, Pauli strings)               │
│  ├── Quantum States (circuit-based states)                 │
│  ├── Quantum Circuits (gate sequences)                     │
│  ├── Samplers (measurement simulation)                     │
│  └── Backends (qulacs, etc.)                               │
└─────────────────────────────────────────────────────────────┘
```

## Design Principles

The architecture follows five core principles:

### 1. Modularity
Each layer has well-defined responsibilities and clean interfaces. The separation of concerns allows for independent development and testing of each component.

### 2. Extensibility
The system is designed to easily accommodate new TE-QSCI variants and algorithm improvements. The abstract base classes provide clear extension points.

### 3. Compatibility
Full integration with the QURI ecosystem ensures that existing tools and workflows can leverage TE-QSCI algorithms without modification.

### 4. Performance
Efficient implementation strategies, including sparse matrix representations and optimized sampling, enable scaling to larger molecular systems.

### 5. Testability
Comprehensive testing at each layer, with particular emphasis on exact diagonalization validation, ensures correctness and reliability.

## Project Organization

```
te-qsci/
├── src/                           # Core implementation
│   ├── __init__.py               # Package exports
│   ├── qsci_algorithms.py        # Core QSCI algorithms
│   ├── qsci_algo_interface.py    # QURI-algo interfaces
│   └── qsci_vm_analysis.py       # QURI VM integration
├── examples/                      # Usage examples
│   ├── h6_molecule_example.py    # Complete H6 study
│   ├── simplified_h6_example.py  # Simplified demo
│   ├── quick_test.py            # Basic functionality test
│   └── main.py                  # Entry point
├── tests/                        # Test suite
│   ├── exact_diagonalizations/   # Exact validation tests
│   ├── utils/                    # Test utilities
│   └── test_implementation.py    # Implementation tests
├── requirements.txt              # Dependencies
├── pyproject.toml               # Project configuration
├── README.md                    # User documentation
└── DESIGN_DOCUMENT.md           # Technical documentation
```

## Module Dependencies

The dependency structure is designed to minimize coupling while maximizing functionality:

```
┌─────────────────────────────────────────┐
│                Examples                 │
│  h6_molecule_example.py                 │
│  simplified_h6_example.py               │
│  quick_test.py                         │
└─────────────────┬───────────────────────┘
                  │ imports
┌─────────────────▼───────────────────────┐
│           src/__init__.py               │
│        (Package Interface)              │
└─────────────────┬───────────────────────┘
                  │ exports
┌─────────────────▼───────────────────────┐
│        qsci_algo_interface.py           │
│     (QURI-Algo Compatibility)           │
└─────────────────┬───────────────────────┘
                  │ imports
┌─────────────────▼───────────────────────┐
│         qsci_algorithms.py              │
│      (Core Implementation)              │
└─────────────────┬───────────────────────┘
                  │ imported by
┌─────────────────▼───────────────────────┐
│        qsci_vm_analysis.py              │
│       (QURI VM Integration)             │
└─────────────────────────────────────────┘
```

## Class Hierarchy

The class hierarchy provides clear inheritance relationships and extension points:

```
Algorithm (quri-algo)
├── QSCIAlgorithmBase
│   ├── VanillaQSCIAlgorithm
│   ├── SingleTimeTeQSCIAlgorithm
│   ├── TimeAverageTeQSCIAlgorithm
│   └── StateVectorTeQSCIAlgorithm
│
QSCIBase (ABC)
├── VanillaQSCI
└── TimeEvolvedQSCI

Analysis (quri-algo)
└── QSCIAnalysis
    └── VMEnabledQSCIAnalysis
```

## Data Flow Architecture

The data flow follows a clear path from user input to final results:

```
User Input (Hamiltonian, Parameters)
        │
        ▼
Factory Function (create_qsci_algorithm)
        │
        ▼
Algorithm Interface (QSCIAlgorithmBase)
        │
        ▼
Core Implementation (QSCIBase subclass)
        │
        ▼
QURI Parts Operations (States, Operators, Sampling)
        │
        ▼
QURI VM Analysis (Resource Estimation)
        │
        ▼
Results (QSCIResult, Analysis Data)
```

## Key Components

### Core Algorithm Layer
- **QSCIBase**: Abstract foundation defining the common interface
- **VanillaQSCI**: Standard QSCI implementation using computational basis sampling
- **TimeEvolvedQSCI**: Base class for all TE-QSCI variants with time evolution capabilities

### Algorithm Interface Layer
- **QSCIAlgorithmBase**: quri-algo compatible interface
- **Algorithm Variants**: Specific implementations for each TE-QSCI variant
- **Factory Functions**: Easy algorithm creation and configuration

### Integration Layer
- **QURI Parts Integration**: Native operator and state management
- **QURI VM Analysis**: Circuit resource estimation and optimization
- **Backend Compatibility**: Support for various quantum simulators

## Performance Considerations

The architecture is optimized for performance through:

- **Sparse Matrix Operations**: Efficient handling of large Hamiltonian matrices
- **Concurrent Sampling**: Parallel measurement processing
- **Memory Optimization**: Smart memory usage for large systems
- **Adaptive Algorithms**: Dynamic parameter adjustment for optimal performance

## Navigation

- **Next**: [QURI Integration](quri_integration.md) - Learn about ecosystem integration
- **Related**: [Time Evolution](time_evolution.md) - Understand time evolution implementation
- **See Also**: [Performance](performance.md) - Detailed performance analysis