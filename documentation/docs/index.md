# QSCI - Quantum Selected Configuration Interaction

Welcome to the comprehensive documentation for **QSCI** (Quantum Selected Configuration Interaction) and **TE-QSCI** (Time-Evolved QSCI) algorithms implemented using the QURI ecosystem.

## What is QSCI?

QSCI is a quantum algorithm for quantum chemistry that leverages time evolution to systematically explore electronic configurations without requiring variational optimization. This implementation provides multiple algorithm variants optimized for different use cases and quantum hardware architectures.

## Key Features

- **Multiple Algorithm Variants**: VanillaQSCI, SingleTimeTE_QSCI, TimeAverageTE_QSCI, and StateVectorTE_QSCI
- **QURI Ecosystem Integration**: Full compatibility with QURI Parts, QURI Algo, and QURI VM
- **Time Evolution Methods**: Support for both Trotter decomposition and exact unitary evolution
- **Quantum Hardware Analysis**: Powered by [QURI-VM](https://github.com/QunaSys/quri-vm) , it supports circuit resource estimation for different quantum architectures, including ideal quantum circuit and FTQC devices


## Algorithm Variants

### VanillaQSCI
Standard QSCI algorithm using computational basis state sampling from quantum circuits.

### SingleTimeTE_QSCI  
Time-evolved QSCI at a single evolution time point for systematic configuration exploration.

### TimeAverageTE_QSCI
Time-evolved QSCI averaged over multiple evolution times to improve sampling diversity.

### StateVectorTE_QSCI
Time-evolved QSCI with direct state vector calculation for exact simulation studies.


## Scientific Background

TE-QSCI generates important electronic configurations through time evolution:

$$|\psi(t)\rangle = e^{-i\hat{H}t}|\psi_I\rangle = |\psi_I\rangle - i\hat{H}t|\psi_I\rangle + \frac{(-i\hat{H}t)^2}{2!}|\psi_I\rangle + \ldots$$

Where the $k$-th order term includes up to $2k$-th order excitations, systematically exploring the configuration space needed for accurate quantum chemistry calculations.

## Getting Started

1. **[Installation Guide](installation.md)** - Set up QSCI with all dependencies
2. **[Quick Start](quickstart.md)** - Run your first QSCI calculation  
3. **[Algorithm Overview](algorithms/overview.md)** - Understand the different QSCI variants
4. **[Examples](examples/basic_usage.md)** - Learn from practical examples

## API Reference

Explore the complete API documentation:

- **[Core Algorithms](api/qsci_algorithms.md)** - VanillaQSCI, TimeEvolvedQSCI, and variants
- **[Algorithm Interfaces](api/qsci_algo_interface.md)** - QURI-algo compatible interfaces  
- **[VM Analysis](api/qsci_vm_analysis.md)** - Circuit analysis and resource estimation


## Contributing

We welcome contributions! Please see our [Contributing Guide](contributing.md) for details on:

- Adding new QSCI variants
- Extending the testing framework
- Improving documentation
- Reporting issues

## Citation

If you use QSCI in your research, please cite:

```bibtex
@article{Kanno:2023rfr,
    author = "Kanno, Keita and Kohda, Masaya and Imai, Ryosuke and Koh, Sho and Mitarai, Kosuke and Mizukami, Wataru and Nakagawa, Yuya O.",
    title = "{Quantum-Selected Configuration Interaction: classical diagonalization of Hamiltonians in subspaces selected by quantum computers}",
    eprint = "2302.11320",
    archivePrefix = "arXiv",
    primaryClass = "quant-ph",
    month = "2",
    year = "2023"
}

@article{Mikkelsen:2024onz,
    author = "Mikkelsen, Mathias and Nakagawa, Yuya O.",
    title = "{Quantum-selected configuration interaction with time-evolved state}",
    eprint = "2412.13839",
    archivePrefix = "arXiv",
    primaryClass = "quant-ph",
    month = "12",
    year = "2024"
}
```

## License

This project is licensed under the MIT License - see the [License](license.md) file for details.