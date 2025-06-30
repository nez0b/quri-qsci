# Changelog

All notable changes to the QSCI project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive MkDocs documentation site
- API reference with auto-generated documentation
- Interactive examples and tutorials

### Changed
- Improved error handling in time evolution methods
- Enhanced performance for large molecular systems

### Fixed
- Memory optimization for sparse matrix operations

## [0.0.1] - 2025-06-28

### Added
- **Complete QSCI Algorithm Suite**: VanillaQSCI and all TE-QSCI variants
- **Uniform Superposition Breakthrough**: Machine precision testing methodology
- **Comprehensive Testing Framework**: 41/41 exact diagonalization tests passing
- **QURI Ecosystem Integration**: Full compatibility with QURI Parts, Algo, and VM
- **Time Evolution Support**: Trotter decomposition using TrotterTimeEvolutionCircuitFactory
- **Scientific Validation**: H6 molecule study with 919-term Hamiltonian
- **Architecture Analysis**: Circuit resource estimation for quantum hardware

#### Algorithm Implementations
- **VanillaQSCI**: Standard QSCI algorithm with computational basis sampling
- **SingleTimeTE_QSCI**: Time-evolved QSCI at single evolution time
- **TimeAverageTE_QSCI**: Time-evolved QSCI averaged over multiple times  
- **StateVectorTE_QSCI**: Exact state vector TE-QSCI for validation

#### Testing
- Machine precision accuracy (< 1e-15) with uniform superposition states
- Complete subspace coverage for deterministic validation
- Trotter and exact unitary time evolution methods
- QURI VM integration for circuit analysis
- Comprehensive Hamiltonian support (TFIM, Heisenberg, molecular)

#### Testing Achievements
- **100% Test Success Rate**: All 41 exact diagonalization tests passing
- **Machine Precision Validation**: Uniform superposition + complete subspace coverage
- **Cross-Algorithm Validation**: All QSCI variants achieve consistent results
- **Molecular System Validation**: H2, H6 molecular studies
- **Model System Coverage**: TFIM, Heisenberg, random sparse Hamiltonians

### Technical Implementation
- **QSCIBase**: Abstract foundation with common functionality
- **TimeEvolvedQSCI**: TE-QSCI base implementation with three execution modes
- **QURI Parts Integration**: Native operator and state management
- **Algorithm Interfaces**: quri-algo compatible interfaces
- **Resource Estimation**: QURI VM integration for circuit analysis

### Performance Optimizations
- Sparse matrix operations for large Hamiltonians
- Concurrent sampling for improved measurement efficiency
- Adaptive eigenvalue computation based on subspace size
- Memory-efficient state vector calculations




## Future Roadmap

### Version 1.1.0 (Planned)
- [ ] Enhanced error mitigation integration
- [ ] Higher-order Trotter decomposition methods
- [ ] Adaptive configuration selection algorithms
- [ ] Hardware-specific circuit optimizations



##

## Contributors

- 

## Acknowledgments


---

For detailed technical changes, see the [GitHub releases](https://github.com/nez0b/quri-qsci/releases) page.