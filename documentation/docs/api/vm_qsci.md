# VM-Aware QSCI API Reference

This page provides comprehensive API documentation for the VM-aware QSCI functionality, including quantum virtual machine integration, authentic time evolution, and resource estimation capabilities.

## Core VM Components

### VMFactory

::: qsci_vm_analysis.VMFactory

The `VMFactory` class provides factory methods for creating pre-configured QURI VM instances optimized for different quantum architectures.

**Key Methods:**
- `create_ideal_vm()`: Creates abstract ideal VM for LogicalCircuit analysis
- `create_star_vm(error_rate)`: Creates STAR architecture VM with specified error rate
- `create_star_vm_for_h6(error_rate)`: Creates STAR VM optimized for H6 calculations

### VMSampler

::: qsci_vm_analysis.VMSampler

The `VMSampler` class implements the quri-parts ConcurrentSampler interface, delegating circuit sampling to QURI VM instances with proper error handling.

**Key Features:**
- Direct delegation to `VM.sample()` 
- Comprehensive error wrapping
- Support for multiple circuit-shot pairs
- Clean integration with quri-parts sampling interface

### QURIVMInterface

::: qsci_vm_analysis.QURIVMInterface

The `QURIVMInterface` class provides a unified interface for QURI VM analysis and sampling across different quantum architectures.

**Key Capabilities:**
- Circuit analysis at LogicalCircuit and ArchLogicalCircuit levels
- Automatic VM instance management
- Sampler creation for different analysis levels
- Error propagation and handling

## VM-Aware Algorithms

### VMAwareSingleTimeTE_QSCI

::: qsci_vm_interface.VMAwareSingleTimeTE_QSCI

The `VMAwareSingleTimeTE_QSCI` class extends the standard TE-QSCI algorithm with comprehensive VM analysis capabilities, including authentic time evolution and resource estimation.

**Enhanced Features:**
- Authentic Trotter decomposition using quri-algo
- Fallback to simplified Trotter when quri-algo unavailable
- VM metrics calculation (overhead and fidelity impact)
- Integration with molecular Hamiltonian systems
- Multi-error-rate comparison support

### VMAwareTimeEvolutionCircuit

::: qsci_vm_interface.VMAwareTimeEvolutionCircuit

The `VMAwareTimeEvolutionCircuit` class wraps quantum circuits with VM analysis capabilities, providing caching and comprehensive circuit characterization.

**Key Features:**
- Analysis result caching for performance
- Summary generation with timing information
- Integration with VM interfaces
- Support for both LogicalCircuit and ArchLogicalCircuit analysis

## Result Classes

### TEQSCIVMResult

::: qsci_vm_interface.TEQSCIVMResult

The `TEQSCIVMResult` dataclass encapsulates enhanced TE-QSCI results with comprehensive VM analysis data.

**Included Data:**
- Standard QSCI results (energies, states)
- VM analysis results at multiple levels
- Performance metrics (overhead, fidelity impact)
- Circuit characterization data
- Algorithm metadata

### VMAnalysisResult

::: qsci_vm_interface.VMAnalysisResult

The `VMAnalysisResult` dataclass contains detailed VM circuit analysis information.

**Analysis Data:**
- LogicalCircuit and ArchLogicalCircuit analysis
- Circuit description and metadata
- Analysis timing information
- Resource estimation data

## Factory Functions

### create_vm_aware_te_qsci

::: qsci_vm_interface.create_vm_aware_te_qsci

Factory function for creating VM-aware TE-QSCI algorithms with VM sampling capabilities.

**Parameters:**
- `hamiltonian`: Target molecular Hamiltonian
- `evolution_time`: Time evolution parameter
- `vm_interface`: Pre-configured VM interface (optional)
- `sampling_level`: Analysis level (LogicalCircuit or ArchLogicalCircuit)
- `use_vm_sampling`: Enable VM-based sampling
- Additional TE-QSCI parameters

### create_star_vm_for_h6

::: qsci_vm_analysis.create_star_vm_for_h6

Factory function for creating STAR architecture VM interface optimized for H6 TE-QSCI calculations.

**Parameters:**
- `error_rate`: Physical error rate for STAR architecture (0.0 for ideal)

**Returns:**
- `QURIVMInterface`: Configured interface with ideal and STAR VMs

## Integration Points

### QURI VM Integration

The VM-aware QSCI functionality integrates seamlessly with QURI VM for realistic quantum hardware simulation:

- **Device Property Generation**: Uses `star_device.generate_device_property()` for STAR architecture configuration
- **VM Analysis**: Leverages `VM.analyze()` for circuit resource estimation
- **VM Sampling**: Utilizes `VM.sample()` for quantum sampling with noise models

### quri-algo Integration

Authentic time evolution leverages quri-algo for real quantum circuit generation:

- **TrotterTimeEvolutionCircuitFactory**: Generates authentic Trotter decomposition circuits
- **QubitHamiltonian**: Converts quri-parts Operators to quri-algo format
- **Fallback Mechanisms**: Graceful degradation when quri-algo unavailable

### Error Handling

Comprehensive error handling throughout the VM stack:

- **Import Error Handling**: Graceful fallback when dependencies unavailable
- **Runtime Error Wrapping**: Clear error messages with context
- **Validation**: Parameter validation and boundary checking
- **Resource Management**: Proper cleanup and error propagation

## Performance Considerations

### Analysis Caching

The VM-aware components implement intelligent caching for performance:

- **Circuit Analysis Caching**: Avoids redundant VM analysis calls
- **Level-Specific Caching**: Separate caches for LogicalCircuit and ArchLogicalCircuit
- **Cache Invalidation**: Proper cache management for accuracy

### Resource Optimization

Optimized for efficient quantum resource usage:

- **Minimal VM Overhead**: Efficient delegation patterns
- **Batch Processing**: Support for multiple circuit analysis
- **Memory Management**: Proper resource cleanup and management

### Scalability

Designed for molecular system scalability:

- **Hamiltonian Complexity**: Scales with molecular system size
- **Circuit Depth**: Handles realistic molecular time evolution
- **Error Rate Sensitivity**: Efficient multi-error-rate analysis