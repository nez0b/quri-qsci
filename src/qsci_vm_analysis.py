"""
QURI VM integration for QSCI algorithms.

This module provides circuit analysis capabilities at both LogicalCircuit and 
ArchLogicalCircuit levels using QURI VM, enabling detailed resource estimation
for QSCI algorithms on different quantum architectures.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union
import numpy as np

from quri_parts.circuit import NonParametricQuantumCircuit, QuantumCircuit
from quri_parts.core.state import CircuitQuantumState
from quri_parts.backend.units import TimeValue, TimeUnit

from .qsci_algo_interface import (
    LoweringLevel, Analysis, QSCIAnalysis, QSCIAlgorithmBase
)

# QURI VM related imports
from quri_vm import VM
from quri_parts.backend.devices import star_device
from quri_parts.backend.units import TimeUnit, TimeValue

QURI_VM_AVAILABLE = True
print("QURI VM successfully imported and available")


@dataclass
class ArchitectureInfo:
    """Information about quantum architecture."""
    name: str
    num_physical_qubits: int
    connectivity: str  # "star", "linear", "grid", etc.
    gate_fidelities: Dict[str, float]
    decoherence_times: Dict[str, float]  # T1, T2 times in microseconds
    gate_times: Dict[str, float]  # Gate execution times in microseconds


@dataclass
class CircuitResourceEstimate:
    """Resource estimates for quantum circuit execution."""
    gate_count: Dict[str, int]
    circuit_depth: int
    physical_qubit_count: int
    logical_qubit_count: int
    execution_time: TimeValue
    fidelity_estimate: float
    swap_overhead: int = 0


class VMCircuitAnalyzer:
    """Circuit analyzer using QURI VM capabilities."""
    
    def __init__(self, architecture: Optional[ArchitectureInfo] = None):
        """Initialize VM circuit analyzer.
        
        Args:
            architecture: Quantum architecture information
        """
        self.architecture = architecture or self._get_default_star_architecture()
    
    def analyze_logical_circuit(
        self, 
        circuit: NonParametricQuantumCircuit
    ) -> CircuitResourceEstimate:
        """Analyze circuit at LogicalCircuit level.
        
        Args:
            circuit: Quantum circuit to analyze
            
        Returns:
            CircuitResourceEstimate with logical-level analysis
        """
        # Logical circuit analysis - no architecture-specific optimizations
        gate_count = self._count_gates(circuit)
        circuit_depth = circuit.depth
        logical_qubit_count = circuit.qubit_count
        
        # Estimate execution time based on gate counts and types
        execution_time_us = sum(
            count * self._get_gate_time(gate_type)
            for gate_type, count in gate_count.items()
        )
        
        # Estimate fidelity based on gate fidelities
        fidelity_estimate = self._estimate_circuit_fidelity(gate_count)
        
        return CircuitResourceEstimate(
            gate_count=gate_count,
            circuit_depth=circuit_depth,
            physical_qubit_count=logical_qubit_count,  # Same at logical level
            logical_qubit_count=logical_qubit_count,
            execution_time=TimeValue(execution_time_us, TimeUnit.MICROSECOND),
            fidelity_estimate=fidelity_estimate
        )
    
    def analyze_arch_logical_circuit(
        self, 
        circuit: NonParametricQuantumCircuit
    ) -> CircuitResourceEstimate:
        """Analyze circuit at ArchLogicalCircuit level with architecture mapping.
        
        Args:
            circuit: Quantum circuit to analyze
            
        Returns:
            CircuitResourceEstimate with architecture-aware analysis
        """
        # Start with logical analysis
        logical_estimate = self.analyze_logical_circuit(circuit)
        
        # Apply architecture-specific optimizations and overhead
        if self.architecture.connectivity == "star":
            arch_estimate = self._analyze_star_architecture(circuit, logical_estimate)
        else:
            # For other architectures, use logical estimate with connectivity overhead
            arch_estimate = self._apply_connectivity_overhead(logical_estimate)
        
        return arch_estimate
    
    def _get_default_star_architecture(self) -> ArchitectureInfo:
        """Get default STAR architecture configuration for H6 calculations."""
        return ArchitectureInfo(
            name="STAR",
            num_physical_qubits=25,  # H6 needs 12 logical qubits + overhead
            connectivity="star",
            gate_fidelities={
                "single": 0.9999,     # High fidelity for STAR logical gates
                "two_qubit": 0.999,   # CNOT between logical qubits
                "measurement": 0.99   # Logical measurement fidelity
            },
            decoherence_times={
                "T1": 1000.0,  # microseconds (STAR has long coherence)
                "T2": 500.0    # microseconds (dephasing time)
            },
            gate_times={
                "single": 1.0,      # microseconds (includes QEC overhead)
                "two_qubit": 5.0,   # microseconds (CNOT with error correction)
                "measurement": 10.0 # microseconds (includes syndrome extraction)
            }
        )
    
    def _count_gates(self, circuit: NonParametricQuantumCircuit) -> Dict[str, int]:
        """Count gates by type in the circuit."""
        gate_count = {"single": 0, "two_qubit": 0, "measurement": 0}
        
        # This is a simplified implementation
        # In practice, would iterate through circuit gates
        gate_count["single"] = len(circuit.gates) // 2  # Rough estimate
        gate_count["two_qubit"] = len(circuit.gates) // 4  # Rough estimate
        
        return gate_count
    
    def _get_gate_time(self, gate_type: str) -> float:
        """Get execution time for gate type in microseconds."""
        return self.architecture.gate_times.get(gate_type, 0.1)
    
    def _estimate_circuit_fidelity(self, gate_count: Dict[str, int]) -> float:
        """Estimate overall circuit fidelity based on gate counts."""
        fidelity = 1.0
        
        for gate_type, count in gate_count.items():
            gate_fidelity = self.architecture.gate_fidelities.get(gate_type, 0.99)
            fidelity *= gate_fidelity ** count
        
        return fidelity
    
    def _analyze_star_architecture(
        self, 
        circuit: NonParametricQuantumCircuit,
        logical_estimate: CircuitResourceEstimate
    ) -> CircuitResourceEstimate:
        """Analyze circuit for STAR architecture connectivity."""
        # STAR architecture has one central qubit connected to all others
        # This reduces SWAP overhead for certain circuit patterns
        
        # Estimate SWAP overhead based on circuit structure
        swap_overhead = self._estimate_star_swap_overhead(circuit)
        
        # Update gate counts with SWAP gates
        arch_gate_count = logical_estimate.gate_count.copy()
        arch_gate_count["two_qubit"] += swap_overhead
        
        # Update execution time
        swap_time_us = swap_overhead * self._get_gate_time("two_qubit")
        arch_execution_time = TimeValue(
            logical_estimate.execution_time.in_ns() / 1000 + swap_time_us,
            TimeUnit.MICROSECOND
        )
        
        # Update fidelity with SWAP overhead
        swap_fidelity = self.architecture.gate_fidelities["two_qubit"] ** swap_overhead
        arch_fidelity = logical_estimate.fidelity_estimate * swap_fidelity
        
        return CircuitResourceEstimate(
            gate_count=arch_gate_count,
            circuit_depth=logical_estimate.circuit_depth + swap_overhead // 2,
            physical_qubit_count=logical_estimate.logical_qubit_count,
            logical_qubit_count=logical_estimate.logical_qubit_count,
            execution_time=arch_execution_time,
            fidelity_estimate=arch_fidelity,
            swap_overhead=swap_overhead
        )
    
    def _estimate_star_swap_overhead(self, circuit: NonParametricQuantumCircuit) -> int:
        """Estimate SWAP overhead for STAR architecture."""
        # Simplified heuristic: assume some fraction of gates need SWAPs
        total_gates = len(circuit.gates)
        swap_fraction = 0.1  # 10% of gates require SWAPs in STAR topology
        return int(total_gates * swap_fraction)
    
    def _apply_connectivity_overhead(
        self, 
        logical_estimate: CircuitResourceEstimate
    ) -> CircuitResourceEstimate:
        """Apply generic connectivity overhead."""
        # Add 20% overhead for non-star architectures
        overhead_factor = 1.2
        
        arch_gate_count = logical_estimate.gate_count.copy()
        arch_gate_count["two_qubit"] = int(arch_gate_count["two_qubit"] * overhead_factor)
        
        arch_execution_time = TimeValue(
            (logical_estimate.execution_time.in_ns() / 1000) * overhead_factor,
            TimeUnit.MICROSECOND
        )
        
        return CircuitResourceEstimate(
            gate_count=arch_gate_count,
            circuit_depth=int(logical_estimate.circuit_depth * overhead_factor),
            physical_qubit_count=logical_estimate.logical_qubit_count,
            logical_qubit_count=logical_estimate.logical_qubit_count,
            execution_time=arch_execution_time,
            fidelity_estimate=logical_estimate.fidelity_estimate * 0.95  # 5% fidelity loss
        )


class QSCIVMAnalysis(QSCIAnalysis):
    """Enhanced QSCI analysis using QURI VM capabilities."""
    
    def __init__(
        self,
        lowering_level: LoweringLevel,
        circuit_estimates: Dict[str, CircuitResourceEstimate],
        total_shots: int,
        architecture: Optional[ArchitectureInfo] = None
    ):
        """Initialize VM-enhanced QSCI analysis.
        
        Args:
            lowering_level: Analysis level
            circuit_estimates: Resource estimates for circuits
            total_shots: Total measurement shots
            architecture: Architecture information
        """
        self.circuit_estimates = circuit_estimates
        self.architecture = architecture
        self.total_shots = total_shots
        
        # Convert circuit estimates to base class format
        circuit_gate_count = {
            name: est.gate_count for name, est in circuit_estimates.items()
        }
        circuit_depth = {
            name: est.circuit_depth for name, est in circuit_estimates.items()
        }
        circuit_latency = {
            name: est.execution_time for name, est in circuit_estimates.items()
        }
        circuit_execution_count = {"total": total_shots}
        circuit_fidelities = {
            name: est.fidelity_estimate for name, est in circuit_estimates.items()
        }
        circuit_qubit_count = {
            name: est.physical_qubit_count for name, est in circuit_estimates.items()
        }
        
        # Create dummy QSCI result for base class
        from .qsci_algorithms import QSCIResult, QSCIVariant
        dummy_qsci_result = QSCIResult(
            eigenvalues=[0.0],
            eigenstates=[],
            selected_states=[],
            subspace_dimension=100,
            total_shots=total_shots,
            algorithm_variant=QSCIVariant.VANILLA
        )
        
        super().__init__(lowering_level, dummy_qsci_result)
        
        # Override with VM analysis
        self.circuit_gate_count = circuit_gate_count
        self.circuit_depth = circuit_depth
        self.circuit_latency = circuit_latency
        self.circuit_execution_count = circuit_execution_count
        self.circuit_fidelities = circuit_fidelities
        self.circuit_qubit_count = circuit_qubit_count
    
    @property
    def total_latency(self) -> TimeValue:
        """Total latency including all circuit executions."""
        total_us = 0.0
        
        for name, estimate in self.circuit_estimates.items():
            # Add execution time multiplied by number of shots
            circuit_shots = self.total_shots // len(self.circuit_estimates)
            total_us += (estimate.execution_time.in_ns() / 1000) * circuit_shots
        
        return TimeValue(total_us, TimeUnit.MICROSECOND)
    
    @property
    def max_physical_qubit_count(self) -> int:
        """Maximum physical qubit count across all circuits."""
        if self.circuit_estimates:
            return max(est.physical_qubit_count for est in self.circuit_estimates.values())
        return 0
    
    @property
    def total_swap_overhead(self) -> int:
        """Total SWAP gate overhead across all circuits."""
        return sum(est.swap_overhead for est in self.circuit_estimates.values())
    
    @property
    def average_fidelity(self) -> float:
        """Average fidelity across all circuits."""
        if self.circuit_estimates:
            return np.mean([est.fidelity_estimate for est in self.circuit_estimates.values()])
        return 1.0


class VMEnabledQSCIAlgorithm(QSCIAlgorithmBase):
    """QSCI algorithm with QURI VM analysis capabilities."""
    
    def __init__(
        self,
        base_algorithm: QSCIAlgorithmBase,
        architecture: Optional[ArchitectureInfo] = None
    ):
        """Initialize VM-enabled QSCI algorithm.
        
        Args:
            base_algorithm: Base QSCI algorithm
            architecture: Quantum architecture for analysis
        """
        # Copy properties from base algorithm
        super().__init__(
            base_algorithm.hamiltonian,
            base_algorithm.sampler,
            base_algorithm.num_states_pick_out
        )
        self.base_algorithm = base_algorithm
        self.vm_analyzer = VMCircuitAnalyzer(architecture)
    
    @property
    def name(self) -> str:
        return f"VM-enabled {self.base_algorithm.name}"
    
    def _create_qsci_implementation(self):
        return self.base_algorithm._create_qsci_implementation()
    
    def run(self, *args, **kwargs):
        """Run the base algorithm."""
        return self.base_algorithm.run(*args, **kwargs)
    
    def analyze(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        lowering_level: LoweringLevel = LoweringLevel.LogicalCircuit,
        **kwargs
    ) -> QSCIVMAnalysis:
        """Perform VM-enhanced analysis of QSCI algorithm."""
        circuit_estimates = {}
        
        for i, state in enumerate(input_states):
            circuit_name = f"input_circuit_{i}"
            
            if lowering_level == LoweringLevel.LogicalCircuit:
                estimate = self.vm_analyzer.analyze_logical_circuit(state.circuit)
            elif lowering_level == LoweringLevel.ArchLogicalCircuit:
                estimate = self.vm_analyzer.analyze_arch_logical_circuit(state.circuit)
            else:
                # For higher levels, use ArchLogicalCircuit analysis
                estimate = self.vm_analyzer.analyze_arch_logical_circuit(state.circuit)
            
            circuit_estimates[circuit_name] = estimate
        
        return QSCIVMAnalysis(
            lowering_level,
            circuit_estimates,
            total_shots,
            self.vm_analyzer.architecture
        )


def create_vm_enabled_algorithm(
    base_algorithm: QSCIAlgorithmBase,
    architecture_name: str = "STAR"
) -> VMEnabledQSCIAlgorithm:
    """Create VM-enabled QSCI algorithm with specified architecture.
    
    Args:
        base_algorithm: Base QSCI algorithm
        architecture_name: Name of quantum architecture ("STAR", etc.)
        
    Returns:
        VM-enabled QSCI algorithm
    """
    if architecture_name.upper() == "STAR":
        analyzer = VMCircuitAnalyzer()  # Uses STAR by default
    else:
        # For other architectures, create custom ArchitectureInfo
        architecture = ArchitectureInfo(
            name=architecture_name,
            num_physical_qubits=100,
            connectivity="custom",
            gate_fidelities={"single": 0.999, "two_qubit": 0.99, "measurement": 0.95},
            decoherence_times={"T1": 100.0, "T2": 50.0},
            gate_times={"single": 0.1, "two_qubit": 0.5, "measurement": 1.0}
        )
        analyzer = VMCircuitAnalyzer(architecture)
    
    return VMEnabledQSCIAlgorithm(base_algorithm, analyzer.architecture)


# Example usage functions
def demonstrate_star_architecture_analysis():
    """Demonstrate QSCI analysis with STAR architecture."""
    from .qsci_algo_interface import create_qsci_algorithm, QSCIVariant
    from quri_parts.core.operator import pauli_label
    
    # Create example Hamiltonian (H2 molecule)
    hamiltonian = pauli_label("Z0 Z1") + 0.5 * pauli_label("X0 X1")
    
    # Create TE-QSCI algorithm
    base_algorithm = create_qsci_algorithm(
        QSCIVariant.SINGLE_TIME_TE,
        hamiltonian,
        evolution_time=1.0
    )
    
    # Enable VM analysis
    vm_algorithm = create_vm_enabled_algorithm(base_algorithm, "STAR")
    
    # Create example input state (would be Hartree-Fock in practice)
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.core.state import CircuitQuantumState
    
    circuit = QuantumCircuit(2)
    circuit.add_X_gate(0)  # Simple example state
    input_state = CircuitQuantumState(2, circuit)
    
    # Analyze at different levels
    logical_analysis = vm_algorithm.analyze(
        [input_state], 1000, LoweringLevel.LogicalCircuit
    )
    
    arch_analysis = vm_algorithm.analyze(
        [input_state], 1000, LoweringLevel.ArchLogicalCircuit
    )
    
    print(f"Logical circuit analysis:")
    print(f"  Total latency: {logical_analysis.total_latency}")
    print(f"  Max qubits: {logical_analysis.max_physical_qubit_count}")
    
    print(f"ArchLogical circuit analysis:")
    print(f"  Total latency: {arch_analysis.total_latency}")
    print(f"  Max qubits: {arch_analysis.max_physical_qubit_count}")
    print(f"  SWAP overhead: {arch_analysis.total_swap_overhead}")
    print(f"  Average fidelity: {arch_analysis.average_fidelity:.4f}")


class VMSampler:
    """VM-based sampler implementing quri-parts ConcurrentSampler interface."""
    
    def __init__(self, vm_instance):
        """Initialize VM sampler.
        
        Args:
            vm_instance: Pre-configured QURI VM instance
        """
        self.vm = vm_instance
    
    def __call__(self, circuit_shot_pairs):
        """Sample from circuits using VM.sample() interface.
        
        Args:
            circuit_shot_pairs: List of (circuit, shots) tuples
            
        Returns:
            List of MeasurementCounts for each circuit
        """
        results = []
        
        for circuit, shots in circuit_shot_pairs:
            try:
                # Delegate directly to VM - all error modeling is VM-internal
                result = self.vm.sample(circuit, shots)
                # VM.sample() should return Counter format expected by quri-parts
                results.append(result)
            except Exception as e:
                raise RuntimeError(f"VM sampling failed: {e}. Ensure QURI VM is properly configured.") from e
        
        return results



class VMFactory:
    """Factory for creating pre-configured QURI VM instances."""
    
    @staticmethod
    def create_ideal_vm() -> 'VM':
        """Create ideal VM for LogicalCircuit level sampling.
        
        Returns:
            VM instance configured for ideal quantum computation
        """
        # Create abstract ideal VM (no device-specific constraints) 
        # Following tutorial: VM() creates abstract VM
        return VM()
    
    @staticmethod
    def create_star_vm(error_rate: float = 0.0) -> 'VM':
        """Create STAR architecture VM with specified error rate.
        
        Args:
            error_rate: Physical error rate (0.0 to 1.0)
            
        Returns:
            VM instance configured with STAR device properties for given error rate
        """
        # Following tutorial: Create STAR VM with specific physical error rate
        # QEC cycle time (1.0 microsecond is standard)
        qec_cycle = TimeValue(value=1.0, unit=TimeUnit.MICROSECOND)
        
        # Generate STAR device property with physical error rate
        # Following tutorial pattern for error rate configuration
        device_property = star_device.generate_device_property(
            qubit_count=16,  # Sufficient for H6 (12 qubits) + overhead
            code_distance=7,  # Standard distance for fault tolerance
            qec_cycle=qec_cycle,
            physical_error_rate=error_rate  # This is the key parameter from tutorial
        )
        
        return VM.from_device_prop(device_property)
    
    @staticmethod
    def create_star_vm_for_h6(error_rate: float = 0.0) -> 'VM':
        """Create STAR VM optimized for H6 TE-QSCI calculations.
        
        Args:
            error_rate: Physical error rate for STAR architecture
            
        Returns:
            VM instance optimized for H6 molecular calculations
        """
        return VMFactory.create_star_vm(error_rate)



class QURIVMInterface:
    """Interface to QURI VM for TE-QSCI analysis and sampling."""
    
    def __init__(self, logical_vm=None, arch_vm=None):
        """Initialize QURI VM interface with pre-configured VM instances.
        
        Args:
            logical_vm: Pre-configured VM for LogicalCircuit level (ideal)
            arch_vm: Pre-configured VM for ArchLogicalCircuit level (with architecture)
        """
        if logical_vm is None:
            logical_vm = VMFactory.create_ideal_vm()
        if arch_vm is None:
            arch_vm = VMFactory.create_star_vm(0.0)  # Default to no error
            
        self.logical_vm = logical_vm
        self.arch_vm = arch_vm
    
    def analyze_circuit_at_level(
        self,
        circuit: NonParametricQuantumCircuit,
        level: LoweringLevel
    ) -> Dict[str, Any]:
        """Analyze circuit at specified lowering level.
        
        Args:
            circuit: Circuit to analyze
            level: Analysis level (LogicalCircuit or ArchLogicalCircuit)
            
        Returns:
            Analysis results dictionary
        """
        if level == LoweringLevel.LogicalCircuit:
            return self._analyze_logical_level(circuit)
        elif level == LoweringLevel.ArchLogicalCircuit:
            return self._analyze_arch_level(circuit)
        else:
            # For higher levels, default to arch analysis
            return self._analyze_arch_level(circuit)
    
    def _analyze_logical_level(self, circuit: NonParametricQuantumCircuit) -> Dict[str, Any]:
        """Analyze at LogicalCircuit level (ideal)."""
        analysis = self.logical_vm.analyze(circuit)
        return analysis
    
    def _analyze_arch_level(self, circuit: NonParametricQuantumCircuit) -> Dict[str, Any]:
        """Analyze at ArchLogicalCircuit level (with architecture)."""
        analysis = self.arch_vm.analyze(circuit)
        return analysis
    
    def create_logical_sampler(self):
        """Create sampler for LogicalCircuit level (ideal VM).
        
        Returns:
            VMSampler instance for ideal quantum sampling
        """
        return VMSampler(self.logical_vm)
    
    def create_arch_sampler(self):
        """Create sampler for ArchLogicalCircuit level (STAR architecture).
        
        Returns:
            VMSampler instance for STAR architecture sampling
        """
        return VMSampler(self.arch_vm)
    
    def create_sampler(self, level: LoweringLevel):
        """Create sampler for specified level.
        
        Args:
            level: LoweringLevel for sampling (LogicalCircuit or ArchLogicalCircuit)
            
        Returns:
            VMSampler instance for the specified level
        """
        if level == LoweringLevel.LogicalCircuit:
            return self.create_logical_sampler()
        elif level == LoweringLevel.ArchLogicalCircuit:
            return self.create_arch_sampler()
        else:
            # Default to arch level for higher levels
            return self.create_arch_sampler()


def create_star_vm_for_h6(error_rate: float = 0.0) -> QURIVMInterface:
    """Create STAR architecture VM interface optimized for H6 TE-QSCI.
    
    Args:
        error_rate: Physical error rate for STAR architecture (0.0 for ideal)
    """
    logical_vm = VMFactory.create_ideal_vm()
    arch_vm = VMFactory.create_star_vm_for_h6(error_rate)
    return QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)


if __name__ == "__main__":
    demonstrate_star_architecture_analysis()