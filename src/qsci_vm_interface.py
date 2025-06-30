"""
VM-Aware TE-QSCI Interface
=========================

This module provides enhanced interfaces for running TE-QSCI algorithms with integrated
QURI VM analysis capabilities. It extends the base QSCI algorithms to capture and analyze
quantum circuits at both LogicalCircuit and ArchLogicalCircuit levels.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import time

from quri_parts.circuit import NonParametricQuantumCircuit, QuantumCircuit
from quri_parts.backend.units import TimeValue, TimeUnit
from quri_parts.core.state import (
    CircuitQuantumState,
    ComputationalBasisState,
    ComputationalBasisSuperposition
)
from quri_parts.core.operator import Operator
from quri_parts.backend.units import TimeValue, TimeUnit

from .qsci_algorithms import QSCIBase, QSCIResult, QSCIVariant, SingleTimeTE_QSCI
from .qsci_algo_interface import LoweringLevel, QSCIAlgorithmBase
from .qsci_vm_analysis import (
    QURIVMInterface,
    VMCircuitAnalyzer,
    QSCIVMAnalysis,
    ArchitectureInfo,
    VMSampler
)


@dataclass
class VMAnalysisResult:
    """Result of VM-enhanced circuit analysis."""
    logical_analysis: Dict[str, Any]
    arch_analysis: Dict[str, Any]
    circuit_description: str
    analysis_time: float


@dataclass
class TEQSCIVMResult:
    """Enhanced TE-QSCI result with VM analysis."""
    # Standard QSCI results
    qsci_result: QSCIResult
    
    # VM Analysis results
    time_evolution_analysis: VMAnalysisResult
    sampling_circuits_analysis: List[VMAnalysisResult]
    
    # Performance metrics
    total_logical_time: TimeValue
    total_arch_time: TimeValue
    architecture_overhead: float
    fidelity_impact: float
    
    # Algorithm metadata
    evolution_time: float
    trotter_steps: int
    architecture_name: str


class VMAwareTimeEvolutionCircuit:
    """Wrapper for time evolution circuits with VM analysis capabilities."""
    
    def __init__(
        self,
        base_circuit: NonParametricQuantumCircuit,
        evolution_time: float,
        trotter_steps: int,
        vm_interface: QURIVMInterface,
        hamiltonian: Optional[Operator] = None
    ):
        """Initialize VM-aware time evolution circuit.
        
        Args:
            base_circuit: Base quantum circuit
            evolution_time: Time evolution parameter
            trotter_steps: Number of Trotter decomposition steps
            vm_interface: QURI VM interface for analysis
            hamiltonian: Target Hamiltonian for authentic analysis
        """
        self.base_circuit = base_circuit
        self.evolution_time = evolution_time
        self.trotter_steps = trotter_steps
        self.vm_interface = vm_interface
        self.hamiltonian = hamiltonian
        
        # Cache for analysis results
        self._logical_analysis: Optional[Dict] = None
        self._arch_analysis: Optional[Dict] = None
    
    @property
    def circuit(self) -> NonParametricQuantumCircuit:
        """Get the underlying circuit."""
        return self.base_circuit
    
    def analyze_at_level(self, level: LoweringLevel) -> Dict[str, Any]:
        """Analyze circuit at specified level with caching."""
        if level == LoweringLevel.LogicalCircuit:
            if self._logical_analysis is None:
                self._logical_analysis = self.vm_interface.analyze_circuit_at_level(
                    self.base_circuit, level
                )
            return self._logical_analysis
        else:
            if self._arch_analysis is None:
                self._arch_analysis = self.vm_interface.analyze_circuit_at_level(
                    self.base_circuit, level
                )
            return self._arch_analysis
    
    def get_analysis_summary(self) -> VMAnalysisResult:
        """Get comprehensive analysis summary."""
        start_time = time.time()
        
        logical_analysis = self.analyze_at_level(LoweringLevel.LogicalCircuit)
        arch_analysis = self.analyze_at_level(LoweringLevel.ArchLogicalCircuit)
        
        analysis_time = time.time() - start_time
        
        circuit_description = f"Time evolution (t={self.evolution_time}, steps={self.trotter_steps})"
        
        return VMAnalysisResult(
            logical_analysis=logical_analysis,
            arch_analysis=arch_analysis,
            circuit_description=circuit_description,
            analysis_time=analysis_time
        )


class VMAwareSingleTimeTE_QSCI(SingleTimeTE_QSCI):
    """TE-QSCI with integrated VM analysis capabilities."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_time: float,
        sampler=None,
        num_states_pick_out: Optional[int] = None,
        trotter_steps: int = 10,
        vm_interface: Optional[QURIVMInterface] = None,
        sampling_level: LoweringLevel = LoweringLevel.ArchLogicalCircuit,
        use_vm_sampling: bool = True
    ):
        """Initialize VM-aware TE-QSCI.
        
        Args:
            hamiltonian: Target Hamiltonian
            evolution_time: Time evolution parameter
            sampler: Quantum sampler (used as fallback if use_vm_sampling=False)
            num_states_pick_out: Number of states for QSCI subspace
            trotter_steps: Number of Trotter steps for time evolution
            vm_interface: Pre-configured QURIVMInterface with appropriate VMs
            sampling_level: LoweringLevel for VM sampling (LogicalCircuit or ArchLogicalCircuit)
            use_vm_sampling: Whether to use VM sampling or fallback to provided sampler
        """
        # Use provided VM interface or create default one
        if vm_interface is None:
            # Create default VM interface with ideal and zero-error STAR VMs
            vm_interface = QURIVMInterface()
        
        self.vm_interface = vm_interface
        
        # Set up sampling strategy
        self.sampling_level = sampling_level
        self.use_vm_sampling = use_vm_sampling
        self.fallback_sampler = sampler
        
        # Create appropriate sampler
        if use_vm_sampling:
            vm_sampler = self.vm_interface.create_sampler(sampling_level)
        else:
            vm_sampler = sampler
        
        super().__init__(
            hamiltonian=hamiltonian,
            sampler=vm_sampler,
            num_states_pick_out=num_states_pick_out,
            evolution_time=evolution_time
        )
        
        self.trotter_steps = trotter_steps
        self.architecture = "STAR"  # Default architecture for VM analysis
        
        # Storage for VM analysis
        self.circuit_analyses: List[VMAnalysisResult] = []
        self.time_evolution_analysis: Optional[VMAnalysisResult] = None
    
    @property
    def name(self) -> str:
        """Algorithm name with VM capabilities."""
        sampling_mode = "VM-sampling" if self.use_vm_sampling else "Qulacs-sampling"
        level_name = "Ideal" if self.sampling_level == LoweringLevel.LogicalCircuit else "STAR"
        return f"VM-aware TE-QSCI ({level_name} {sampling_mode})"
    
    def _create_time_evolution_circuit(
        self, 
        initial_state: CircuitQuantumState
    ) -> VMAwareTimeEvolutionCircuit:
        """Create time evolution circuit with VM analysis using real Trotter decomposition.
        
        Args:
            initial_state: Initial quantum state
            
        Returns:
            VM-aware time evolution circuit
        """
        # Create authentic time evolution circuit using quri-algo
        n_qubits = initial_state.qubit_count
        circuit = QuantumCircuit(n_qubits)
        
        # Add initial state preparation
        if hasattr(initial_state, 'circuit'):
            for gate in initial_state.circuit.gates:
                circuit.add_gate(gate)
        
        # Generate real time evolution circuit using quri-algo
        try:
            time_evo_circuit = self._generate_authentic_time_evolution_circuit(n_qubits)
            
            # Combine initial state with time evolution
            for gate in time_evo_circuit.gates:
                circuit.add_gate(gate)
                
        except Exception as e:
            print(f"   Warning: Failed to generate authentic time evolution circuit: {e}")
            print(f"   Falling back to simplified Trotter approximation")
            
            # Fallback to improved (but still simplified) Trotter decomposition
            self._add_simplified_trotter_steps(circuit, n_qubits)
        
        return VMAwareTimeEvolutionCircuit(
            base_circuit=circuit,
            evolution_time=self.evolution_time,
            trotter_steps=self.trotter_steps,
            vm_interface=self.vm_interface,
            hamiltonian=self.hamiltonian  # Pass Hamiltonian for authentic analysis
        )
    
    def _generate_authentic_time_evolution_circuit(self, n_qubits: int) -> QuantumCircuit:
        """Generate authentic time evolution circuit using quri-algo Trotter decomposition.
        
        Args:
            n_qubits: Number of qubits
            
        Returns:
            Quantum circuit implementing real Trotter decomposition
        """
        try:
            from quri_algo.circuit.time_evolution.trotter_time_evo import TrotterTimeEvolutionCircuitFactory
            from quri_algo.problem.operators.hamiltonian import QubitHamiltonian
            
            # Convert quri-parts Operator to quri-algo QubitHamiltonian
            qubit_hamiltonian = self._convert_operator_to_qubit_hamiltonian(self.hamiltonian, n_qubits)
            
            # Create Trotter time evolution factory
            factory = TrotterTimeEvolutionCircuitFactory(
                qubit_hamiltonian, 
                n_trotter=self.trotter_steps
            )
            
            # Generate time evolution circuit
            time_evo_circuit = factory(self.evolution_time)
            
            return time_evo_circuit
            
        except ImportError as e:
            raise ImportError(f"quri-algo not available for authentic time evolution: {e}")
        except Exception as e:
            raise RuntimeError(f"Failed to create authentic time evolution circuit: {e}")
    
    def _convert_operator_to_qubit_hamiltonian(self, operator: Operator, n_qubits: int):
        """Convert quri-parts Operator to quri-algo QubitHamiltonian.
        
        Args:
            operator: quri-parts Operator
            n_qubits: Number of qubits
            
        Returns:
            quri-algo QubitHamiltonian
        """
        try:
            from quri_algo.problem.operators.hamiltonian import QubitHamiltonian
            
            # QubitHamiltonian constructor takes (n_qubit, quri_parts.Operator)
            # This is much simpler than trying to convert Pauli labels
            return QubitHamiltonian(n_qubits, operator)
            
        except Exception as e:
            raise RuntimeError(f"Failed to convert Operator to QubitHamiltonian: {e}")
    
    def _add_simplified_trotter_steps(self, circuit: QuantumCircuit, n_qubits: int):
        """Add simplified but improved Trotter steps as fallback.
        
        This is more realistic than the previous random implementation but still
        simplified compared to full Hamiltonian decomposition.
        
        Args:
            circuit: Quantum circuit to modify
            n_qubits: Number of qubits
        """
        dt = self.evolution_time / self.trotter_steps
        
        # Try to extract some information from the Hamiltonian for more realistic circuits
        try:
            # Get rough idea of Hamiltonian structure
            max_qubit = max(max(pauli.qubit_indices) if pauli.qubit_indices else 0 
                           for pauli in self.hamiltonian.keys())
            hamiltonian_terms = len(self.hamiltonian)
            
            print(f"   Using simplified Trotter with {hamiltonian_terms} Hamiltonian terms on {max_qubit+1} qubits")
            
        except Exception:
            max_qubit = n_qubits - 1
            hamiltonian_terms = n_qubits
        
        for step in range(self.trotter_steps):
            # Add single-qubit rotations proportional to evolution time
            for i in range(min(n_qubits, max_qubit + 1)):
                # Use smaller angles proportional to actual evolution time
                circuit.add_RZ_gate(i, dt * 0.1)  # More realistic than π/8
            
            # Add two-qubit interactions based on typical molecular Hamiltonians
            for i in range(0, min(n_qubits - 1, max_qubit), 2):
                circuit.add_CNOT_gate(i, i + 1)
                circuit.add_RZ_gate(i + 1, dt * 0.05)  # Smaller coupling terms
                circuit.add_CNOT_gate(i, i + 1)
            
            # Add additional entangling layer with reduced complexity
            if self.trotter_steps > 5:  # Only for deeper circuits
                for i in range(1, min(n_qubits - 1, max_qubit), 2):
                    circuit.add_CNOT_gate(i, i + 1)
                    circuit.add_RY_gate(i + 1, dt * 0.02)
                    circuit.add_CNOT_gate(i, i + 1)
    
    def run_with_vm_analysis(
        self, 
        initial_state: CircuitQuantumState,
        analyze_sampling: bool = True,
        num_shots: int = 1000
    ) -> TEQSCIVMResult:
        """Run TE-QSCI with comprehensive VM analysis.
        
        Args:
            initial_state: Initial quantum state
            analyze_sampling: Whether to analyze sampling circuits
            
        Returns:
            Enhanced TE-QSCI result with VM analysis
        """
        print(f" Running {self.name}")
        print(f"   • Evolution time: {self.evolution_time}")
        print(f"   • Trotter steps: {self.trotter_steps}")
        print(f"   • Architecture: {self.architecture}")
        print(f"   • Sampling mode: {'VM (' + self.sampling_level.name + ')' if self.use_vm_sampling else 'Qulacs (fallback)'}")
        
        # Validate sampling setup
        if self.use_vm_sampling and not isinstance(self.sampler, VMSampler):
            print(f"   Warning: Expected VMSampler but got {type(self.sampler)}")
            if self.fallback_sampler:
                print(f"   Falling back to Qulacs sampler...")
                self.sampler = self.fallback_sampler
                self.use_vm_sampling = False
        
        # Create and analyze time evolution circuit
        print("   Analyzing time evolution circuit...")
        te_circuit = self._create_time_evolution_circuit(initial_state)
        time_evolution_analysis = te_circuit.get_analysis_summary()
        self.time_evolution_analysis = time_evolution_analysis
        
        # Run standard TE-QSCI algorithm
        print("   Running TE-QSCI computation...")
        qsci_result = self.run([initial_state], num_shots)
        
        # Analyze sampling circuits if requested
        sampling_analyses = []
        if analyze_sampling:
            print("   Analyzing sampling circuits...")
            # For demonstration, create a few representative sampling circuits
            sampling_analyses = self._analyze_sampling_circuits(initial_state)
        
        # Calculate performance metrics
        logical_time = time_evolution_analysis.logical_analysis.latency or TimeValue(0, TimeUnit.MICROSECOND)
        arch_time = time_evolution_analysis.arch_analysis.latency or TimeValue(0, TimeUnit.MICROSECOND)
        
        total_logical_time = logical_time
        total_arch_time = arch_time
        
        # Add sampling time if analyzed
        for analysis in sampling_analyses:
            logical_latency = analysis.logical_analysis.latency
            arch_latency = analysis.arch_analysis.latency
            if logical_latency:
                total_logical_time = TimeValue(
                    total_logical_time.in_ns() + logical_latency.in_ns(),
                    TimeUnit.NANOSECOND
                )
            if arch_latency:
                total_arch_time = TimeValue(
                    total_arch_time.in_ns() + arch_latency.in_ns(),
                    TimeUnit.NANOSECOND
                )
        
        # Calculate overhead and fidelity impact
        if total_logical_time.in_ns() > 0:
            architecture_overhead = total_arch_time.in_ns() / total_logical_time.in_ns()
        else:
            architecture_overhead = 1.0  # No overhead if no time measured
        
        logical_fidelity = time_evolution_analysis.logical_analysis.fidelity or 1.0
        arch_fidelity = time_evolution_analysis.arch_analysis.fidelity or 1.0
        fidelity_impact = logical_fidelity - arch_fidelity
        
        print(f"   VM analysis complete!")
        print(f"      • Architecture overhead: {architecture_overhead:.1f}x")
        print(f"      • Fidelity impact: {fidelity_impact:.1%}")
        
        return TEQSCIVMResult(
            qsci_result=qsci_result,
            time_evolution_analysis=time_evolution_analysis,
            sampling_circuits_analysis=sampling_analyses,
            total_logical_time=total_logical_time,
            total_arch_time=total_arch_time,
            architecture_overhead=architecture_overhead,
            fidelity_impact=fidelity_impact,
            evolution_time=self.evolution_time,
            trotter_steps=self.trotter_steps,
            architecture_name=self.architecture
        )
    
    def _analyze_sampling_circuits(
        self, 
        initial_state: CircuitQuantumState,
        n_sample_circuits: int = 3
    ) -> List[VMAnalysisResult]:
        """Analyze representative sampling circuits.
        
        Args:
            initial_state: Initial quantum state
            n_sample_circuits: Number of sample circuits to analyze
            
        Returns:
            List of sampling circuit analyses
        """
        analyses = []
        
        for i in range(n_sample_circuits):
            # Create representative sampling circuit
            n_qubits = initial_state.qubit_count
            circuit = QuantumCircuit(n_qubits)
            
            # Add some typical measurement preparation gates
            for j in range(n_qubits):
                if np.random.random() < 0.3:  # 30% chance of rotation
                    circuit.add_RY_gate(j, np.pi / 4)
            
            # Add measurements (would be done by sampling in practice)
            # For analysis purposes, we just count the preparation overhead
            
            # Analyze this circuit
            start_time = time.time()
            logical_analysis = self.vm_interface.analyze_circuit_at_level(
                circuit, LoweringLevel.LogicalCircuit
            )
            arch_analysis = self.vm_interface.analyze_circuit_at_level(
                circuit, LoweringLevel.ArchLogicalCircuit
            )
            analysis_time = time.time() - start_time
            
            analysis = VMAnalysisResult(
                logical_analysis=logical_analysis,
                arch_analysis=arch_analysis,
                circuit_description=f"Sampling circuit {i+1}",
                analysis_time=analysis_time
            )
            
            analyses.append(analysis)
        
        return analyses


def create_vm_aware_te_qsci(
    hamiltonian: Operator,
    evolution_time: float = 1.5,
    vm_interface: Optional[QURIVMInterface] = None,
    sampling_level: LoweringLevel = LoweringLevel.ArchLogicalCircuit,
    use_vm_sampling: bool = True,
    **kwargs
) -> VMAwareSingleTimeTE_QSCI:
    """Create VM-aware TE-QSCI algorithm with VM sampling.
    
    Args:
        hamiltonian: Target Hamiltonian
        evolution_time: Time evolution parameter
        vm_interface: Pre-configured QURIVMInterface with appropriate VMs
        sampling_level: LoweringLevel for VM sampling (LogicalCircuit for ideal, ArchLogicalCircuit for STAR)
        use_vm_sampling: Whether to use VM sampling or fallback to Qulacs
        **kwargs: Additional arguments for TE-QSCI
        
    Returns:
        VM-aware TE-QSCI algorithm with VM sampling capabilities
    """
    return VMAwareSingleTimeTE_QSCI(
        hamiltonian=hamiltonian,
        evolution_time=evolution_time,
        vm_interface=vm_interface,
        sampling_level=sampling_level,
        use_vm_sampling=use_vm_sampling,
        **kwargs
    )


def demonstrate_h6_vm_analysis():
    """Demonstrate H6 TE-QSCI with VM analysis."""
    from quri_parts.core.operator import pauli_label
    
    # Create simplified H6 Hamiltonian for demonstration
    hamiltonian = Operator()
    
    # Add typical terms found in H6 molecular Hamiltonian
    n_qubits = 12  # H6 requires 12 qubits
    
    # Single-qubit terms
    for i in range(n_qubits):
        hamiltonian += 0.1 * pauli_label(f"Z{i}")
    
    # Two-qubit terms
    for i in range(0, n_qubits - 1, 2):
        hamiltonian += 0.05 * pauli_label(f"Z{i} Z{i+1}")
        hamiltonian += 0.02 * pauli_label(f"X{i} X{i+1}")
    
    # Create VM-aware TE-QSCI
    te_qsci = create_vm_aware_te_qsci(
        hamiltonian=hamiltonian,
        evolution_time=1.5,
        architecture="STAR",
        num_states_pick_out=100,
        trotter_steps=10
    )
    
    # Create simple initial state
    from quri_parts.circuit import QuantumCircuit
    from quri_parts.core.state import CircuitQuantumState
    
    circuit = QuantumCircuit(n_qubits)
    # Create HF-like state for H6 (6 electrons)
    for i in range(6):
        circuit.add_X_gate(i)
    
    initial_state = CircuitQuantumState(n_qubits, circuit)
    
    # Run with VM analysis
    result = te_qsci.run_with_vm_analysis(initial_state)
    
    print("\n VM Analysis Results:")
    print(f"   Architecture overhead: {result.architecture_overhead:.1f}x")
    print(f"   Fidelity impact: {result.fidelity_impact:.1%}")
    print(f"   Logical execution time: {result.total_logical_time}")
    print(f"   STAR execution time: {result.total_arch_time}")
    
    return result


if __name__ == "__main__":
    demonstrate_h6_vm_analysis()