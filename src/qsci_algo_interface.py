"""
QSCI and TE-QSCI algorithm implementations compatible with quri-algo interface.

This module provides algorithm objects that conform to quri-algo's algorithm interface,
enabling seamless integration with QURI SDK ecosystem.
"""

from abc import abstractmethod
from typing import Any, Dict, List, Mapping, Optional, Sequence, TypeVar
import time
import numpy as np

from quri_parts.core.operator import Operator
from quri_parts.core.sampling import ConcurrentSampler
from quri_parts.core.state import CircuitQuantumState
from quri_parts.circuit import NonParametricQuantumCircuit
from quri_parts.backend.units import TimeValue, TimeUnit

# Import the algorithm interface types (based on the fetched interface)
# We'll create our own minimal version since we couldn't fetch the exact interface
from abc import ABC, abstractmethod
from enum import Enum

class LoweringLevel(Enum):
    """Analysis lowering levels for quantum circuits."""
    LogicalCircuit = 0
    ArchLogicalCircuit = 1
    ArchInstruction = 2
    DeviceInstruction = 3

class AlgorithmResult(ABC):
    """Base class for algorithm results."""
    
    def __init__(self, algorithm: "Algorithm", elapsed_time: Optional[float] = None):
        self.algorithm = algorithm
        self.elapsed_time = elapsed_time
    
    @property
    def name(self) -> str:
        """The name of the algorithm."""
        return self.algorithm.name

class Algorithm(ABC):
    """Base class for all algorithms."""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """The name of the algorithm."""
        pass
    
    @abstractmethod
    def run(self, *args: Any, **kwargs: Any) -> AlgorithmResult:
        """Run the algorithm itself."""
        pass
    
    def __str__(self) -> str:
        """Basic information about the algorithm should be returned."""
        return self.name

class Analysis(ABC):
    """Analysis of algorithm resource requirements."""
    
    def __init__(
        self,
        lowering_level: LoweringLevel,
        circuit_gate_count: Mapping[Any, int],
        circuit_depth: Mapping[Any, int],
        circuit_latency: Mapping[Any, TimeValue | None],
        circuit_execution_count: Mapping[Any, int],
        circuit_fidelities: Mapping[Any, float | None],
        circuit_qubit_count: Mapping[Any, int],
    ) -> None:
        self.lowering_level = lowering_level
        self.circuit_gate_count = circuit_gate_count
        self.circuit_depth = circuit_depth
        self.circuit_latency = circuit_latency
        self.circuit_execution_count = circuit_execution_count
        self.circuit_fidelities = circuit_fidelities
        self.circuit_qubit_count = circuit_qubit_count
    
    @property
    @abstractmethod
    def total_latency(self) -> TimeValue:
        """Total latency of the circuit is algorithm dependent."""
        pass
    
    @property
    @abstractmethod
    def max_physical_qubit_count(self) -> int:
        """Maximum physical qubit count is algorithm dependent."""
        pass

class CircuitAnalysisMixin:
    """Mixin for algorithms that support circuit analysis."""
    
    @abstractmethod
    def analyze(self, *args: Any, **kwargs: Any) -> Analysis:
        """The quantum resource analysis of the algorithm."""
        pass

class QuantumAlgorithm(Algorithm, CircuitAnalysisMixin, ABC):
    """Base class for quantum algorithms with circuit analysis."""
    pass

# Import our QSCI implementations
from .qsci_algorithms import (
    QSCIResult, QSCIVariant, VanillaQSCI, TimeEvolvedQSCI
)

T = TypeVar("T")

class QSCIAlgorithmResult(AlgorithmResult):
    """Result wrapper for QSCI algorithms compatible with quri-algo interface."""
    
    def __init__(
        self, 
        algorithm: "QSCIAlgorithmBase", 
        qsci_result: QSCIResult,
        elapsed_time: Optional[float] = None
    ):
        super().__init__(algorithm, elapsed_time or qsci_result.execution_time)
        self.qsci_result = qsci_result
    
    @property
    def ground_state_energy(self) -> float:
        """Ground state energy from QSCI calculation."""
        return self.qsci_result.ground_state_energy
    
    @property
    def eigenvalues(self) -> Sequence[float]:
        """All computed eigenvalues."""
        return self.qsci_result.eigenvalues
    
    @property
    def subspace_dimension(self) -> int:
        """Dimension of the selected subspace."""
        return self.qsci_result.subspace_dimension


class QSCIAnalysis(Analysis):
    """Analysis implementation for QSCI algorithms."""
    
    def __init__(
        self,
        lowering_level: LoweringLevel,
        qsci_result: QSCIResult,
        circuit_analysis: Optional[Dict[str, Any]] = None
    ):
        # Initialize with QSCI-specific circuit analysis
        circuit_gate_count = circuit_analysis.get("gate_count", {}) if circuit_analysis else {}
        circuit_depth = circuit_analysis.get("depth", {}) if circuit_analysis else {}
        circuit_latency = circuit_analysis.get("latency", {}) if circuit_analysis else {}
        circuit_execution_count = {"sampling": qsci_result.total_shots}
        circuit_fidelities = circuit_analysis.get("fidelities", {}) if circuit_analysis else {}
        circuit_qubit_count = circuit_analysis.get("qubit_count", {}) if circuit_analysis else {}
        
        super().__init__(
            lowering_level,
            circuit_gate_count,
            circuit_depth,
            circuit_latency,
            circuit_execution_count,
            circuit_fidelities,
            circuit_qubit_count
        )
        self.qsci_result = qsci_result
    
    @property
    def total_latency(self) -> TimeValue:
        """Total latency for QSCI execution."""
        # For QSCI, latency is dominated by sampling time
        if self.circuit_latency:
            total_us = sum(
                (lat.value if lat else 0) 
                for lat in self.circuit_latency.values()
            )
            return TimeValue(total_us, TimeUnit.MICROSECOND)
        return TimeValue(0, TimeUnit.MICROSECOND)
    
    @property
    def max_physical_qubit_count(self) -> int:
        """Maximum physical qubit count required."""
        if self.circuit_qubit_count:
            return max(self.circuit_qubit_count.values())
        # Estimate from selected states
        if self.qsci_result.selected_states:
            return self.qsci_result.selected_states[0].qubit_count
        return 0


class QSCIAlgorithmBase(QuantumAlgorithm):
    """Base class for QSCI algorithms compatible with quri-algo interface."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None
    ):
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        self.num_states_pick_out = num_states_pick_out
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Algorithm name."""
        pass
    
    @abstractmethod
    def _create_qsci_implementation(self):
        """Create the underlying QSCI implementation."""
        pass
    
    def analyze(
        self, 
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        lowering_level: LoweringLevel = LoweringLevel.LogicalCircuit,
        **kwargs
    ) -> QSCIAnalysis:
        """Analyze quantum resource requirements for QSCI algorithm."""
        # For analysis, we can estimate without actually running
        circuit_analysis = {}
        
        if input_states:
            # Analyze input circuits
            total_gates = sum(len(state.circuit) for state in input_states)
            max_depth = max(state.circuit.depth for state in input_states) if input_states else 0
            qubit_count = input_states[0].qubit_count if input_states else 0
            
            circuit_analysis = {
                "gate_count": {"total": total_gates},
                "depth": {"max": max_depth},
                "qubit_count": {"required": qubit_count},
                "latency": {},  # Would need backend-specific estimation
                "fidelities": {}  # Would need noise model
            }
        
        # Create dummy result for analysis
        dummy_result = QSCIResult(
            eigenvalues=[0.0],
            eigenstates=[],
            selected_states=[],
            subspace_dimension=self.num_states_pick_out or 100,
            total_shots=total_shots,
            algorithm_variant=QSCIVariant.VANILLA
        )
        
        return QSCIAnalysis(lowering_level, dummy_result, circuit_analysis)


class VanillaQSCIAlgorithm(QSCIAlgorithmBase):
    """Vanilla QSCI algorithm compatible with quri-algo interface."""
    
    @property
    def name(self) -> str:
        return "Vanilla QSCI"
    
    def _create_qsci_implementation(self):
        return VanillaQSCI(
            self.hamiltonian, 
            self.sampler, 
            self.num_states_pick_out
        )
    
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run vanilla QSCI algorithm."""
        start_time = time.time()
        
        qsci_impl = self._create_qsci_implementation()
        qsci_result = qsci_impl.run(input_states, total_shots, start_time=start_time, **kwargs)
        
        elapsed_time = time.time() - start_time
        return QSCIAlgorithmResult(self, qsci_result, elapsed_time)


class SingleTimeTeQSCIAlgorithm(QSCIAlgorithmBase):
    """Single-time TE-QSCI algorithm compatible with quri-algo interface."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_time: float,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        trotter_steps: Optional[int] = None,
        time_evolution_method: str = "trotter"
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
        self.trotter_steps = trotter_steps
        self.time_evolution_method = time_evolution_method
    
    @property
    def name(self) -> str:
        return f"Single-time TE-QSCI (t={self.evolution_time}, method={self.time_evolution_method})"
    
    def _create_qsci_implementation(self):
        return TimeEvolvedQSCI(
            self.hamiltonian,
            self.sampler,
            self.num_states_pick_out,
            self.time_evolution_method
        )
    
    def run(
        self,
        initial_state: CircuitQuantumState,
        total_shots: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run single-time TE-QSCI algorithm."""
        start_time = time.time()
        
        # Use trotter_steps from constructor, or from kwargs if not set in constructor
        trotter_steps = self.trotter_steps if self.trotter_steps is not None else kwargs.get("trotter_steps")
        
        # Remove trotter_steps from kwargs to avoid duplicate argument
        filtered_kwargs = {k: v for k, v in kwargs.items() if k != "trotter_steps"}
        
        te_qsci_impl = self._create_qsci_implementation()
        qsci_result = te_qsci_impl.run_single_time(
            initial_state,
            self.evolution_time,
            total_shots,
            trotter_steps,
            start_time=start_time,
            **filtered_kwargs
        )
        
        elapsed_time = time.time() - start_time
        return QSCIAlgorithmResult(self, qsci_result, elapsed_time)


class TimeAverageTeQSCIAlgorithm(QSCIAlgorithmBase):
    """Time-average TE-QSCI algorithm compatible with quri-algo interface."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_times: Sequence[float],
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        trotter_steps: Optional[int] = None,
        time_evolution_method: str = "trotter"
    ):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_times = evolution_times
        self.trotter_steps = trotter_steps
        self.time_evolution_method = time_evolution_method
    
    @property
    def name(self) -> str:
        return f"Time-average TE-QSCI ({len(self.evolution_times)} times, method={self.time_evolution_method})"
    
    def _create_qsci_implementation(self):
        return TimeEvolvedQSCI(
            self.hamiltonian,
            self.sampler,
            self.num_states_pick_out,
            self.time_evolution_method
        )
    
    def run(
        self,
        initial_state: CircuitQuantumState,
        shots_per_time: int,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run time-average TE-QSCI algorithm."""
        start_time = time.time()
        
        # Use trotter_steps from constructor, or from kwargs if not set in constructor
        trotter_steps = self.trotter_steps if self.trotter_steps is not None else kwargs.get("trotter_steps")
        
        # Remove trotter_steps from kwargs to avoid duplicate argument
        filtered_kwargs = {k: v for k, v in kwargs.items() if k != "trotter_steps"}
        
        te_qsci_impl = self._create_qsci_implementation()
        qsci_result = te_qsci_impl.run_time_average(
            initial_state,
            self.evolution_times,
            shots_per_time,
            trotter_steps,
            start_time=start_time,
            **filtered_kwargs
        )
        
        elapsed_time = time.time() - start_time
        return QSCIAlgorithmResult(self, qsci_result, elapsed_time)


class StateVectorTeQSCIAlgorithm(QSCIAlgorithmBase):
    """State vector TE-QSCI algorithm compatible with quri-algo interface."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        evolution_time: float,
        num_states_pick_out: Optional[int] = None,
        time_evolution_method: str = "exact"
    ):
        super().__init__(hamiltonian, None, num_states_pick_out)  # No sampler needed
        self.evolution_time = evolution_time
        self.time_evolution_method = time_evolution_method
    
    @property
    def name(self) -> str:
        return f"State vector TE-QSCI (t={self.evolution_time}, method={self.time_evolution_method})"
    
    def _create_qsci_implementation(self):
        return TimeEvolvedQSCI(
            self.hamiltonian,
            None,  # No sampler for state vector
            self.num_states_pick_out,
            self.time_evolution_method
        )
    
    def run(
        self,
        initial_state: CircuitQuantumState,
        num_eigenstates: int = 1,
        **kwargs
    ) -> QSCIAlgorithmResult:
        """Run state vector TE-QSCI algorithm."""
        start_time = time.time()
        
        te_qsci_impl = self._create_qsci_implementation()
        qsci_result = te_qsci_impl.run_state_vector(
            initial_state,
            self.evolution_time,
            num_eigenstates,
            start_time=start_time,
            **kwargs
        )
        
        elapsed_time = time.time() - start_time
        return QSCIAlgorithmResult(self, qsci_result, elapsed_time)


# Factory function for easy algorithm creation
def create_qsci_algorithm(
    variant: QSCIVariant,
    hamiltonian: Operator,
    sampler: Optional[ConcurrentSampler] = None,
    num_states_pick_out: Optional[int] = None,
    **variant_kwargs
) -> QSCIAlgorithmBase:
    """Factory function to create QSCI algorithm variants.
    
    Args:
        variant: QSCI algorithm variant
        hamiltonian: Target Hamiltonian
        sampler: Quantum sampler
        num_states_pick_out: Number of states to select
        **variant_kwargs: Variant-specific parameters
        
    Returns:
        QSCI algorithm instance
    """
    if variant == QSCIVariant.VANILLA:
        return VanillaQSCIAlgorithm(hamiltonian, sampler, num_states_pick_out)
    
    elif variant == QSCIVariant.SINGLE_TIME_TE:
        evolution_time = variant_kwargs.get("evolution_time", 1.0)
        trotter_steps = variant_kwargs.get("trotter_steps")
        time_evolution_method = variant_kwargs.get("time_evolution_method", "trotter")
        return SingleTimeTeQSCIAlgorithm(
            hamiltonian, evolution_time, sampler, num_states_pick_out, trotter_steps, time_evolution_method
        )
    
    elif variant == QSCIVariant.TIME_AVERAGE_TE:
        evolution_times = variant_kwargs.get("evolution_times", [0.5, 1.0, 1.5])
        trotter_steps = variant_kwargs.get("trotter_steps")
        time_evolution_method = variant_kwargs.get("time_evolution_method", "trotter")
        return TimeAverageTeQSCIAlgorithm(
            hamiltonian, evolution_times, sampler, num_states_pick_out, trotter_steps, time_evolution_method
        )
    
    elif variant == QSCIVariant.STATE_VECTOR:
        evolution_time = variant_kwargs.get("evolution_time", 1.0)
        time_evolution_method = variant_kwargs.get("time_evolution_method", "exact")
        return StateVectorTeQSCIAlgorithm(hamiltonian, evolution_time, num_states_pick_out, time_evolution_method)
    
    else:
        raise ValueError(f"Unknown QSCI variant: {variant}")