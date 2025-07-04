"""
Implementation of QSCI and TE-QSCI algorithms using quri-parts framework.

This module provides implementations of:
1. Vanilla QSCI algorithm
2. Time-Evolved QSCI (TE-QSCI) variants
3. Algorithm interfaces compatible with quri-algo

Key algorithms:
- QSCI: Quantum Selected Configuration Interaction
- TE-QSCI: Time-Evolved QSCI (single-time and time-average variants)
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union
import numpy as np
import numpy.typing as npt

from quri_parts.core.operator import Operator, is_hermitian
from quri_parts.core.sampling import ConcurrentSampler, MeasurementCounts
from quri_parts.core.state import (
    CircuitQuantumState,
    ComputationalBasisState,
    ComputationalBasisSuperposition,
    GeneralCircuitQuantumState,
    QuantumState
)
from quri_parts.circuit import NonParametricQuantumCircuit, QuantumCircuit
from quri_parts.algo.optimizer import OptimizerState
from quri_parts.backend.units import TimeValue


class QSCIVariant(Enum):
    """Enumeration of QSCI algorithm variants."""
    VANILLA = "vanilla"
    SINGLE_TIME_TE = "single_time_te"
    TIME_AVERAGE_TE = "time_average_te"
    STATE_VECTOR = "state_vector"


@dataclass
class QSCIResult:
    """Result of QSCI algorithm execution."""
    eigenvalues: Sequence[float]
    eigenstates: Sequence[ComputationalBasisSuperposition]
    selected_states: Sequence[ComputationalBasisState]
    subspace_dimension: int
    total_shots: int
    algorithm_variant: QSCIVariant
    execution_time: Optional[float] = None
    computational_basis_probabilities: Optional[Dict[int, float]] = None
    
    @property
    def ground_state_energy(self) -> float:
        """Ground state energy (lowest eigenvalue)."""
        return min(self.eigenvalues)
    
    @property
    def ground_state(self) -> ComputationalBasisSuperposition:
        """Ground state (eigenstate corresponding to lowest eigenvalue)."""
        min_idx = np.argmin(self.eigenvalues)
        return self.eigenstates[min_idx]
    
    @property
    def ground_state_probability(self) -> float:
        """Probability of measuring the ground state in computational basis."""
        if self.computational_basis_probabilities is None:
            return 0.0
        ground_state = self.ground_state
        if hasattr(ground_state, 'bits'):
            return self.computational_basis_probabilities.get(ground_state.bits, 0.0)
        return 0.0


class QSCIBase(ABC):
    """Base class for QSCI algorithms."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None
    ):
        """Initialize QSCI algorithm.
        
        Args:
            hamiltonian: Target Hamiltonian to diagonalize
            sampler: Quantum sampler for measurement
            num_states_pick_out: Number of states to select for subspace
        """
        if not is_hermitian(hamiltonian):
            raise ValueError("Hamiltonian must be Hermitian")
        
        self.hamiltonian = hamiltonian
        self.sampler = sampler
        self.num_states_pick_out = num_states_pick_out
    
    @abstractmethod
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """Run the QSCI algorithm.
        
        Args:
            input_states: Input quantum states for sampling
            total_shots: Total number of measurement shots
            **kwargs: Additional algorithm-specific parameters
            
        Returns:
            QSCIResult containing eigenvalues, eigenstates, and metadata
        """
        pass
    
    def _pick_out_states(
        self,
        qubit_count: int,
        meas_counts: MeasurementCounts,
        num_states_pick_out: Optional[int] = None
    ) -> Sequence[ComputationalBasisState]:
        """Pick the most frequent states from measurement counts."""
        num_states = num_states_pick_out or self.num_states_pick_out
        keys_sorted = sorted(meas_counts.keys(), key=lambda x: meas_counts[x], reverse=True)
        
        if num_states:
            keys_sorted = keys_sorted[:num_states]
        
        states = [ComputationalBasisState(qubit_count, bits=key) for key in keys_sorted]
        return states
    
    def _generate_truncated_hamiltonian(
        self,
        states: Sequence[ComputationalBasisState]
    ) -> npt.NDArray[np.complex128]:
        """Generate truncated Hamiltonian matrix in the selected subspace."""
        from quri_parts.core.operator.representation import (
            transition_amp_comp_basis,
            transition_amp_representation,
        )
        from scipy.sparse import coo_matrix
        
        dim = len(states)
        values = []
        row_ids = []
        column_ids = []
        
        h_transition_amp_repr = transition_amp_representation(self.hamiltonian)
        
        for m in range(dim):
            for n in range(m, dim):
                mn_val = transition_amp_comp_basis(
                    h_transition_amp_repr, states[m].bits, states[n].bits
                )
                if mn_val:
                    values.append(mn_val)
                    row_ids.append(m)
                    column_ids.append(n)
                    if m != n:
                        values.append(mn_val.conjugate())
                        row_ids.append(n)
                        column_ids.append(m)
        
        truncated_hamiltonian = coo_matrix(
            (values, (row_ids, column_ids)), shape=(dim, dim)
        ).tocsc(copy=False)
        truncated_hamiltonian.eliminate_zeros()
        
        return truncated_hamiltonian.toarray()
    
    def _diagonalize_truncated_hamiltonian(
        self,
        truncated_hamiltonian: npt.NDArray[np.complex128],
        k: int
    ) -> Tuple[List[float], List[List[float]]]:
        """Diagonalize truncated Hamiltonian and return k eigenvalues and eigenvectors."""
        import scipy.linalg
        import scipy.sparse.linalg
        
        n = truncated_hamiltonian.shape[0]
        
        if n == 1:
            eigvals = [truncated_hamiltonian[0, 0].real]
            eigvecs = [[1.0]]
        elif k >= n - 1:
            # Use dense diagonalization
            eigvals_all, eigvecs_all = scipy.linalg.eigh(truncated_hamiltonian)
            eigvals = eigvals_all[:k].tolist()
            eigvecs = eigvecs_all.T[:k].tolist()
        else:
            # Use sparse diagonalization
            from scipy.sparse import csc_matrix
            sparse_H = csc_matrix(truncated_hamiltonian)
            eigvals_all, eigvecs_all = scipy.sparse.linalg.eigsh(sparse_H, k, which="SA")
            eigvals = eigvals_all.tolist()
            eigvecs = eigvecs_all.T[:k].tolist()  # Fix: Take only k eigenvectors to match dense branch
        
        return eigvals, eigvecs
    
    def _calculate_computational_basis_probabilities(
        self,
        quantum_state: CircuitQuantumState,
        backend: Optional[str] = "qulacs_vector"
    ) -> Dict[int, float]:
        """Calculate probability of measuring each computational basis state.
        
        This computes P(|x⟩) = |⟨x|ψ⟩|² for all computational basis states |x⟩.
        Essential for Figure 1 reproduction in TE-QSCI paper.
        
        Args:
            quantum_state: The quantum state to analyze
            backend: Backend for exact simulation ("qulacs_vector" or "numpy")
            
        Returns:
            Dictionary mapping computational basis state (as int) to probability
        """
        try:
            if backend == "qulacs_vector":
                # Use Qulacs vector simulator for exact state vector computation
                from quri_parts.qulacs.simulator import evaluate_state_to_vector
                
                # Get the exact state vector
                state_vector = evaluate_state_to_vector(quantum_state)
                
                # Calculate probabilities for each computational basis state
                probabilities = {}
                n_qubits = quantum_state.qubit_count
                
                # Handle QuantumStateVector object from quri-parts
                if hasattr(state_vector, 'vector'):
                    # quri-parts QuantumStateVector object
                    vector = state_vector.vector
                    for i in range(len(vector)):
                        probability = abs(vector[i])**2
                        if probability > 1e-12:  # Only store non-negligible probabilities
                            probabilities[i] = probability
                elif hasattr(state_vector, '__getitem__'):
                    # Direct indexing (numpy array)
                    for i in range(2**n_qubits):
                        probability = abs(state_vector[i])**2
                        if probability > 1e-12:
                            probabilities[i] = probability
                else:
                    # Try to convert to numpy array
                    import numpy as np
                    arr = np.asarray(state_vector)
                    if arr.ndim > 0:  # Check if array has dimensions
                        for i in range(len(arr)):
                            probability = abs(arr[i])**2
                            if probability > 1e-12:
                                probabilities[i] = probability
                        
                return probabilities
                
            else:
                raise NotImplementedError(f"Backend {backend} not implemented for probability calculation")
                
        except Exception as e:
            print(f"Warning: Could not calculate exact probabilities: {e}")
            # Fallback to sampling-based probability estimation
            return self._estimate_probabilities_from_sampling(quantum_state)
    
    def _estimate_probabilities_from_sampling(
        self,
        quantum_state: CircuitQuantumState,
        shots: int = 10000
    ) -> Dict[int, float]:
        """Estimate probabilities from sampling (fallback method)."""
        if not self.sampler:
            return {}
            
        # Sample the quantum state
        meas_counts = self.sampler([(quantum_state.circuit, shots)])[0]
        
        # Convert to probabilities
        probabilities = {}
        for bits, counts in meas_counts.items():
            probabilities[bits] = counts / shots
            
        return probabilities


class VanillaQSCI(QSCIBase):
    """Vanilla QSCI algorithm implementation."""
    
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """Run vanilla QSCI algorithm."""
        if not self.sampler:
            raise ValueError("Sampler is required for vanilla QSCI")
        
        start_time = kwargs.get("start_time")
        qubit_count = input_states[0].qubit_count
        num_eigs_calc = len(input_states)
        
        # Validate num_states_pick_out
        if self.num_states_pick_out and self.num_states_pick_out < num_eigs_calc:
            raise ValueError(
                "num_states_pick_out must be larger than or equal to the number of input_states"
            )
        
        # Sample from input states
        circuits = [state.circuit for state in input_states]
        meas_counts_list = self.sampler(
            [(circuit, total_shots // num_eigs_calc) for circuit in circuits]
        )
        
        # Merge measurement counts
        merged_meas_counts: Dict[int, Union[int, float]] = {}
        for meas_counts in meas_counts_list:
            for bits, counts in meas_counts.items():
                merged_meas_counts[bits] = merged_meas_counts.get(bits, 0) + counts
        
        # Select important states
        selected_states = self._pick_out_states(
            qubit_count, merged_meas_counts, self.num_states_pick_out
        )
        
        # Generate and diagonalize truncated Hamiltonian
        truncated_hamiltonian = self._generate_truncated_hamiltonian(selected_states)
        # Fix: Use adaptive eigenvalue computation instead of fixed num_eigs_calc=1
        k = min(len(selected_states), max(1, len(selected_states) // 4))  # 25% of subspace
        eigvals, eigvecs = self._diagonalize_truncated_hamiltonian(
            truncated_hamiltonian, k
        )
        
        # Construct eigenstates
        eigenstates = [(eigvecs[i], selected_states) for i in range(min(k, len(eigvecs)))]
        
        execution_time = None
        if start_time:
            import time
            execution_time = time.time() - start_time
        
        return QSCIResult(
            eigenvalues=eigvals,
            eigenstates=eigenstates,
            selected_states=selected_states,
            subspace_dimension=len(selected_states),
            total_shots=total_shots,
            algorithm_variant=QSCIVariant.VANILLA,
            execution_time=execution_time
        )


class TimeEvolvedQSCI(QSCIBase):
    """Time-Evolved QSCI (TE-QSCI) algorithm implementation."""
    
    def __init__(
        self,
        hamiltonian: Operator,
        sampler: Optional[ConcurrentSampler] = None,
        num_states_pick_out: Optional[int] = None,
        time_evolution_method: str = "trotter"
    ):
        """Initialize TE-QSCI algorithm.
        
        Args:
            hamiltonian: Target Hamiltonian
            sampler: Quantum sampler (None for state vector variant)
            num_states_pick_out: Number of states to select
            time_evolution_method: Method for time evolution ("trotter", "exact")
        """
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.time_evolution_method = time_evolution_method
        
        # Store hamiltonian input for lazy factory creation
        self._hamiltonian_input = None  # Will be created when needed
        self._trotter_factory = None    # Lazy instantiation
        self._exact_factory = None      # Lazy instantiation
    
    def run_single_time(
        self,
        initial_state: CircuitQuantumState,
        evolution_time: float,
        total_shots: int,
        trotter_steps: Optional[int] = None,
        **kwargs
    ) -> QSCIResult:
        """Run single-time TE-QSCI algorithm.
        
        Args:
            initial_state: Initial quantum state
            evolution_time: Time for evolution
            total_shots: Number of measurement shots
            trotter_steps: Number of Trotter steps (if using Trotter method)
            
        Returns:
            QSCIResult with TE-QSCI results
        """
        if not self.sampler:
            raise ValueError("Sampler is required for single-time TE-QSCI")
        
        start_time = kwargs.get("start_time")
        
        # Create time-evolved state
        time_evolved_state = self._create_time_evolved_state(
            initial_state, evolution_time, trotter_steps
        )
        
        # Calculate computational basis probabilities for Figure 1
        probabilities = self._calculate_computational_basis_probabilities(time_evolved_state)
        
        # Run QSCI on time-evolved state
        result = self._run_qsci_on_states(
            [time_evolved_state], 
            total_shots, 
            QSCIVariant.SINGLE_TIME_TE,
            start_time
        )
        
        # Add probability information to result
        result.computational_basis_probabilities = probabilities
        
        return result
    
    def run_time_average(
        self,
        initial_state: CircuitQuantumState,
        evolution_times: Sequence[float],
        shots_per_time: int,
        trotter_steps: Optional[int] = None,
        **kwargs
    ) -> QSCIResult:
        """Run time-average TE-QSCI algorithm.
        
        Args:
            initial_state: Initial quantum state
            evolution_times: Times for evolution
            shots_per_time: Number of shots per time point
            trotter_steps: Number of Trotter steps
            
        Returns:
            QSCIResult with time-average TE-QSCI results
        """
        if not self.sampler:
            raise ValueError("Sampler is required for time-average TE-QSCI")
        
        start_time = kwargs.get("start_time")
        
        # Create time-evolved states
        time_evolved_states = []
        for t in evolution_times:
            evolved_state = self._create_time_evolved_state(
                initial_state, t, trotter_steps
            )
            time_evolved_states.append(evolved_state)
        
        # Sample from all time-evolved states
        qubit_count = initial_state.qubit_count
        circuits = [state.circuit for state in time_evolved_states]
        meas_counts_list = self.sampler(
            [(circuit, shots_per_time) for circuit in circuits]
        )
        
        # Merge all measurement counts
        merged_meas_counts: Dict[int, Union[int, float]] = {}
        for meas_counts in meas_counts_list:
            for bits, counts in meas_counts.items():
                merged_meas_counts[bits] = merged_meas_counts.get(bits, 0) + counts
        
        # Select states and diagonalize
        selected_states = self._pick_out_states(
            qubit_count, merged_meas_counts, self.num_states_pick_out
        )
        
        truncated_hamiltonian = self._generate_truncated_hamiltonian(selected_states)
        eigvals, eigvecs = self._diagonalize_truncated_hamiltonian(
            truncated_hamiltonian, 1  # Typically only ground state for time-average
        )
        
        eigenstates = [(eigvecs[i], selected_states) for i in range(len(eigvals))]
        
        execution_time = None
        if start_time:
            import time
            execution_time = time.time() - start_time
        
        return QSCIResult(
            eigenvalues=eigvals,
            eigenstates=eigenstates,
            selected_states=selected_states,
            subspace_dimension=len(selected_states),
            total_shots=len(evolution_times) * shots_per_time,
            algorithm_variant=QSCIVariant.TIME_AVERAGE_TE,
            execution_time=execution_time
        )
    
    def run_state_vector(
        self,
        initial_state: QuantumState,
        evolution_time: float,
        num_eigenstates: int = 1,
        **kwargs
    ) -> QSCIResult:
        """Run TE-QSCI with direct state vector calculation.
        
        Args:
            initial_state: Initial quantum state
            evolution_time: Time for evolution
            num_eigenstates: Number of eigenstates to compute
            
        Returns:
            QSCIResult with state vector TE-QSCI results
        """
        start_time = kwargs.get("start_time")
        
        # Create time-evolved state using exact simulation
        time_evolved_state = self._create_time_evolved_state(
            initial_state, evolution_time, kwargs.get("trotter_steps")
        )
        
        # Calculate exact computational basis probabilities
        probabilities = self._calculate_computational_basis_probabilities(time_evolved_state)
        
        # For state vector method, select the most probable states for QSCI subspace
        sorted_probs = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
        num_states_to_select = self.num_states_pick_out or min(100, len(sorted_probs))
        
        # Create computational basis states from most probable outcomes
        qubit_count = initial_state.qubit_count
        selected_states = []
        for bits_int, prob in sorted_probs[:num_states_to_select]:
            state = ComputationalBasisState(qubit_count, bits=bits_int)
            selected_states.append(state)
        
        # Generate and diagonalize truncated Hamiltonian
        truncated_hamiltonian = self._generate_truncated_hamiltonian(selected_states)
        eigvals, eigvecs = self._diagonalize_truncated_hamiltonian(
            truncated_hamiltonian, num_eigenstates
        )
        
        # Construct eigenstates
        eigenstates = [(eigvecs[i], selected_states) for i in range(num_eigenstates)]
        
        execution_time = None
        if start_time:
            import time
            execution_time = time.time() - start_time
        
        return QSCIResult(
            eigenvalues=eigvals,
            eigenstates=eigenstates,
            selected_states=selected_states,
            subspace_dimension=len(selected_states),
            total_shots=0,  # No sampling in state vector method
            algorithm_variant=QSCIVariant.STATE_VECTOR,
            execution_time=execution_time,
            computational_basis_probabilities=probabilities
        )
    
    def _create_time_evolved_state(
        self,
        initial_state: CircuitQuantumState,
        evolution_time: float,
        trotter_steps: Optional[int] = None
    ) -> CircuitQuantumState:
        """Create time-evolved quantum state using quri-algo time evolution.
        
        Supports both Trotter and exact unitary evolution methods.
        """
        from quri_parts.core.state import GeneralCircuitQuantumState
        from quri_parts.circuit import QuantumCircuit
        
        # Create QubitHamiltonian if not already created
        if self._hamiltonian_input is None:
            from quri_algo.problem import QubitHamiltonian
            self._hamiltonian_input = QubitHamiltonian(initial_state.qubit_count, self.hamiltonian)
        
        try:
            if self.time_evolution_method == "exact":
                # Use exact unitary time evolution
                if self._exact_factory is None:
                    from quri_algo.circuit.time_evolution.exact_unitary import ExactUnitaryTimeEvolutionCircuitFactory
                    self._exact_factory = ExactUnitaryTimeEvolutionCircuitFactory(self._hamiltonian_input)
                factory = self._exact_factory
                
            elif self.time_evolution_method == "trotter":
                # Use Trotter approximation
                if trotter_steps is None:
                    trotter_steps = max(1, int(evolution_time * 10))  # Default heuristic
                
                if self._trotter_factory is None:
                    from quri_algo.circuit.time_evolution.trotter_time_evo import TrotterTimeEvolutionCircuitFactory
                    self._trotter_factory = TrotterTimeEvolutionCircuitFactory(
                        self._hamiltonian_input, n_trotter=trotter_steps
                    )
                factory = self._trotter_factory
                
            else:
                raise ValueError(f"Unknown time evolution method: {self.time_evolution_method}")
            
            # Generate time evolution circuit using the selected factory
            time_evolution_circuit = factory(evolution_time)
            
            # Combine initial state circuit with time evolution
            evolved_circuit = QuantumCircuit(initial_state.qubit_count)
            
            # Add initial state preparation
            for gate in initial_state.circuit.gates:
                evolved_circuit.add_gate(gate)
            
            # Add time evolution
            for gate in time_evolution_circuit.gates:
                evolved_circuit.add_gate(gate)
                
        except ImportError as e:
            # Fallback to simplified Trotter implementation for backwards compatibility
            print(f"Warning: quri-algo not available ({e}), attempting fallback implementation")
            if trotter_steps is None:
                trotter_steps = max(1, int(evolution_time * 10))
            try:
                evolved_circuit = self._create_fallback_trotter_circuit(
                    initial_state.circuit, evolution_time, trotter_steps
                )
            except NotImplementedError as nie:
                raise NotImplementedError(
                    f"Time evolution failed: {nie}. "
                    f"To use time evolution, please install quri-algo or use 'exact' method "
                    f"for small systems."
                ) from nie
        
        return GeneralCircuitQuantumState(initial_state.qubit_count, evolved_circuit)
    
    def _create_fallback_trotter_circuit(
        self,
        initial_circuit: NonParametricQuantumCircuit,
        evolution_time: float,
        trotter_steps: int
    ) -> NonParametricQuantumCircuit:
        """Fallback Trotter implementation when quri-algo is not available."""
        from quri_parts.circuit import QuantumCircuit
        
        # Create new circuit starting with the initial circuit
        evolved_circuit = QuantumCircuit(initial_circuit.qubit_count)
        
        # Add all gates from initial circuit
        for gate in initial_circuit.gates:
            evolved_circuit.add_gate(gate)
        
        # Time step for each Trotter slice
        dt = evolution_time / trotter_steps
        
        # Apply simplified Trotter steps
        for step in range(trotter_steps):
            # For each Trotter step, apply evolution for each Pauli term
            self._add_fallback_trotter_step(evolved_circuit, dt)
        
        return evolved_circuit
    
    def _add_fallback_trotter_step(self, circuit: QuantumCircuit, dt: float) -> None:
        """Add a single Trotter step using simplified implementation."""
        raise NotImplementedError(
            "Fallback Trotter implementation is not available. "
            "This would require a full Pauli operator decomposition and proper "
            "time evolution implementation. Please ensure quri-algo is properly "
            "installed and use the 'exact' time evolution method for testing, "
            "or implement a proper Pauli decomposition for your Hamiltonian."
        )
    
    
    def _run_qsci_on_states(
        self,
        states: Sequence[CircuitQuantumState],
        total_shots: int,
        variant: QSCIVariant,
        start_time: Optional[float] = None
    ) -> QSCIResult:
        """Helper method to run QSCI on given states."""
        qubit_count = states[0].qubit_count
        num_eigs_calc = len(states)
        
        # Sample from states
        circuits = [state.circuit for state in states]
        meas_counts_list = self.sampler(
            [(circuit, total_shots // num_eigs_calc) for circuit in circuits]
        )
        
        # Merge measurement counts
        merged_meas_counts: Dict[int, Union[int, float]] = {}
        for meas_counts in meas_counts_list:
            for bits, counts in meas_counts.items():
                merged_meas_counts[bits] = merged_meas_counts.get(bits, 0) + counts
        
        # Select states and diagonalize
        selected_states = self._pick_out_states(
            qubit_count, merged_meas_counts, self.num_states_pick_out
        )
        
        truncated_hamiltonian = self._generate_truncated_hamiltonian(selected_states)
        eigvals, eigvecs = self._diagonalize_truncated_hamiltonian(
            truncated_hamiltonian, num_eigs_calc
        )
        
        eigenstates = [(eigvecs[i], selected_states) for i in range(num_eigs_calc)]
        
        execution_time = None
        if start_time:
            import time
            execution_time = time.time() - start_time
        
        return QSCIResult(
            eigenvalues=eigvals,
            eigenstates=eigenstates,
            selected_states=selected_states,
            subspace_dimension=len(selected_states),
            total_shots=total_shots,
            algorithm_variant=variant,
            execution_time=execution_time
        )

    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """General run method for TE-QSCI (defaults to single-time)."""
        if len(input_states) != 1:
            raise ValueError("TE-QSCI expects exactly one initial state")
        
        evolution_time = kwargs.get("evolution_time", 1.0)
        trotter_steps = kwargs.get("trotter_steps")
        
        return self.run_single_time(
            input_states[0], evolution_time, total_shots, trotter_steps, **kwargs
        )


# Wrapper classes for testing compatibility
class SingleTimeTE_QSCI(TimeEvolvedQSCI):
    """Single-time TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_time, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(self, input_states, total_shots, **kwargs):
        if len(input_states) != 1:
            raise ValueError("SingleTimeTE_QSCI expects exactly one initial state")
        return self.run_single_time(
            input_states[0], self.evolution_time, total_shots, **kwargs
        )


class TimeAverageTE_QSCI(TimeEvolvedQSCI):
    """Time-average TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_times, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_times = evolution_times
    
    def run(self, input_states, total_shots, **kwargs):
        if len(input_states) != 1:
            raise ValueError("TimeAverageTE_QSCI expects exactly one initial state")
        shots_per_time = total_shots // len(self.evolution_times)
        return self.run_time_average(
            input_states[0], self.evolution_times, shots_per_time, **kwargs
        )


class FfsimStateVectorQSCI(QSCIBase):
    """QSCI algorithm for direct ffsim state vector input.
    
    This class provides direct integration with ffsim's LUCJ/UCJ ansatz state vectors,
    bypassing circuit sampling and using exact computational basis probabilities
    for more accurate QSCI analysis.
    """
    
    def __init__(
        self,
        hamiltonian: Operator,
        num_states_pick_out: Optional[int] = None
    ):
        """Initialize FfsimStateVectorQSCI.
        
        Args:
            hamiltonian: Target Hamiltonian to diagonalize (quri-parts Operator)
            num_states_pick_out: Number of states to select for subspace
        """
        # No sampler needed for state vector method
        super().__init__(hamiltonian, sampler=None, num_states_pick_out=num_states_pick_out)
    
    def run(
        self,
        input_states: Sequence[CircuitQuantumState],
        total_shots: int,
        **kwargs
    ) -> QSCIResult:
        """Run QSCI algorithm - maintained for interface compatibility.
        
        Note: This method maintains the standard QSCI interface but is not the primary
        entry point for ffsim integration. Use run_ffsim_state_vector() instead.
        """
        raise NotImplementedError(
            "FfsimStateVectorQSCI is designed for direct state vector input. "
            "Use run_ffsim_state_vector() method instead."
        )
    
    def run_ffsim_state_vector(
        self,
        ffsim_state_vector: np.ndarray,
        n_qubits: int,
        num_eigenstates: int = 1,
        nelec: tuple[int, int] = None,
        **kwargs
    ) -> QSCIResult:
        """Run QSCI with direct ffsim state vector input.
        
        Args:
            ffsim_state_vector: State vector from ffsim LUCJ/UCJ ansatz (numpy array)
            n_qubits: Number of qubits in the system (typically 2 * norb)
            num_eigenstates: Number of eigenstates to compute
            **kwargs: Additional parameters
            
        Returns:
            QSCIResult with exact QSCI analysis results
        """
        import time
        start_time = kwargs.get("start_time", time.time())
        
        print(f"Running FfsimStateVectorQSCI...")
        print(f"  Input state vector dimension: {len(ffsim_state_vector)}")
        print(f"  Target qubit space dimension: 2^{n_qubits} = {2**n_qubits}")
        
        # Convert ffsim state vector to computational basis probabilities
        probabilities = self._ffsim_to_computational_probabilities(
            ffsim_state_vector, n_qubits, nelec
        )
        
        print(f"  Extracted {len(probabilities)} non-zero computational basis states")
        
        # Select the most probable states for QSCI subspace
        sorted_probs = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
        num_states_to_select = self.num_states_pick_out or min(100, len(sorted_probs))
        
        print(f"  Selecting top {num_states_to_select} states for QSCI subspace")
        
        # Create computational basis states from most probable outcomes  
        selected_states = []
        for bits_int, prob in sorted_probs[:num_states_to_select]:
            state = ComputationalBasisState(n_qubits, bits=bits_int)
            selected_states.append(state)
            
        # Show selected states for debugging
        print(f"  Top selected states:")
        for i, (bits_int, prob) in enumerate(sorted_probs[:min(5, num_states_to_select)]):
            bits_str = format(bits_int, f'0{n_qubits}b')
            print(f"    {i+1}: |{bits_str}⟩ = {bits_int}, prob = {prob:.6f}")
        
        # Generate and diagonalize truncated Hamiltonian
        print(f"  Generating truncated Hamiltonian in {len(selected_states)}-dimensional subspace...")
        truncated_hamiltonian = self._generate_truncated_hamiltonian(selected_states)
        
        print(f"  Diagonalizing for {num_eigenstates} eigenstate(s)...")
        eigvals, eigvecs = self._diagonalize_truncated_hamiltonian(
            truncated_hamiltonian, num_eigenstates
        )
        
        # Construct eigenstates
        eigenstates = [(eigvecs[i], selected_states) for i in range(num_eigenstates)]
        
        execution_time = time.time() - start_time
        
        print(f"  QSCI ground state energy: {eigvals[0]:.6f} Ha")
        print(f"  Execution time: {execution_time:.3f} seconds")
        
        return QSCIResult(
            eigenvalues=eigvals,
            eigenstates=eigenstates,
            selected_states=selected_states,
            subspace_dimension=len(selected_states),
            total_shots=0,  # No sampling in state vector method
            algorithm_variant=QSCIVariant.STATE_VECTOR,
            execution_time=execution_time,
            computational_basis_probabilities=probabilities
        )
    
    def _ffsim_to_computational_probabilities(
        self,
        ffsim_state_vector: np.ndarray,
        n_qubits: int,
        nelec: tuple[int, int] = None
    ) -> Dict[int, float]:
        """Convert ffsim state vector to computational basis probabilities.
        
        Uses the corrected fermionic-to-computational mapping with proper
        Jordan-Wigner bit ordering.
        """
        # Import from ffsim module
        try:
            from ffsim_integration.state_conversion import _map_fermionic_amplitudes_directly
            
            # Convert fermionic state vector to computational basis amplitudes
            computational_amplitudes = _map_fermionic_amplitudes_directly(
                ffsim_state_vector, n_qubits, nelec
            )
            
            # Calculate probabilities |ψ_i|²
            probabilities = {}
            for i, amplitude in enumerate(computational_amplitudes):
                prob = abs(amplitude) ** 2
                if prob > 1e-12:  # Only store significant probabilities
                    probabilities[i] = prob
                    
            # Verify normalization
            total_prob = sum(probabilities.values())
            print(f"  Total probability: {total_prob:.6f}")
            
            if abs(total_prob - 1.0) > 1e-6:
                print(f"  Warning: State not properly normalized, renormalizing...")
                for i in probabilities:
                    probabilities[i] /= total_prob
                    
            return probabilities
            
        except ImportError as e:
            raise ImportError(
                f"Cannot import ffsim state conversion utilities: {e}. "
                f"Install ffsim integration with: pip install quri-qsci[ffsim]"
            )


class StateVectorTE_QSCI(TimeEvolvedQSCI):
    """State vector TE-QSCI wrapper for testing compatibility."""
    
    def __init__(self, hamiltonian, sampler, evolution_time, num_states_pick_out=None):
        super().__init__(hamiltonian, sampler, num_states_pick_out)
        self.evolution_time = evolution_time
    
    def run(self, input_states, total_shots, **kwargs):
        if len(input_states) != 1:
            raise ValueError("StateVectorTE_QSCI expects exactly one initial state")
        
        initial_state = input_states[0]
        
        # Validate that the state has the required circuit attribute
        if not hasattr(initial_state, 'circuit'):
            raise TypeError(
                f"StateVectorTE_QSCI requires a GeneralCircuitQuantumState with a 'circuit' attribute, "
                f"but received {type(initial_state)}. Please ensure you're using "
                f"create_circuit_state_from_amplitudes() instead of create_custom_state_vector()."
            )
        
        # For state vector method, we don't need shots but keeping interface consistent
        return self.run_state_vector(
            initial_state, self.evolution_time, **kwargs
        )