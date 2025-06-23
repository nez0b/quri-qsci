"""
Dual exact/sampling probability calculator for TE-QSCI excitation analysis.

This module implements both exact and sampling-based probability calculation
methods for analyzing excitation probabilities in time-evolved quantum states.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from quri_parts.core.state import GeneralCircuitQuantumState
from quri_parts.qulacs.simulator import evaluate_state_to_vector
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
from itertools import combinations
import time

class ProbabilityCalculator:
    """Dual exact/sampling probability calculator for quantum states."""
    
    def __init__(self, method: str = "auto", verbose: bool = True):
        """Initialize probability calculator.
        
        Args:
            method: "exact", "sampling", or "auto" for automatic selection
            verbose: Whether to print progress information
        """
        self.method = method
        self.verbose = verbose
        self.sampler = None
        
        # Performance thresholds for auto-selection
        self.exact_qubit_threshold = 14  # Use exact for ≤14 qubits
        self.exact_state_threshold = 16384  # Use exact for ≤16K states
        
    def select_states_at_t1(
        self, 
        evolved_state_t1: GeneralCircuitQuantumState,
        n_electrons: int,
        n_qubits: int,
        R: int = 850
    ) -> Tuple[np.ndarray, Dict[int, List[int]]]:
        """Select top R states at t=1 and classify by excitation order.
        
        Always uses exact calculation for state selection regardless of method.
        
        Args:
            evolved_state_t1: Quantum state evolved to t=1
            n_electrons: Number of electrons in the system
            n_qubits: Number of qubits in the system
            R: Number of top states to select
            
        Returns:
            Tuple of (selected_state_indices, state_orders_dict)
        """
        if self.verbose:
            print(f"Selecting top R={R} states at t=1 using exact calculation...")
        
        # Get exact state vector
        state_vector = evaluate_state_to_vector(evolved_state_t1)
        if hasattr(state_vector, 'vector'):
            probabilities = np.abs(state_vector.vector)**2
        else:
            probabilities = np.abs(state_vector)**2
        
        # Filter to only valid electron configurations
        valid_states = []
        valid_probabilities = []
        
        for i in range(len(probabilities)):
            if bin(i).count('1') == n_electrons:
                valid_states.append(i)
                valid_probabilities.append(probabilities[i])
        
        valid_states = np.array(valid_states)
        valid_probabilities = np.array(valid_probabilities)
        
        if self.verbose:
            print(f"✓ Found {len(valid_states)} valid electron configurations")
        
        # Select top R states from valid configurations
        if len(valid_states) <= R:
            selected_indices = valid_states
            if self.verbose:
                print(f"✓ Using all {len(valid_states)} valid states (fewer than R={R})")
        else:
            top_indices = np.argsort(valid_probabilities)[-R:]
            selected_indices = valid_states[top_indices]
            if self.verbose:
                selected_probs = valid_probabilities[top_indices]
                print(f"✓ Selected top {R} states")
                print(f"✓ Probability range: {selected_probs.min():.6f} to {selected_probs.max():.6f}")
        
        # Classify selected states by excitation order
        hf_state_idx = (1 << n_electrons) - 1  # |11111...000...⟩
        selected_state_orders = {}
        
        for state_idx in selected_indices:
            order = self._classify_excitation_order(state_idx, hf_state_idx, n_electrons)
            if order >= 0:
                if order not in selected_state_orders:
                    selected_state_orders[order] = []
                selected_state_orders[order].append(state_idx)
        
        if self.verbose:
            print("✓ Selected state distribution by excitation order:")
            for order in sorted(selected_state_orders.keys()):
                count = len(selected_state_orders[order])
                print(f"   Order {order}: {count} states")
        
        return selected_indices, selected_state_orders
    
    def calculate_probabilities(
        self,
        evolved_state: GeneralCircuitQuantumState,
        selected_states: np.ndarray,
        n_qubits: int,
        method_override: Optional[str] = None
    ) -> np.ndarray:
        """Calculate probabilities for selected states using chosen method.
        
        Args:
            evolved_state: Time-evolved quantum state
            selected_states: Array of selected state indices
            n_qubits: Number of qubits
            method_override: Override the default method selection
            
        Returns:
            Array of probabilities for selected states
        """
        # Determine method to use
        if method_override:
            use_method = method_override
        elif self.method == "auto":
            use_method = self._auto_select_method(n_qubits, len(selected_states))
        else:
            use_method = self.method
        
        if use_method == "exact":
            return self._exact_probabilities(evolved_state, selected_states)
        elif use_method == "sampling":
            return self._sampling_probabilities(evolved_state, selected_states, n_qubits)
        else:
            raise ValueError(f"Unknown method: {use_method}")
    
    def _auto_select_method(self, n_qubits: int, n_selected_states: int) -> str:
        """Automatically select calculation method based on system size."""
        if n_qubits <= self.exact_qubit_threshold and 2**n_qubits <= self.exact_state_threshold:
            if self.verbose:
                print(f"Auto-selected: exact method ({n_qubits} qubits, {2**n_qubits} states)")
            return "exact"
        else:
            if self.verbose:
                print(f"Auto-selected: sampling method ({n_qubits} qubits, {2**n_qubits} states)")
            return "sampling"
    
    def _exact_probabilities(
        self,
        evolved_state: GeneralCircuitQuantumState,
        selected_states: np.ndarray
    ) -> np.ndarray:
        """Calculate exact probabilities from state vector."""
        state_vector = evaluate_state_to_vector(evolved_state)
        if hasattr(state_vector, 'vector'):
            full_probabilities = np.abs(state_vector.vector)**2
        else:
            full_probabilities = np.abs(state_vector)**2
        
        return full_probabilities[selected_states]
    
    def _sampling_probabilities(
        self,
        evolved_state: GeneralCircuitQuantumState,
        selected_states: np.ndarray,
        n_qubits: int,
        base_samples: int = 100000
    ) -> np.ndarray:
        """Calculate probabilities using direct state vector sampling.
        
        For time-evolved states, we need to sample from the state vector directly
        rather than executing circuits.
        
        Args:
            evolved_state: Time-evolved quantum state
            selected_states: Array of selected state indices
            n_qubits: Number of qubits
            base_samples: Base number of samples to use
            
        Returns:
            Array of estimated probabilities for selected states
        """
        # Adaptive sampling: more samples for systems with more rare events
        n_samples = min(base_samples * max(1, len(selected_states) // 100), 10**7)
        
        if self.verbose:
            print(f"Sampling with {n_samples:,} shots...")
        
        # Get exact state vector for sampling
        from quri_parts.qulacs.simulator import evaluate_state_to_vector
        state_vector = evaluate_state_to_vector(evolved_state)
        
        if hasattr(state_vector, 'vector'):
            full_probabilities = np.abs(state_vector.vector)**2
        else:
            full_probabilities = np.abs(state_vector)**2
        
        # Normalize probabilities (should already be normalized, but just in case)
        full_probabilities = full_probabilities / np.sum(full_probabilities)
        
        # Sample from the probability distribution
        # Generate random samples according to the state vector probabilities
        state_indices = np.arange(len(full_probabilities))
        sampled_states = np.random.choice(
            state_indices, 
            size=n_samples, 
            p=full_probabilities
        )
        
        # Count occurrences of selected states
        selected_state_set = set(selected_states)
        state_counts = {state: 0 for state in selected_states}
        
        for sampled_state in sampled_states:
            if sampled_state in selected_state_set:
                state_counts[sampled_state] += 1
        
        # Convert counts to probabilities
        probabilities = np.array([state_counts[state] / n_samples for state in selected_states])
        
        if self.verbose:
            non_zero_count = np.sum(probabilities > 0)
            print(f"✓ Measured {non_zero_count}/{len(selected_states)} selected states")
            
            # Compare with exact probabilities for validation
            exact_probabilities = full_probabilities[selected_states]
            max_error = np.max(np.abs(probabilities - exact_probabilities))
            mean_error = np.mean(np.abs(probabilities - exact_probabilities))
            print(f"✓ Sampling accuracy: max error = {max_error:.6f}, mean error = {mean_error:.6f}")
        
        return probabilities
    
    def _classify_excitation_order(self, state_idx: int, hf_state_idx: int, n_electrons: int) -> int:
        """Classify excitation order of a computational basis state."""
        # Check for correct number of electrons
        if bin(state_idx).count('1') != n_electrons:
            return -1
        
        # Find differing bits
        diff = state_idx ^ hf_state_idx
        num_flips = bin(diff).count('1')
        
        # Each excitation causes two flips (particle-hole pair)
        excitation_order = num_flips // 2
        
        return excitation_order
    
    def calculate_grouped_probabilities(
        self,
        probabilities: np.ndarray,
        selected_state_orders: Dict[int, List[int]],
        selected_states: np.ndarray
    ) -> Tuple[float, float, float]:
        """Calculate grouped probabilities for Figure 1 reproduction.
        
        Args:
            probabilities: Probabilities for selected states
            selected_state_orders: Dictionary mapping excitation order to state indices
            selected_states: Array of selected state indices
            
        Returns:
            Tuple of (prob_1_2, prob_3_4, prob_5_6) for grouped excitation orders
        """
        # Create mapping from state index to probability
        state_to_prob = {selected_states[i]: probabilities[i] for i in range(len(selected_states))}
        
        # Group 1: One/Two-electron excitations (orders 1,2)
        prob_1_2 = 0.0
        for order in [1, 2]:
            if order in selected_state_orders:
                prob_1_2 += sum(state_to_prob.get(state_idx, 0.0) for state_idx in selected_state_orders[order])
        
        # Group 2: Three/Four-electron excitations (orders 3,4)
        prob_3_4 = 0.0
        for order in [3, 4]:
            if order in selected_state_orders:
                prob_3_4 += sum(state_to_prob.get(state_idx, 0.0) for state_idx in selected_state_orders[order])
        
        # Group 3: Five/Six-electron excitations (orders 5,6)
        prob_5_6 = 0.0
        for order in [5, 6]:
            if order in selected_state_orders:
                prob_5_6 += sum(state_to_prob.get(state_idx, 0.0) for state_idx in selected_state_orders[order])
        
        return prob_1_2, prob_3_4, prob_5_6

class H6FigureOneAnalyzer:
    """Specialized analyzer for reproducing Figure 1 with H6 molecule."""
    
    def __init__(self, hamiltonian, hartree_fock_state, active_space, method: str = "auto"):
        """Initialize H6 Figure 1 analyzer.
        
        Args:
            hamiltonian: H6 Hamiltonian operator
            hartree_fock_state: Hartree-Fock initial state
            active_space: Active space information
            method: Probability calculation method
        """
        self.hamiltonian = hamiltonian
        self.hartree_fock_state = hartree_fock_state
        self.active_space = active_space
        self.calculator = ProbabilityCalculator(method=method)
        
        # H6 system parameters
        self.n_electrons = active_space.n_active_ele if active_space else 6
        self.n_qubits = 2 * active_space.n_active_orb if active_space else 12
        
    def run_figure_one_analysis(self, R: int = 850, time_points: Optional[np.ndarray] = None) -> Dict:
        """Run complete Figure 1 analysis for H6.
        
        Args:
            R: Number of states to select (paper uses R=850)
            time_points: Custom time points (default: focus on small-t regime)
            
        Returns:
            Dictionary containing all analysis results
        """
        if time_points is None:
            # Focus on small-t regime for scaling analysis, then broader range
            small_times = np.logspace(-2, 0, 20)  # 0.01 to 1.0
            large_times = np.linspace(1.2, 3.0, 10)  # 1.2 to 3.0
            time_points = np.concatenate([small_times, large_times])
            time_points = np.sort(time_points)
        
        print(f"=== H6 FIGURE 1 ANALYSIS ===")
        print(f"System: {self.n_electrons} electrons, {self.n_qubits} qubits")
        print(f"States to select: R = {R}")
        print(f"Time points: {len(time_points)} (range {time_points.min():.3f} to {time_points.max():.3f})")
        print(f"Method: {self.calculator.method}")
        
        # Step 1: State selection at t=1
        print(f"\n=== STEP 1: STATE SELECTION AT t=1 ===")
        start_time = time.time()
        
        evolved_state_t1 = self._evolve_exact(1.0)
        selected_states, selected_state_orders = self.calculator.select_states_at_t1(
            evolved_state_t1, self.n_electrons, self.n_qubits, R
        )
        
        selection_time = time.time() - start_time
        print(f"✓ State selection completed in {selection_time:.2f} seconds")
        
        # Step 2: Time evolution analysis
        print(f"\n=== STEP 2: TIME EVOLUTION ANALYSIS ===")
        start_time = time.time()
        
        results_1_2 = []
        results_3_4 = []
        results_5_6 = []
        
        for i, t in enumerate(time_points):
            if i % max(1, len(time_points) // 10) == 0:
                print(f"  Progress: {i+1}/{len(time_points)} (t = {t:.4f})")
            
            # Evolve state to time t
            evolved_state = self._evolve_exact(t)
            
            # Calculate probabilities for selected states
            probabilities = self.calculator.calculate_probabilities(
                evolved_state, selected_states, self.n_qubits
            )
            
            # Group by excitation orders
            prob_1_2, prob_3_4, prob_5_6 = self.calculator.calculate_grouped_probabilities(
                probabilities, selected_state_orders, selected_states
            )
            
            results_1_2.append(prob_1_2)
            results_3_4.append(prob_3_4)
            results_5_6.append(prob_5_6)
        
        evolution_time = time.time() - start_time
        print(f"✓ Time evolution analysis completed in {evolution_time:.2f} seconds")
        
        # Package results
        results = {
            'times': time_points,
            'prob_1_2': np.array(results_1_2),
            'prob_3_4': np.array(results_3_4),
            'prob_5_6': np.array(results_5_6),
            'selected_states': selected_states,
            'selected_state_orders': selected_state_orders,
            'system_info': {
                'n_electrons': self.n_electrons,
                'n_qubits': self.n_qubits,
                'R': R,
                'method': self.calculator.method
            },
            'timing': {
                'selection_time': selection_time,
                'evolution_time': evolution_time,
                'total_time': selection_time + evolution_time
            }
        }
        
        print(f"\n=== ANALYSIS COMPLETE ===")
        print(f"Total time: {results['timing']['total_time']:.2f} seconds")
        
        return results
    
    def _evolve_exact(self, evolution_time: float) -> GeneralCircuitQuantumState:
        """Evolve HF state using exact time evolution."""
        from src import TimeEvolvedQSCI
        
        te_qsci = TimeEvolvedQSCI(
            hamiltonian=self.hamiltonian,
            sampler=None,
            time_evolution_method="exact"
        )
        
        evolved_state = te_qsci._create_time_evolved_state(
            self.hartree_fock_state,
            evolution_time
        )
        
        return evolved_state