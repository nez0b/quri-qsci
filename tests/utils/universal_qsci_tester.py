"""
Universal QSCI Tester

This module provides a comprehensive testing framework for all QSCI variants,
enabling systematic precision validation and comparative analysis across
different algorithms and initial states.
"""

import numpy as np
import time
from typing import Dict, List, Optional, Any, Union, Tuple
import sys
import os

# Add parent directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from .exact_ground_state_utils import (
    extract_exact_ground_state,
    create_exact_circuit_state_from_vector,
    validate_exact_state_preparation
)

try:
    from src.qsci_algorithms import (
        VanillaQSCI, SingleTimeTE_QSCI, TimeAverageTE_QSCI, StateVectorTE_QSCI
    )
    from quri_parts.core.operator import Operator
    from quri_parts.core.state import GeneralCircuitQuantumState
    from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
    QSCI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: QSCI implementation not available: {e}")
    QSCI_AVAILABLE = False

from scipy.sparse import csc_matrix


class UniversalQSCITester:
    """
    Universal tester for all QSCI algorithm variants.
    
    This class provides a unified interface for testing and comparing
    VanillaQSCI, SingleTimeTE_QSCI, TimeAverageTE_QSCI, and StateVectorTE_QSCI
    with various initial states and precision requirements.
    """
    
    def __init__(self, 
                 default_shots: int = 2000,
                 default_num_states: int = 4,
                 precision_tolerance: float = 1e-6,
                 verbose: bool = True):
        """
        Initialize the Universal QSCI Tester.
        
        Args:
            default_shots: Default number of shots for sampling
            default_num_states: Default number of states to pick out
            precision_tolerance: Expected precision tolerance for validation
            verbose: Whether to print detailed progress information
        """
        self.default_shots = default_shots
        self.default_num_states = default_num_states
        self.precision_tolerance = precision_tolerance
        self.verbose = verbose
        
        if not QSCI_AVAILABLE:
            raise ImportError("QSCI implementation not available")
    
    def test_all_variants_with_exact_ground_state(self,
                                                  hamiltonian_matrix: csc_matrix,
                                                  quri_hamiltonian: Operator,
                                                  n_qubits: int,
                                                  test_config: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Test all QSCI variants using exact ground state as initial state.
        
        This is the core validation method that proves QSCI algorithm correctness
        by starting with the exact ground state and expecting machine precision.
        
        Args:
            hamiltonian_matrix: Sparse Hamiltonian matrix
            quri_hamiltonian: QURI-Parts Operator format Hamiltonian
            n_qubits: Number of qubits in the system
            test_config: Optional configuration overrides
            
        Returns:
            Dictionary with results for each QSCI variant
        """
        if self.verbose:
            print(f"=== UNIVERSAL QSCI TESTING WITH EXACT GROUND STATE ===")
            print(f"System: {n_qubits} qubits, {2**n_qubits} dimensional Hilbert space")
        
        # Extract exact ground state
        start_time = time.time()
        exact_ground_vector, exact_ground_energy = extract_exact_ground_state(hamiltonian_matrix)
        
        if self.verbose:
            print(f"Exact ground state energy: {exact_ground_energy:.12f}")
            print(f"Ground state extraction time: {time.time() - start_time:.3f}s")
        
        # Create exact circuit state
        exact_state = create_exact_circuit_state_from_vector(exact_ground_vector, n_qubits)
        
        # Validate state preparation
        validation = validate_exact_state_preparation(
            exact_ground_vector, exact_state, tolerance=1e-10
        )
        
        if self.verbose:
            print(f"State preparation fidelity: {validation['fidelity']:.8f}")
        
        # Configuration
        config = test_config or {}
        shots = config.get('shots', self.default_shots)
        num_states = config.get('num_states', self.default_num_states)
        
        # Test results storage
        results = {
            'exact_ground_energy': exact_ground_energy,
            'exact_ground_vector': exact_ground_vector,
            'state_preparation_validation': validation,
            'variants': {},
            'summary': {}
        }
        
        # Test each QSCI variant
        variants_to_test = [
            ('VanillaQSCI', self._test_vanilla_qsci),
            ('SingleTimeTE_QSCI', self._test_single_time_te_qsci),
            ('StateVectorTE_QSCI', self._test_state_vector_te_qsci),
            ('TimeAverageTE_QSCI', self._test_time_average_te_qsci)
        ]
        
        successful_variants = []
        failed_variants = []
        
        for variant_name, test_method in variants_to_test:
            if self.verbose:
                print(f"\n--- Testing {variant_name} ---")
            
            try:
                variant_result = test_method(
                    quri_hamiltonian, exact_state, exact_ground_energy, config
                )
                results['variants'][variant_name] = variant_result
                
                if variant_result['success']:
                    successful_variants.append(variant_name)
                    if self.verbose:
                        energy = variant_result['energy']
                        error = variant_result['error']
                        print(f"✅ {variant_name}: {energy:.12f} (error: {error:.2e})")
                else:
                    failed_variants.append(variant_name)
                    if self.verbose:
                        print(f"❌ {variant_name}: FAILED - {variant_result.get('error_msg', 'Unknown error')}")
                        
            except Exception as e:
                failed_variants.append(variant_name)
                results['variants'][variant_name] = {
                    'success': False,
                    'error_msg': str(e),
                    'execution_time': 0.0
                }
                if self.verbose:
                    print(f"❌ {variant_name}: EXCEPTION - {e}")
        
        # Summary analysis
        results['summary'] = {
            'successful_variants': successful_variants,
            'failed_variants': failed_variants,
            'success_rate': len(successful_variants) / len(variants_to_test),
            'best_variant': self._find_best_variant(results['variants']),
            'precision_achieved': self._check_precision_achieved(results['variants']),
            'total_test_time': time.time() - start_time
        }
        
        if self.verbose:
            self._print_summary(results['summary'])
        
        return results
    
    def test_variant_with_approximate_states(self,
                                           hamiltonian_matrix: csc_matrix,
                                           quri_hamiltonian: Operator,
                                           n_qubits: int,
                                           approximate_states: List[GeneralCircuitQuantumState],
                                           variant_name: str = 'VanillaQSCI',
                                           test_config: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Test a specific QSCI variant with approximate initial states.
        
        This method tests how QSCI performs with realistic (approximate) initial states,
        providing comparison against the exact ground state baseline.
        
        Args:
            hamiltonian_matrix: Sparse Hamiltonian matrix
            quri_hamiltonian: QURI-Parts Operator format Hamiltonian  
            n_qubits: Number of qubits
            approximate_states: List of approximate initial states to test
            variant_name: Which QSCI variant to test
            test_config: Optional configuration overrides
            
        Returns:
            Dictionary with results for each approximate state
        """
        if self.verbose:
            print(f"=== TESTING {variant_name} WITH APPROXIMATE STATES ===")
        
        # Get exact ground energy for comparison
        _, exact_ground_energy = extract_exact_ground_state(hamiltonian_matrix)
        
        config = test_config or {}
        results = {
            'exact_ground_energy': exact_ground_energy,
            'variant_name': variant_name,
            'approximate_state_results': {},
            'summary': {}
        }
        
        # Test each approximate state
        for i, approx_state in enumerate(approximate_states):
            state_name = f"approx_state_{i}"
            
            if self.verbose:
                print(f"\n--- Testing {state_name} ---")
            
            try:
                if variant_name == 'VanillaQSCI':
                    result = self._test_vanilla_qsci(
                        quri_hamiltonian, approx_state, exact_ground_energy, config
                    )
                elif variant_name == 'SingleTimeTE_QSCI':
                    result = self._test_single_time_te_qsci(
                        quri_hamiltonian, approx_state, exact_ground_energy, config
                    )
                elif variant_name == 'StateVectorTE_QSCI':
                    result = self._test_state_vector_te_qsci(
                        quri_hamiltonian, approx_state, exact_ground_energy, config
                    )
                else:
                    raise ValueError(f"Unsupported variant: {variant_name}")
                
                results['approximate_state_results'][state_name] = result
                
                if self.verbose and result['success']:
                    energy = result['energy']
                    error = result['error']
                    print(f"  Energy: {energy:.12f} (error: {error:.2e})")
                    
            except Exception as e:
                results['approximate_state_results'][state_name] = {
                    'success': False,
                    'error_msg': str(e)
                }
                if self.verbose:
                    print(f"  FAILED: {e}")
        
        # Analyze results
        results['summary'] = self._analyze_approximate_state_results(
            results['approximate_state_results']
        )
        
        return results
    
    def comparative_precision_analysis(self,
                                     hamiltonian_matrix: csc_matrix,
                                     quri_hamiltonian: Operator, 
                                     n_qubits: int,
                                     approximate_states: Optional[List[GeneralCircuitQuantumState]] = None,
                                     test_config: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Perform comprehensive comparative analysis between exact and approximate initial states.
        
        This method runs the full comparison to demonstrate the impact of state preparation
        accuracy on QSCI algorithm performance.
        
        Args:
            hamiltonian_matrix: Sparse Hamiltonian matrix
            quri_hamiltonian: QURI-Parts Operator format Hamiltonian
            n_qubits: Number of qubits
            approximate_states: Optional list of approximate states for comparison
            test_config: Optional configuration overrides
            
        Returns:
            Comprehensive analysis results
        """
        if self.verbose:
            print(f"=== COMPARATIVE PRECISION ANALYSIS ===")
        
        # Test with exact ground state
        exact_results = self.test_all_variants_with_exact_ground_state(
            hamiltonian_matrix, quri_hamiltonian, n_qubits, test_config
        )
        
        # Test with approximate states if provided
        approximate_results = {}
        if approximate_states:
            for variant_name in exact_results['variants']:
                if exact_results['variants'][variant_name]['success']:
                    approx_result = self.test_variant_with_approximate_states(
                        hamiltonian_matrix, quri_hamiltonian, n_qubits,
                        approximate_states, variant_name, test_config
                    )
                    approximate_results[variant_name] = approx_result
        
        # Comparative analysis
        comparison = {
            'exact_results': exact_results,
            'approximate_results': approximate_results,
            'improvement_analysis': self._analyze_improvement_factors(
                exact_results, approximate_results
            ),
            'conclusions': self._generate_conclusions(exact_results, approximate_results)
        }
        
        if self.verbose:
            self._print_comparative_analysis(comparison)
        
        return comparison
    
    def _test_vanilla_qsci(self, hamiltonian, initial_state, exact_ground_energy, config):
        """Test VanillaQSCI with given initial state."""
        start_time = time.time()
        
        sampler = create_qulacs_vector_concurrent_sampler()
        vanilla_qsci = VanillaQSCI(
            hamiltonian=hamiltonian,
            sampler=sampler,
            num_states_pick_out=config.get('num_states', self.default_num_states)
        )
        
        result = vanilla_qsci.run([initial_state], total_shots=config.get('shots', self.default_shots))
        
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        execution_time = time.time() - start_time
        
        return {
            'success': True,
            'energy': result.ground_state_energy,
            'error': energy_error,
            'eigenvalues': result.eigenvalues,
            'execution_time': execution_time,
            'precision_achieved': energy_error < self.precision_tolerance
        }
    
    def _test_single_time_te_qsci(self, hamiltonian, initial_state, exact_ground_energy, config):
        """Test SingleTimeTE_QSCI with given initial state."""
        start_time = time.time()
        
        evolution_time = config.get('evolution_time', 0.5)
        sampler = create_qulacs_vector_concurrent_sampler()
        
        st_te_qsci = SingleTimeTE_QSCI(
            hamiltonian=hamiltonian,
            sampler=sampler,
            evolution_time=evolution_time,
            num_states_pick_out=config.get('num_states', self.default_num_states)
        )
        
        result = st_te_qsci.run([initial_state], total_shots=config.get('shots', self.default_shots))
        
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        execution_time = time.time() - start_time
        
        return {
            'success': True,
            'energy': result.ground_state_energy,
            'error': energy_error,
            'eigenvalues': result.eigenvalues,
            'execution_time': execution_time,
            'evolution_time': evolution_time,
            'precision_achieved': energy_error < self.precision_tolerance
        }
    
    def _test_state_vector_te_qsci(self, hamiltonian, initial_state, exact_ground_energy, config):
        """Test StateVectorTE_QSCI with given initial state."""
        start_time = time.time()
        
        evolution_time = config.get('evolution_time', 0.5)
        sampler = create_qulacs_vector_concurrent_sampler()
        
        sv_te_qsci = StateVectorTE_QSCI(
            hamiltonian=hamiltonian,
            sampler=sampler,
            evolution_time=evolution_time,
            num_states_pick_out=config.get('num_states', self.default_num_states)
        )
        
        result = sv_te_qsci.run([initial_state], total_shots=config.get('shots', self.default_shots))
        
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        execution_time = time.time() - start_time
        
        return {
            'success': True,
            'energy': result.ground_state_energy,
            'error': energy_error,
            'eigenvalues': result.eigenvalues,
            'execution_time': execution_time,
            'evolution_time': evolution_time,
            'precision_achieved': energy_error < self.precision_tolerance
        }
    
    def _test_time_average_te_qsci(self, hamiltonian, initial_state, exact_ground_energy, config):
        """Test TimeAverageTE_QSCI with given initial state."""
        start_time = time.time()
        
        evolution_times = config.get('evolution_times', [0.3, 0.7, 1.0])
        shots_per_time = config.get('shots_per_time', self.default_shots // len(evolution_times))
        sampler = create_qulacs_vector_concurrent_sampler()
        
        ta_te_qsci = TimeAverageTE_QSCI(
            hamiltonian=hamiltonian,
            sampler=sampler,
            evolution_times=evolution_times,
            num_states_pick_out=config.get('num_states', self.default_num_states)
        )
        
        result = ta_te_qsci.run([initial_state], shots_per_time=shots_per_time)
        
        energy_error = abs(result.ground_state_energy - exact_ground_energy)
        execution_time = time.time() - start_time
        
        return {
            'success': True,
            'energy': result.ground_state_energy,
            'error': energy_error,
            'eigenvalues': result.eigenvalues,
            'execution_time': execution_time,
            'evolution_times': evolution_times,
            'shots_per_time': shots_per_time,
            'precision_achieved': energy_error < self.precision_tolerance
        }
    
    def _find_best_variant(self, variants_results):
        """Find the QSCI variant with the best (lowest) error."""
        best_variant = None
        best_error = float('inf')
        
        for variant_name, result in variants_results.items():
            if result.get('success', False):
                error = result.get('error', float('inf'))
                if error < best_error:
                    best_error = error
                    best_variant = variant_name
        
        return {
            'name': best_variant,
            'error': best_error
        }
    
    def _check_precision_achieved(self, variants_results):
        """Check which variants achieved the target precision."""
        precision_achieved = {}
        
        for variant_name, result in variants_results.items():
            if result.get('success', False):
                precision_achieved[variant_name] = result.get('precision_achieved', False)
        
        return precision_achieved
    
    def _analyze_approximate_state_results(self, approximate_results):
        """Analyze results from approximate state testing."""
        successful_tests = [r for r in approximate_results.values() if r.get('success', False)]
        
        if not successful_tests:
            return {'num_successful': 0, 'avg_error': float('inf')}
        
        errors = [r['error'] for r in successful_tests]
        
        return {
            'num_successful': len(successful_tests),
            'avg_error': np.mean(errors),
            'min_error': np.min(errors),
            'max_error': np.max(errors),
            'std_error': np.std(errors)
        }
    
    def _analyze_improvement_factors(self, exact_results, approximate_results):
        """Analyze improvement factors between exact and approximate states."""
        improvement_analysis = {}
        
        for variant_name in exact_results['variants']:
            if (exact_results['variants'][variant_name].get('success', False) and
                variant_name in approximate_results):
                
                exact_error = exact_results['variants'][variant_name]['error']
                approx_summary = approximate_results[variant_name]['summary']
                
                if approx_summary.get('num_successful', 0) > 0:
                    avg_approx_error = approx_summary['avg_error']
                    improvement_factor = avg_approx_error / max(exact_error, 1e-16)
                    
                    improvement_analysis[variant_name] = {
                        'exact_error': exact_error,
                        'avg_approximate_error': avg_approx_error,
                        'improvement_factor': improvement_factor
                    }
        
        return improvement_analysis
    
    def _generate_conclusions(self, exact_results, approximate_results):
        """Generate high-level conclusions from the analysis."""
        conclusions = []
        
        # Check if any variant achieved high precision with exact states
        exact_precision_variants = [
            name for name, result in exact_results['variants'].items()
            if result.get('success', False) and result.get('precision_achieved', False)
        ]
        
        if exact_precision_variants:
            conclusions.append(
                f"✅ QSCI algorithms are mathematically correct: "
                f"{exact_precision_variants} achieved high precision with exact ground states"
            )
        else:
            conclusions.append(
                "⚠️ No QSCI variant achieved target precision even with exact ground states"
            )
        
        # Analyze state preparation impact
        if approximate_results:
            conclusions.append(
                "📊 State preparation quality significantly impacts QSCI accuracy"
            )
        
        return conclusions
    
    def _print_summary(self, summary):
        """Print test summary."""
        print(f"\n=== TEST SUMMARY ===")
        print(f"Successful variants: {summary['successful_variants']}")
        print(f"Failed variants: {summary['failed_variants']}")
        print(f"Success rate: {summary['success_rate']:.1%}")
        
        if summary['best_variant']['name']:
            print(f"Best variant: {summary['best_variant']['name']} "
                  f"(error: {summary['best_variant']['error']:.2e})")
        
        print(f"Total test time: {summary['total_test_time']:.2f}s")
    
    def _print_comparative_analysis(self, comparison):
        """Print comparative analysis results."""
        print(f"\n=== COMPARATIVE ANALYSIS CONCLUSIONS ===")
        for conclusion in comparison['conclusions']:
            print(f"{conclusion}")
        
        if comparison['improvement_analysis']:
            print(f"\n=== IMPROVEMENT FACTORS ===")
            for variant, analysis in comparison['improvement_analysis'].items():
                print(f"{variant}: {analysis['improvement_factor']:.1f}x improvement "
                      f"with exact vs approximate states")