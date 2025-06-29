"""
Main interface for ffsim/UCJ ansatz integration with QSCI framework.

This module provides high-level functions to run QSCI using UCJ/LUCJ ansatz
as input states, combining all the components into a seamless workflow.

This module requires ffsim to be installed: pip install quri-qsci[ffsim]
"""

from typing import List, Dict, Optional, Tuple
import numpy as np
from dataclasses import dataclass
import time

# ffsim imports - with proper error handling
try:
    import ffsim
    FFSIM_AVAILABLE = True
except ImportError:
    FFSIM_AVAILABLE = False
    
def _require_ffsim():
    """Check if ffsim is available and raise helpful error if not."""
    if not FFSIM_AVAILABLE:
        raise ImportError(
            "ffsim is required for this functionality. "
            "Install with: pip install quri-qsci[ffsim]"
        )

# QSCI imports from parent package
from qsci_algorithms import VanillaQSCI, QSCIResult

# Local imports
from .molecular_systems import MolecularSystem, create_h2_molecule, create_n2_molecule
from .integration import create_ucj_ansatz, create_lucj_ansatz, UCJResult
from .state_conversion import ucj_result_to_quri_state, ConversionMetrics


@dataclass
class LUCJQSCIResult:
    """Result from LUCJ/UCJ + QSCI workflow."""
    molecule_name: str
    ansatz_type: str
    ansatz_energy: float
    ansatz_optimization_success: bool
    
    # QSCI results for different subspace sizes
    qsci_results: Dict[int, QSCIResult]
    
    # State conversion metrics
    conversion_metrics: ConversionMetrics
    
    # Benchmark data
    hartree_fock_energy: float
    fci_energy: float
    target_energy: Optional[float]
    
    # Performance metrics
    total_time: float
    ansatz_time: float
    conversion_time: float
    qsci_time: float
    
    @property
    def best_qsci_energy(self) -> float:
        """Best (lowest) energy from QSCI calculations."""
        if not self.qsci_results:
            return self.ansatz_energy
        return min(result.ground_state_energy for result in self.qsci_results.values())
    
    @property
    def best_subspace_size(self) -> int:
        """Subspace size that gave the best energy."""
        if not self.qsci_results:
            return 0
        best_size = min(self.qsci_results.keys(), 
                       key=lambda k: self.qsci_results[k].ground_state_energy)
        return best_size
    
    @property
    def energy_improvement_vs_ansatz(self) -> float:
        """Energy improvement of best QSCI result vs ansatz."""
        return self.ansatz_energy - self.best_qsci_energy
    
    @property
    def energy_error_vs_fci(self) -> float:
        """Energy error of best QSCI result vs FCI reference."""
        return self.best_qsci_energy - self.fci_energy


def run_lucj_qsci(
    molecule: str,
    ansatz_type: str = "LUCJ",
    subspace_sizes: List[int] = [50, 100, 150],
    basis: str = "sto-3g",
    bond_length: Optional[float] = None,
    n_reps: int = 1,
    max_optimization_iterations: int = 50,
    total_shots: int = 5000,
    conversion_method: str = "sampling_circuit",
    active_space: Optional[Tuple[int, int]] = None,
    use_homo_lumo: bool = True
) -> LUCJQSCIResult:
    """Run complete LUCJ/UCJ + QSCI workflow.
    
    Args:
        molecule: Molecule name ("H2" or "N2")
        ansatz_type: "UCJ" or "LUCJ" 
        subspace_sizes: List of subspace dimensions to test
        basis: Basis set for quantum chemistry
        bond_length: Bond length (None for default)
        n_reps: Number of ansatz repetitions
        max_optimization_iterations: Max iterations for ansatz optimization
        total_shots: Total measurement shots for QSCI
        conversion_method: State conversion method
        active_space: (n_electrons, n_orbitals) for N2 active space (None for full space)
        use_homo_lumo: If True, use HOMO-LUMO focused active space for N2 (more efficient)
        
    Returns:
        LUCJQSCIResult with complete workflow results
    """
    print(f"=" * 60)
    print(f"Running {ansatz_type} + QSCI workflow for {molecule}")
    print(f"=" * 60)
    
    start_time = time.time()
    
    # 1. Create molecular system
    print(f"Step 1: Creating {molecule} molecular system...")
    ansatz_start = time.time()
    
    if molecule.upper() == "H2":
        bond_length = bond_length or 0.74
        mol_system = create_h2_molecule(basis=basis, bond_length=bond_length)
    elif molecule.upper() == "N2":
        bond_length = bond_length or 1.0
        if use_homo_lumo:
            # Use HOMO-LUMO focused active space for efficient optimization
            mol_system = create_n2_homo_lumo_molecule(
                basis=basis,
                bond_length=bond_length,
                active_space_size=4
            )
        elif active_space is not None:
            n_frozen = active_space[0] if len(active_space) > 0 else 2
            mol_system = create_n2_active_space_molecule(
                n_frozen=n_frozen,
                basis=basis,
                bond_length=bond_length
            )
        else:
            mol_system = create_n2_molecule(basis=basis, bond_length=bond_length)
    else:
        raise ValueError(f"Unsupported molecule: {molecule}")
    
    # 2. Generate ansatz state
    print(f"\nStep 2: Generating {ansatz_type} ansatz...")
    
    if ansatz_type.upper() == "UCJ":
        ansatz_result = create_ucj_ansatz(
            mol_system, 
            n_reps=n_reps, 
            max_iterations=max_optimization_iterations
        )
    elif ansatz_type.upper() == "LUCJ":
        ansatz_result = create_lucj_ansatz(
            mol_system, 
            n_reps=n_reps, 
            max_iterations=max_optimization_iterations
        )
    else:
        raise ValueError(f"Unsupported ansatz type: {ansatz_type}")
    
    ansatz_time = time.time() - ansatz_start
    
    # 3. Convert to QURI Parts format
    print(f"\nStep 3: Converting {ansatz_type} state to QURI Parts format...")
    conversion_start = time.time()
    
    # Jordan-Wigner mapping: each spatial orbital needs 2 qubits (alpha and beta spin)
    # This gives us the full computational basis space needed for QSCI
    n_qubits = 2 * mol_system.ffsim_mol_data.norb
    quri_state, conversion_metrics = ucj_result_to_quri_state(
        ansatz_result, n_qubits, conversion_method=conversion_method
    )
    
    conversion_time = time.time() - conversion_start
    
    # 4. Run QSCI with different subspace sizes
    print(f"\nStep 4: Running VanillaQSCI with different subspace sizes...")
    qsci_start = time.time()
    
    qsci_results = {}
    
    for subspace_size in subspace_sizes:
        print(f"\n  Running QSCI with subspace size {subspace_size}...")
        
        # Create QSCI algorithm
        qsci_algo = VanillaQSCI(
            hamiltonian=mol_system.quri_hamiltonian,
            sampler=None,  # Will be created automatically
            num_states_pick_out=subspace_size
        )
        
        # Create sampler
        from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
        qsci_algo.sampler = create_qulacs_vector_concurrent_sampler()
        
        # Run QSCI
        try:
            qsci_result = qsci_algo.run(
                input_states=[quri_state],
                total_shots=total_shots,
                start_time=time.time()
            )
            qsci_results[subspace_size] = qsci_result
            
            print(f"    ✓ QSCI R={subspace_size}: Energy = {qsci_result.ground_state_energy:.6f} Ha")
            
        except Exception as e:
            print(f"    ✗ QSCI R={subspace_size} failed: {e}")
            continue
    
    qsci_time = time.time() - qsci_start
    total_time = time.time() - start_time
    
    # 5. Create result summary
    print(f"\nStep 5: Summarizing results...")
    
    # Determine target energy for benchmark
    target_energy = None
    if molecule.upper() == "N2" and basis == "6-31g":
        target_energy = -109.0  # Target benchmark energy
    
    result = LUCJQSCIResult(
        molecule_name=molecule.upper(),
        ansatz_type=ansatz_type.upper(),
        ansatz_energy=ansatz_result.final_energy,
        ansatz_optimization_success=ansatz_result.optimization_success,
        qsci_results=qsci_results,
        conversion_metrics=conversion_metrics,
        hartree_fock_energy=mol_system.hartree_fock_energy,
        fci_energy=mol_system.fci_energy,
        target_energy=target_energy,
        total_time=total_time,
        ansatz_time=ansatz_time,
        conversion_time=conversion_time,
        qsci_time=qsci_time
    )
    
    # Print summary
    print_result_summary(result)
    
    return result


def print_result_summary(result: LUCJQSCIResult):
    """Print a comprehensive summary of LUCJ/UCJ + QSCI results."""
    print(f"\n" + "=" * 60)
    print(f"WORKFLOW SUMMARY: {result.ansatz_type} + QSCI for {result.molecule_name}")
    print(f"=" * 60)
    
    print(f"\n📊 ENERGY RESULTS:")
    print(f"  Hartree-Fock energy:     {result.hartree_fock_energy:.6f} Ha")
    print(f"  FCI reference energy:    {result.fci_energy:.6f} Ha")
    print(f"  {result.ansatz_type} ansatz energy:      {result.ansatz_energy:.6f} Ha")
    print(f"  Best QSCI energy:        {result.best_qsci_energy:.6f} Ha")
    
    if result.target_energy:
        print(f"  Target benchmark energy: {result.target_energy:.6f} Ha")
        error_vs_target = result.best_qsci_energy - result.target_energy
        print(f"  Error vs target:         {error_vs_target:.6f} Ha ({error_vs_target*1000:.1f} mHa)")
    
    print(f"\n📈 PERFORMANCE METRICS:")
    print(f"  {result.ansatz_type} vs HF improvement:  {result.hartree_fock_energy - result.ansatz_energy:.6f} Ha")
    print(f"  QSCI vs {result.ansatz_type} improvement: {result.energy_improvement_vs_ansatz:.6f} Ha")
    print(f"  Error vs FCI:            {result.energy_error_vs_fci:.6f} Ha")
    print(f"  Best subspace size:      R = {result.best_subspace_size}")
    
    print(f"\n🔄 CONVERSION QUALITY:")
    print(f"  State conversion fidelity: {result.conversion_metrics.fidelity:.4f}")
    print(f"  Probability overlap:       {result.conversion_metrics.probability_overlap:.4f}")
    
    print(f"\n⏱️  TIMING BREAKDOWN:")
    print(f"  {result.ansatz_type} optimization time:  {result.ansatz_time:.2f} seconds")
    print(f"  State conversion time:    {result.conversion_time:.2f} seconds")
    print(f"  QSCI calculation time:    {result.qsci_time:.2f} seconds")
    print(f"  Total workflow time:      {result.total_time:.2f} seconds")
    
    print(f"\n📊 SUBSPACE CONVERGENCE:")
    if result.qsci_results:
        for size in sorted(result.qsci_results.keys()):
            energy = result.qsci_results[size].ground_state_energy
            error = energy - result.fci_energy
            print(f"  R={size:3d}: {energy:.6f} Ha (error: {error:.6f} Ha)")
    
    print(f"\n" + "=" * 60)


def run_convergence_study(
    molecule: str,
    ansatz_type: str = "LUCJ",
    max_subspace: int = 200,
    step_size: int = 25,
    **kwargs
) -> LUCJQSCIResult:
    """Run a convergence study with increasing subspace sizes.
    
    Args:
        molecule: Molecule name
        ansatz_type: Ansatz type
        max_subspace: Maximum subspace size to test
        step_size: Step size for subspace sizes
        **kwargs: Additional arguments for run_lucj_qsci
        
    Returns:
        LUCJQSCIResult with convergence data
    """
    subspace_sizes = list(range(25, max_subspace + 1, step_size))
    
    print(f"Running convergence study for {ansatz_type} + QSCI")
    print(f"Subspace sizes: {subspace_sizes}")
    
    return run_lucj_qsci(
        molecule=molecule,
        ansatz_type=ansatz_type,
        subspace_sizes=subspace_sizes,
        **kwargs
    )


def benchmark_against_reference(
    result: LUCJQSCIResult, 
    reference_energy: float,
    tolerance: float = 0.001
) -> Dict[str, bool]:
    """Benchmark LUCJ/UCJ + QSCI results against reference energy.
    
    Args:
        result: LUCJQSCIResult to benchmark
        reference_energy: Reference energy in Hartree
        tolerance: Energy tolerance in Hartree
        
    Returns:
        Dictionary with benchmark results
    """
    best_energy = result.best_qsci_energy
    energy_error = abs(best_energy - reference_energy)
    
    benchmark = {
        "energy_within_tolerance": energy_error < tolerance,
        "convergence_demonstrated": len(result.qsci_results) > 1,
        "ansatz_improves_hf": result.ansatz_energy < result.hartree_fock_energy,
        "qsci_improves_ansatz": result.best_qsci_energy < result.ansatz_energy,
        "reasonable_conversion_fidelity": result.conversion_metrics.fidelity > 0.5
    }
    
    print(f"\nBENCHMARK RESULTS vs {reference_energy:.6f} Ha:")
    for criterion, passed in benchmark.items():
        status = "✓" if passed else "✗"
        print(f"  {status} {criterion.replace('_', ' ').title()}")
    
    all_passed = all(benchmark.values())
    print(f"\nOverall benchmark: {'✓ PASSED' if all_passed else '✗ FAILED'}")
    
    return benchmark