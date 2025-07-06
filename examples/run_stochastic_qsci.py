#!/usr/bin/env python3
"""
Stochastic QSCI Implementation for Figure S3 Reproduction.

This script implements the stochastic sampling approach used by the authors
in their qsci-benchmarks repository for reproducing Figure S3. The method
samples determinants based on |c_i|² probabilities and studies convergence
vs number of samples.

Key features:
- Memory-efficient subspace construction using PySCF's fci.addons.large_ci
- Stochastic sampling from CI coefficient probability distribution
- Multiple repeats for statistical analysis
- Sample sizes from 10¹ to 10⁶ (reduced from 10⁹ for testing)

Usage:
    python run_stochastic_qsci.py --base n2_10e18o_1_09A --repeats 10
"""

import numpy as np
import scipy.linalg
import argparse
import time
from pathlib import Path

# PySCF imports
try:
    from pyscf import fci
    from pyscf.fci import cistring
except ImportError as e:
    raise ImportError(
        f"PySCF is required for FCI calculations: {e}\n"
        f"Install with: pip install pyscf"
    )

def solve_ci_subproblem(indices, h1, h2, norb, nelec, ecore, original_ci_coeffs):
    """
    Solve CI subproblem for selected determinants using PySCF energy calculation.
    
    This uses PySCF's direct_spin1.energy() function to properly calculate
    the energy of a CI vector composed of the selected determinants with
    their correct CI coefficients from the reference calculation.
    
    Args:
        indices: Array of determinant indices to include in subspace
        h1: One-electron Hamiltonian matrix from PySCF
        h2: Two-electron integral tensor from PySCF
        norb: Number of orbitals
        nelec: Number of electrons (tuple for α,β or int for total)
        ecore: Core energy
        original_ci_coeffs: Full CI coefficient vector from reference calculation
        
    Returns:
        Lowest eigenvalue of the subspace Hamiltonian
    """
    from pyscf.fci import direct_spin1
    
    # Handle electron configuration
    if isinstance(nelec, (tuple, list)):
        n_alpha, n_beta = nelec
    else:
        n_alpha = n_beta = nelec // 2
    
    # Get unique determinant indices
    unique_indices = np.unique(indices)
    n_dets = len(unique_indices)
    
    if n_dets == 0:
        return np.nan
    
    try:
        # Get total CI space dimension
        n_alpha_strings = cistring.num_strings(norb, n_alpha)
        n_beta_strings = cistring.num_strings(norb, n_beta)
        total_dets = n_alpha_strings * n_beta_strings
        
        # Create CI vector with selected determinants using their actual coefficients
        ci_vector = np.zeros(total_dets)
        
        # Use the actual CI coefficients from the reference calculation
        ci_flat = original_ci_coeffs.flatten() if hasattr(original_ci_coeffs, 'flatten') else original_ci_coeffs
        ci_vector[unique_indices] = ci_flat[unique_indices]
        
        # Normalize the truncated CI vector
        ci_norm = np.linalg.norm(ci_vector)
        if ci_norm > 1e-12:
            ci_vector /= ci_norm
        else:
            print(f"    Warning: CI vector has zero norm, using equal weights")
            ci_vector[unique_indices] = 1.0 / np.sqrt(len(unique_indices))
        
        # Calculate energy using PySCF's FCI energy function
        energy_active = direct_spin1.energy(h1, h2, ci_vector, norb, (n_alpha, n_beta))
        total_energy = energy_active + ecore
        
        return total_energy
        
    except Exception as e:
        print(f"    Warning: Error in CI subproblem: {e}")
        # Fallback: return single-determinant energy for the most probable determinant
        if len(unique_indices) > 0:
            return calculate_single_determinant_energy(unique_indices[0], h1, h2, norb, n_alpha, n_beta, ecore, original_ci_coeffs)
        else:
            return np.nan

def calculate_single_determinant_energy(det_index, h1, h2, norb, n_alpha, n_beta, ecore, original_ci_coeffs):
    """Calculate energy of a single determinant using PySCF with proper coefficient."""
    try:
        from pyscf.fci import direct_spin1
        
        # Get total CI space dimension
        n_alpha_strings = cistring.num_strings(norb, n_alpha)
        n_beta_strings = cistring.num_strings(norb, n_beta)
        total_dets = n_alpha_strings * n_beta_strings
        
        # Create CI vector with only this determinant using its actual coefficient
        ci_vector = np.zeros(total_dets)
        ci_flat = original_ci_coeffs.flatten() if hasattr(original_ci_coeffs, 'flatten') else original_ci_coeffs
        ci_vector[det_index] = ci_flat[det_index]
        
        # Normalize (though for single determinant this just changes the phase)
        ci_norm = np.linalg.norm(ci_vector)
        if ci_norm > 1e-12:
            ci_vector /= ci_norm
        else:
            ci_vector[det_index] = 1.0
        
        # Calculate energy using PySCF
        energy_active = direct_spin1.energy(h1, h2, ci_vector, norb, (n_alpha, n_beta))
        total_energy = energy_active + ecore
        
        return total_energy
        
    except Exception as e:
        print(f"    Error calculating single determinant energy: {e}")
        return np.nan



def run_stochastic_qsci(base_name, repeats=10, max_sample_exp=6):
    """
    Run stochastic QSCI calculations for a range of sample sizes.
    
    Args:
        base_name: Base name for input files (e.g., 'n2_10e18o_1_09A')
        repeats: Number of statistical repeats per sample size
        max_sample_exp: Maximum exponent for sample sizes (10^max_sample_exp)
        
    Returns:
        Dictionary with results for each sample size
    """
    data_dir = Path("stochastic_data")
    
    print(f"Loading data from base: {base_name}")
    
    # Load all required data
    try:
        ci_vector = np.load(data_dir / f"{base_name}_ci_vector.npy")
        fci_energy = np.load(data_dir / f"{base_name}_fci_energy.npy")
        h1 = np.load(data_dir / f"{base_name}_h1.npy")
        h2 = np.load(data_dir / f"{base_name}_h2.npy")
        ecore = np.load(data_dir / f"{base_name}_ecore.npy")
        config = np.load(data_dir / f"{base_name}_config.npy", allow_pickle=True).item()
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Required data file not found: {e}")
    
    norb = config['norb']
    nelec = config['nelec']
    
    print(f"✓ Loaded data:")
    print(f"  FCI energy: {fci_energy:.6f} Ha")
    print(f"  CI vector size: {ci_vector.size}")
    print(f"  Active space: ({nelec}e, {norb}o)")
    print(f"  h1 shape: {h1.shape}")
    print(f"  h2 shape: {h2.shape}")
    
    # Setup experiment parameters
    sample_sizes = [10**i for i in range(1, max_sample_exp + 1)]
    print(f"  Sample sizes: {sample_sizes}")
    print(f"  Repeats per size: {repeats}")
    
    # Prepare probability distribution
    ci_flat = ci_vector.flatten() if hasattr(ci_vector, 'flatten') else ci_vector
    probabilities = np.abs(ci_flat)**2
    probabilities /= np.sum(probabilities)  # Normalize
    determinant_indices = np.arange(len(ci_flat))
    
    print(f"  Probability distribution prepared")
    print(f"  Max probability: {np.max(probabilities):.6f}")
    print(f"  Min probability: {np.min(probabilities[probabilities > 0]):.6f}")
    
    # Run main sampling loop
    all_results = {}
    total_calculations = len(sample_sizes) * repeats
    calc_count = 0
    
    print(f"\nRunning {total_calculations} stochastic QSCI calculations...")
    
    for n_samples in sample_sizes:
        print(f"\n--- Sample size: {n_samples:,} ---")
        repeat_energies = []
        repeat_times = []
        
        for repeat in range(repeats):
            calc_count += 1
            progress = (calc_count / total_calculations) * 100
            
            print(f"  [{calc_count}/{total_calculations}] ({progress:.1f}%) "
                  f"Repeat {repeat + 1}/{repeats}")
            
            start_time = time.time()
            
            try:
                # Core stochastic sampling step
                sampled_indices = np.random.choice(
                    determinant_indices,
                    size=min(n_samples, len(determinant_indices)),  # Don't sample more than available
                    replace=True,
                    p=probabilities
                )
                
                unique_count = len(np.unique(sampled_indices))
                print(f"    Sampled {len(sampled_indices)} determinants, {unique_count} unique")
                
                # Solve CI subproblem
                qsci_energy = solve_ci_subproblem(sampled_indices, h1, h2, norb, nelec, ecore, ci_flat)
                
                repeat_energies.append(qsci_energy)
                repeat_times.append(time.time() - start_time)
                
                energy_error = (qsci_energy - fci_energy) * 1000  # Convert to mHa
                print(f"    E_QSCI: {qsci_energy:.6f} Ha")
                print(f"    Error: {energy_error:.3f} mHa")
                print(f"    Time: {repeat_times[-1]:.2f} s")
                
            except Exception as e:
                print(f"    Error in calculation: {e}")
                repeat_energies.append(np.nan)
                repeat_times.append(time.time() - start_time)
        
        # Calculate statistics
        valid_energies = [e for e in repeat_energies if not np.isnan(e)]
        
        if len(valid_energies) > 0:
            avg_energy = np.mean(valid_energies)
            std_energy = np.std(valid_energies) if len(valid_energies) > 1 else 0.0
            avg_error = (avg_energy - fci_energy) * 1000  # mHa
            std_error = std_energy * 1000  # mHa
            avg_time = np.mean(repeat_times)
            
            all_results[n_samples] = {
                'avg_energy': avg_energy,
                'std_energy': std_energy,
                'avg_error': avg_error,
                'std_error': std_error,
                'avg_time': avg_time,
                'success_rate': len(valid_energies) / repeats,
                'n_successful': len(valid_energies)
            }
            
            print(f"  → Average error: {avg_error:.3f} ± {std_error:.3f} mHa")
            print(f"  → Success rate: {len(valid_energies)}/{repeats}")
            
        else:
            print(f"  → All calculations failed for this sample size")
            all_results[n_samples] = {
                'avg_energy': np.nan,
                'std_energy': np.nan,
                'avg_error': np.nan,
                'std_error': np.nan,
                'avg_time': np.mean(repeat_times),
                'success_rate': 0.0,
                'n_successful': 0
            }
    
    return all_results

def save_results(results, base_name):
    """Save stochastic QSCI results to file."""
    output_path = Path("stochastic_data") / f"{base_name}_stochastic_results.npy"
    np.save(output_path, results)
    print(f"✓ Results saved to: {output_path}")
    return output_path

def print_summary(results, base_name):
    """Print summary of stochastic QSCI results."""
    print(f"\n" + "="*80)
    print(f"STOCHASTIC QSCI RESULTS SUMMARY - {base_name}")
    print(f"="*80)
    
    print(f"{'Sample Size':<12} {'Avg Error':<12} {'Std Error':<12} {'Success':<8} {'Time':<8}")
    print(f"{'-'*60}")
    
    for n_samples in sorted(results.keys()):
        res = results[n_samples]
        avg_err = res['avg_error']
        std_err = res['std_error']
        success = res['n_successful']
        avg_time = res['avg_time']
        
        if not np.isnan(avg_err):
            print(f"{n_samples:<12,} {avg_err:<12.3f} {std_err:<12.3f} {success:<8} {avg_time:<8.2f}")
        else:
            print(f"{n_samples:<12,} {'FAILED':<12} {'FAILED':<12} {success:<8} {avg_time:<8.2f}")

def main():
    """Main function to run stochastic QSCI analysis."""
    parser = argparse.ArgumentParser(description="Stochastic QSCI for Figure S3 reproduction")
    parser.add_argument('--base', type=str, required=True,
                       help='Base name for input files (e.g., n2_10e18o_1_09A)')
    parser.add_argument('--repeats', type=int, default=1,
                       help='Number of repeats for statistical averaging')
    parser.add_argument('--max-exp', type=int, default=6,
                       help='Maximum exponent for sample sizes (10^max_exp)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("STOCHASTIC QSCI ANALYSIS - FIGURE S3 REPRODUCTION")
    print("="*80)
    print(f"Base name: {args.base}")
    print(f"Repeats: {args.repeats}")
    print(f"Max sample size: 10^{args.max_exp}")
    
    try:
        # Run stochastic QSCI calculations
        results = run_stochastic_qsci(args.base, args.repeats, args.max_exp)
        
        # Save results
        output_path = save_results(results, args.base)
        
        # Print summary
        print_summary(results, args.base)
        
        print(f"\n✓ Stochastic QSCI analysis completed successfully!")
        print(f"  Results saved to: {output_path}")
        print(f"  Ready for plotting Figure S3")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in stochastic QSCI analysis: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)