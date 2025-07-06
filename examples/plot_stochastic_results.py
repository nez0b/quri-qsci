#!/usr/bin/env python3
"""
Plotting script for stochastic QSCI results - Figure S3 reproduction.

This script reads the results from run_stochastic_qsci.py and creates
publication-quality plots showing energy error convergence vs number
of samples, reproducing the style of Figure S3 from arxiv:2501.07231.

Usage:
    python plot_stochastic_results.py --files n2_10e18o_1_09A_stochastic_results.npy
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def load_results(filepath):
    """Load stochastic QSCI results from file."""
    results = np.load(filepath, allow_pickle=True).item()
    
    # Extract data for plotting
    sample_sizes = sorted([k for k in results.keys() if isinstance(k, (int, np.integer))])
    avg_errors = []
    std_errors = []
    
    for size in sample_sizes:
        res = results[size]
        avg_errors.append(res['avg_error'])
        std_errors.append(res['std_error'])
    
    return sample_sizes, avg_errors, std_errors

def plot_stochastic_convergence(result_files, labels=None, colors=None, output_dir="figures"):
    """
    Plot stochastic QSCI convergence results.
    
    Args:
        result_files: List of paths to result files
        labels: List of labels for each file (optional)
        colors: List of colors for each file (optional)
        output_dir: Directory to save plots
    """
    # Setup plotting
    plt.figure(figsize=(10, 7))
    
    # Default colors and labels if not provided
    if colors is None:
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    if labels is None:
        labels = [f"System {i+1}" for i in range(len(result_files))]
    
    # Plot each result file
    for i, (filepath, label, color) in enumerate(zip(result_files, labels, colors)):
        try:
            sample_sizes, avg_errors, std_errors = load_results(filepath)
            
            # Convert to numpy arrays for easier handling
            sample_sizes = np.array(sample_sizes)
            avg_errors = np.array(avg_errors)
            std_errors = np.array(std_errors)
            
            # Filter out NaN values
            valid_mask = ~np.isnan(avg_errors)
            if np.sum(valid_mask) == 0:
                print(f"Warning: No valid data in {filepath}")
                continue
                
            sample_sizes_valid = sample_sizes[valid_mask]
            avg_errors_valid = avg_errors[valid_mask]
            std_errors_valid = std_errors[valid_mask]
            
            # Plot with error bars
            plt.errorbar(sample_sizes_valid, avg_errors_valid, yerr=std_errors_valid,
                        label=label, marker='o', capsize=4, capthick=1.5,
                        linewidth=2, markersize=6, color=color)
            
            print(f"✓ Plotted {len(sample_sizes_valid)} data points for {label}")
            
        except Exception as e:
            print(f"❌ Error plotting {filepath}: {e}")
            continue
    
    # Format plot
    plt.xscale('log')
    plt.xlabel('Number of Samples', fontsize=14)
    plt.ylabel('Energy Error (mHa)', fontsize=14)
    plt.title('Stochastic QSCI Convergence Analysis\nN2/cc-pVDZ System', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.3)
    plt.axhline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Improve layout
    plt.tight_layout()
    
    # Save plots
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save as PNG
    png_path = output_path / "stochastic_qsci_convergence.png"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved plot: {png_path}")
    
    # Save as PDF
    pdf_path = output_path / "stochastic_qsci_convergence.pdf"
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"✓ Saved plot: {pdf_path}")
    
    plt.show()
    
    return png_path, pdf_path

def print_convergence_analysis(result_files):
    """Print detailed convergence analysis."""
    print(f"\n" + "="*80)
    print("CONVERGENCE ANALYSIS")
    print("="*80)
    
    for filepath in result_files:
        try:
            sample_sizes, avg_errors, std_errors = load_results(filepath)
            
            print(f"\nFile: {filepath}")
            print(f"{'Sample Size':<12} {'Error (mHa)':<12} {'Std (mHa)':<12} {'Improvement':<12}")
            print("-" * 60)
            
            prev_error = None
            for size, error, std in zip(sample_sizes, avg_errors, std_errors):
                if not np.isnan(error):
                    if prev_error is not None and not np.isnan(prev_error):
                        improvement = abs(prev_error) - abs(error)
                        improvement_str = f"{improvement:+.3f}"
                    else:
                        improvement_str = "N/A"
                    
                    print(f"{size:<12,} {error:<12.3f} {std:<12.3f} {improvement_str:<12}")
                    prev_error = error
                else:
                    print(f"{size:<12,} {'FAILED':<12} {'FAILED':<12} {'N/A':<12}")
                    
        except Exception as e:
            print(f"❌ Error analyzing {filepath}: {e}")

def extract_system_info(filepath):
    """Extract system information from filename for automatic labeling."""
    filename = Path(filepath).stem
    
    # Try to parse standard naming convention: n2_10e18o_1_09A_stochastic_results
    parts = filename.split('_')
    
    if len(parts) >= 4:
        molecule = parts[0].upper()
        active_space = parts[1] + parts[2]  # e.g., "10e18o"
        bond_length = parts[3].replace('_', '.')  # e.g., "1.09A"
        
        return f"{molecule} ({active_space}) at {bond_length}"
    else:
        return Path(filepath).stem

def main():
    """Main function to plot stochastic QSCI results."""
    parser = argparse.ArgumentParser(description="Plot stochastic QSCI convergence results")
    parser.add_argument('--files', type=str, nargs='+', required=True,
                       help='Paths to stochastic results files (.npy)')
    parser.add_argument('--labels', type=str, nargs='*',
                       help='Custom labels for each file')
    parser.add_argument('--output-dir', type=str, default="figures",
                       help='Output directory for plots')
    parser.add_argument('--analysis', action='store_true',
                       help='Print detailed convergence analysis')
    
    args = parser.parse_args()
    
    print("="*80)
    print("STOCHASTIC QSCI PLOTTING - FIGURE S3 REPRODUCTION")
    print("="*80)
    
    # Validate input files
    valid_files = []
    for filepath in args.files:
        path = Path(filepath)
        if not path.exists():
            # Try adding stochastic_data prefix
            alt_path = Path("stochastic_data") / path.name
            if alt_path.exists():
                valid_files.append(str(alt_path))
                print(f"✓ Found file: {alt_path}")
            else:
                print(f"❌ File not found: {filepath}")
        else:
            valid_files.append(filepath)
            print(f"✓ Found file: {filepath}")
    
    if not valid_files:
        print("❌ No valid result files found!")
        return False
    
    # Generate labels if not provided
    if args.labels:
        if len(args.labels) != len(valid_files):
            print("❌ Number of labels must match number of files!")
            return False
        labels = args.labels
    else:
        labels = [extract_system_info(f) for f in valid_files]
    
    print(f"\nPlotting {len(valid_files)} result file(s):")
    for f, l in zip(valid_files, labels):
        print(f"  {f} → {l}")
    
    try:
        # Create plots
        png_path, pdf_path = plot_stochastic_convergence(
            valid_files, labels=labels, output_dir=args.output_dir
        )
        
        # Print analysis if requested
        if args.analysis:
            print_convergence_analysis(valid_files)
        
        print(f"\n✓ Plotting completed successfully!")
        print(f"  PNG: {png_path}")
        print(f"  PDF: {pdf_path}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in plotting: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)