"""
ffsim Integration Module for QSCI

This module provides integration between ffsim's UCJ/LUCJ ansatz and the QSCI framework.
It requires ffsim to be installed: pip install quri-qsci[ffsim]

Key Components:
- molecular_systems: Molecular system creation utilities
- integration: ffsim UCJ/LUCJ ansatz integration
- state_conversion: State format conversion between ffsim and QURI Parts
- qsci_interface: High-level QSCI-ffsim integration interface
"""

# Check if ffsim is available
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

# Conditional imports - only import if ffsim is available
if FFSIM_AVAILABLE:
    from .molecular_systems import (
        MolecularSystem,
        create_h2_molecule,
        create_n2_molecule,
        get_reference_energies,
        validate_molecular_system,
    )
    
    from .integration import (
        UCJResult,
        create_ucj_ansatz,
        create_lucj_ansatz,
    )
    
    from .state_conversion import (
        ConversionMetrics,
        ffsim_to_quri_state,
        ucj_result_to_quri_state,
        _map_fermionic_amplitudes_directly,
    )
    
    from .qsci_interface import (
        LUCJQSCIResult,
        run_lucj_qsci,
        print_result_summary,
        run_convergence_study,
        benchmark_against_reference,
    )
    
    __all__ = [
        # Molecular systems
        "MolecularSystem",
        "create_h2_molecule", 
        "create_n2_molecule",
        "get_reference_energies",
        "validate_molecular_system",
        
        # ffsim integration
        "UCJResult",
        "create_ucj_ansatz",
        "create_lucj_ansatz", 
        
        # State conversion
        "ConversionMetrics",
        "ffsim_to_quri_state",
        "ucj_result_to_quri_state",
        "_map_fermionic_amplitudes_directly",
        
        # QSCI interface
        "LUCJQSCIResult",
        "run_lucj_qsci",
        "print_result_summary",
        "run_convergence_study",
        "benchmark_against_reference",
        
        # Utilities
        "FFSIM_AVAILABLE",
    ]
    
else:
    # Define what happens when ffsim is not available
    def __getattr__(name):
        """Provide helpful error message when ffsim functionality is accessed without ffsim installed."""
        _require_ffsim()
    
    __all__ = ["FFSIM_AVAILABLE"]

__version__ = "0.1.0"