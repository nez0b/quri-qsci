"""
QSCI (Quantum Selected Configuration Interaction) package with optional ffsim integration.

This package implements the QSCI algorithm and its variants using the QURI Parts framework.
Optional ffsim integration provides UCJ/LUCJ ansatz support.
"""

from .qsci_algorithms import (
    QSCIResult,
    QSCIVariant, 
    VanillaQSCI,
    TimeEvolvedQSCI
)

from .qsci_algo_interface import (
    QSCIAlgorithmResult,
    QSCIAlgorithmBase,
    VanillaQSCIAlgorithm,
    SingleTimeTeQSCIAlgorithm,
    TimeAverageTeQSCIAlgorithm,
    StateVectorTeQSCIAlgorithm,
    create_qsci_algorithm,
    LoweringLevel
)

from .qsci_vm_analysis import (
    CircuitResourceEstimate,
    VMCircuitAnalyzer,
    create_vm_enabled_algorithm
)

# Optional ffsim module - only available if ffsim is installed
try:
    from . import ffsim_integration
    FFSIM_AVAILABLE = True
except ImportError:
    ffsim_integration = None
    FFSIM_AVAILABLE = False

__version__ = "0.1.0"
__author__ = "QSCI Development Team"
__all__ = [
    # Core algorithms
    "QSCIResult", "QSCIVariant", "VanillaQSCI", "TimeEvolvedQSCI",
    # Algorithm interfaces
    "QSCIAlgorithmResult", "QSCIAlgorithmBase", 
    "VanillaQSCIAlgorithm", "SingleTimeTeQSCIAlgorithm", 
    "TimeAverageTeQSCIAlgorithm", "StateVectorTeQSCIAlgorithm",
    "create_qsci_algorithm", "LoweringLevel",
    # VM analysis
    "CircuitResourceEstimate", "VMCircuitAnalyzer", "create_vm_enabled_algorithm",
    # Optional ffsim module
    "ffsim_integration", "FFSIM_AVAILABLE"
]