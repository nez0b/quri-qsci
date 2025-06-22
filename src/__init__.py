"""
TE-QSCI (Time-Evolved Quantum Selected Configuration Interaction) package.

This package implements the TE-QSCI algorithm and its variants using the QURI Parts framework.
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

__version__ = "0.1.0"
__author__ = "TE-QSCI Implementation"
__all__ = [
    # Core algorithms
    "QSCIResult", "QSCIVariant", "VanillaQSCI", "TimeEvolvedQSCI",
    # Algorithm interfaces
    "QSCIAlgorithmResult", "QSCIAlgorithmBase", 
    "VanillaQSCIAlgorithm", "SingleTimeTeQSCIAlgorithm", 
    "TimeAverageTeQSCIAlgorithm", "StateVectorTeQSCIAlgorithm",
    "create_qsci_algorithm", "LoweringLevel",
    # VM analysis
    "CircuitResourceEstimate", "VMCircuitAnalyzer", "create_vm_enabled_algorithm"
]