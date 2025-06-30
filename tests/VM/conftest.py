"""
Pytest configuration and shared fixtures for VM tests.

Provides common test utilities, mock configurations, and 
shared fixtures for comprehensive VM module testing.
"""

import pytest
from unittest.mock import Mock, patch
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))


@pytest.fixture(autouse=True)
def mock_external_dependencies():
    """Auto-mock external dependencies to isolate VM tests."""
    with patch.dict('sys.modules', {
        'quri_vm': Mock(),
        'quri_parts.backend.devices.star_device': Mock(),
        'quri_algo.circuit.time_evolution.trotter_time_evo': Mock(),
        'quri_algo.problem.operators.hamiltonian': Mock(),
    }):
        yield


@pytest.fixture
def clean_imports():
    """Ensure clean import state for each test."""
    # Store original modules
    original_modules = sys.modules.copy()
    yield
    # Restore original modules
    sys.modules.clear()
    sys.modules.update(original_modules)


@pytest.fixture
def disable_print_output():
    """Disable print output during tests for cleaner output."""
    with patch('builtins.print'):
        yield