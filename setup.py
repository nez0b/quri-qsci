#!/usr/bin/env python3
"""
Custom setup.py for quri-qsci with conditional ffsim_integration inclusion.

This setup script conditionally includes the ffsim_integration package
only when ffsim is available in the environment.
"""

from setuptools import setup, find_packages
import os
import sys


def get_packages():
    """Get list of packages to include, conditionally including ffsim_integration."""
    all_packages = find_packages(where="src")
    
    # Always include non-ffsim packages
    base_packages = [pkg for pkg in all_packages if not pkg.startswith('ffsim_integration')]
    
    # Check if ffsim should be included
    try:
        import ffsim
        print("ffsim detected - including ffsim_integration package")
        # Include ffsim_integration packages as well
        ffsim_packages = [pkg for pkg in all_packages if pkg.startswith('ffsim_integration')]
        return base_packages + ffsim_packages
    except ImportError:
        print("ffsim not available - excluding ffsim_integration package")
        # Only include base packages
        return base_packages


def get_package_data():
    """Get package data, excluding ffsim_integration if not available."""
    try:
        import ffsim
        # Include all package data
        return {}
    except ImportError:
        # No special package data needed when ffsim_integration is excluded
        return {}


if __name__ == "__main__":
    # Get packages based on ffsim availability
    packages = get_packages()
    package_data = get_package_data()
    
    print(f"Packages to include: {packages}")
    
    # Use setuptools with conditional packages
    setup(
        packages=packages,
        package_dir={"": "src"},
        package_data=package_data,
    )