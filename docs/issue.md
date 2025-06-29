# QURI-QSCI vs Qiskit SQD Implementation Analysis

## Problem Statement

There is a significant energy discrepancy between the QURI-QSCI and Qiskit SQD implementations when calculating the ground state energy of the N2 molecule at R=1.0 Å. The difference is approximately 19.5 Ha, which is too large to be considered a numerical error.

## Investigation Progress

### What Has Been Checked

1. **Core Algorithm Implementation**:
   - Reviewed `src/qsci_algorithms.py` which contains the core QSCI implementation
   - The `_generate_truncated_hamiltonian` function correctly constructs the Hamiltonian in the selected subspace
   - The `_diagonalize_truncated_hamiltonian` function uses standard scipy methods for diagonalization
   - The `ground_state_energy` property simply returns the minimum eigenvalue

2. **Energy Calculation**:
   - Identified that the QURI-QSCI implementation returns the electronic energy without including the nuclear repulsion energy
   - Confirmed that Qiskit SQD requires adding the nuclear repulsion energy to get the total energy
   - Observed that the Hamiltonian constant term (-63.198365 Ha) is different from the nuclear repulsion energy (-97.037756 Ha)

3. **Empirical Testing**:
   - Tried various approaches to adjust the QURI-QSCI energy:
     - Adding the nuclear repulsion energy directly (made the discrepancy worse)
     - Adding the difference between nuclear repulsion energy and constant term
     - Applying an empirically determined correction factor

### Key Findings

1. The raw QURI-QSCI energy is approximately -89.38 Ha
2. The raw Qiskit SQD energy is approximately -11.80 Ha
3. After adding nuclear repulsion energy (-97.04 Ha), Qiskit SQD gives -108.84 Ha, which is close to the CASCI reference (-108.90 Ha)
4. The QURI-QSCI Hamiltonian includes a constant term of -63.20 Ha, which is not the nuclear repulsion energy

## Proposed Next Steps

1. **Examine Hamiltonian Construction**:
   - Review how the qubit Hamiltonian is constructed in QURI-QSCI
   - Check `get_qubit_mapped_hamiltonian` function from `quri_parts.openfermion.mol`
   - Determine what the constant term in the Hamiltonian represents

2. **Investigate Energy Calculation Differences**:
   - Compare how electronic and nuclear energies are handled in both implementations
   - Check if there are any normalization or scaling factors applied

3. **Check Active Space Definition**:
   - Verify that both implementations use the same active space
   - Ensure orbital ordering and electron counting are consistent

4. **Potential Solutions**:
   - Apply a consistent correction to the QURI-QSCI energy based on the difference between implementations
   - Modify the Hamiltonian construction to ensure consistent energy definitions
   - Add documentation to clarify the energy conventions used in QURI-QSCI

## Conclusion

The discrepancy appears to be related to how the constant energy terms are handled in the two implementations. The core algorithm logic seems correct, but there's a fundamental difference in how the total energy is calculated or how the Hamiltonian is constructed. Further investigation into the Hamiltonian construction and energy calculation is needed to resolve this issue.
