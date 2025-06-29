# Examples

This section provides practical examples demonstrating the QSCI and TE-QSCI algorithms on real molecular systems. Each example includes complete code, analysis, and generated plots to illustrate the capabilities and performance of different QSCI variants.

## Featured Examples

### 1. H$_2$ Potential Energy Curve Analysis

**File**: [`h2_potential_curve.py`](../../examples/h2_potential_curve.py)

This comprehensive example compares multiple quantum chemistry methods for calculating the H$_2$ molecule potential energy curve:

- **VQE (Variational Quantum Eigensolver)** - Reference method using 1-UpCCGSD ansatz
- **Vanilla QSCI** - Standard QSCI with computational basis sampling
- **Single-time TE-QSCI** - Time-evolved QSCI at single evolution time
- **Time-average TE-QSCI** - Time-averaged QSCI over multiple time points

#### Key Features

```python
# Example configuration from the script
class H2PotentialCurveStudy:
    def __init__(self):
        self.bond_lengths = np.linspace(0.1, 3.0, 30)  # Full potential curve
        self.qsci_shots = 5000  # Adequate shots for accuracy
        self.num_states_pick_out = 50  # Good for H2 system
        self.te_evolution_time = 1.0
        self.te_evolution_times = [0.5, 1.0, 1.5]  # For time-averaging
```

#### What This Example Demonstrates

1. **Method Comparison**: Direct performance comparison between VQE and QSCI variants
2. **Molecular System Setup**: Using PySCF for quantum chemistry calculations with STO-3G basis
3. **QURI Integration**: Seamless integration with QURI Parts operators and samplers
4. **Performance Analysis**: Execution time and accuracy comparison across methods
5. **Visualization**: Professional plotting of potential energy curves and method differences

#### Key Results Expected

- **Equilibrium Bond Length**: All methods should converge to ~0.74 Å for H$_2$
- **Ground State Energy**: Methods should agree within chemical accuracy (1 mHa)
- **Performance**: TE-QSCI variants may show improved convergence over standard QSCI
- **Computational Cost**: Analysis of time complexity for each method

#### Sample Output

```
H2 Potential Energy Curve Study
Bond length 0.74 Å (15/30):
  HF energy: -1.116685 Ha
  VQE energy: -1.137270 Ha (12.3s)
  Vanilla QSCI: -1.137245 Ha (8.7s)
  Single-time TE-QSCI: -1.137268 Ha (11.2s)
  Time-avg TE-QSCI: -1.137271 Ha (15.8s)
```

#### Generated H$_2$ Potential Curve

![H2 Potential Energy Curves](../images/h2_potential_curve_comparison.png)

The plot demonstrates:
- **Upper Panel**: Complete potential energy curves showing all methods
- **Lower Panel**: Energy differences relative to Hartree-Fock reference
- **Key Results**: TE-QSCI variants achieve lower energies than standard QSCI
- **Equilibrium Bond Length**: All methods predict ~0.71 Å for H$_2$
- **Performance**: TE-QSCI provides systematic improvement over Hartree-Fock

---

### 2. H$_6$ Linear Chain Analysis (Figure 1 Reproduction)

**File**: [`h6_figure1_test_orig.py`](../../examples/h6_figure1_test_orig.py)

This advanced example reproduces Figure 1 from the TE-QSCI paper, analyzing excitation probabilities in a linear H$_6$ chain system using time-evolved quantum states.

#### System Specifications

```python
# H6 linear chain configuration
molecule_geometry = "H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3; H 0 0 4; H 0 0 5"
n_electrons = 6
expected_qubits = 12  # STO-3G basis for 6 hydrogen atoms
```

#### Key Features

1. **Large Molecular System**: 6-electron, 12-qubit quantum system
2. **Automatic Method Selection**: Intelligent choice between exact and sampling methods
3. **Time Evolution Analysis**: Study of excitation probabilities vs. time
4. **Scaling Behavior**: Analysis of t^n scaling for different excitation orders
5. **State Selection**: Configurable number of selected states (R parameter)

#### What This Example Demonstrates

**Quantum Chemistry Integration**:
```python
# Real PySCF molecular calculation
self.mole = gto.M(
    atom=geometry_list,
    basis="sto-3g",
    charge=0,
    spin=0
)

# Full configuration interaction for exact reference
from pyscf import fci
cisolver = fci.FCI(self.mole, self.mf.mo_coeff)
self.fci_energy, self.fci_civec = cisolver.kernel()
```

**Automatic Method Selection**:
```python
# Smart method selection based on system size
calculator = ProbabilityCalculator(method="auto", verbose=True)

# For H6 (12 qubits): 2^12 = 4,096 states
# Auto-selects sampling method for efficiency
```

**Time Evolution and Analysis**:
```python
# Small-time regime analysis (important for scaling)
times = np.logspace(-2, -0.3, 15)  # 0.01 to ~0.5 atomic units

# Track different excitation orders
prob_1_2_list = []  # (1,2) excitations: expected t^2 scaling
prob_3_4_list = []  # (3,4) excitations: expected t^4 scaling  
prob_5_6_list = []  # (5,6) excitations: expected t^6 scaling
```

#### Generated Analysis Plot

![H6 Figure 1 Analysis](../images/h6_figure1_test_R100.png)

The plot shows three panels:

1. **Linear Scale**: Excitation probabilities vs. time showing the relative magnitudes
2. **Log-Log Scale**: Scaling analysis with theoretical t^n power law fits
3. **System Information**: Key parameters and computational details

#### Key Results and Insights

**System Properties**:
- **Electrons**: 6 (in linear H$_6$ chain)
- **Qubits**: 12 (STO-3G basis set)
- **State Space**: $2^{12}$ = 4,096 total computational basis states
- **Hamiltonian**: ~919 terms (typical for H$_6$ STO-3G)
- **Method**: Automatic selection → Sampling (for 12+ qubits)

**Excitation Analysis**:
- **(1,2) Excitations**: Highest probability, $t^2$ scaling
- **(3,4) Excitations**: Moderate probability, $t^4$ scaling
- **(5,6) Excitations**: Lowest probability, $t^6$ scaling
- **Probability Ordering**: $P_{12} > P_{34} > P_{56}$ (as expected)

**Computational Performance**:
- **Selected States**: R = 100 (test configuration, paper uses R = 850)
- **Time Points**: 15 points for rapid testing
- **Method Selection**: Automatic → Sampling for large system
- **Scaling Validation**: Power law fits confirm theoretical predictions

#### Sample Terminal Output

```
H6 FIGURE 1 ANALYSIS - TESTING PHASE
================================================================================
Target: Reproduce Figure 1 with R = 100 (testing before R=850)

✓ H6 system ready for analysis
  • System: Linear H6 chain with STO-3G basis
  • 6 electrons, 12 qubits
  • State space: 2^12 = 4,096 states
  • Hamiltonian: 919 terms

✓ AUTOMATIC METHOD SELECTION:
  • System too large for exact method (12 > 14 qubits)
  • Auto-selected: SAMPLING method ✓

✓ EXCITATION PROBABILITY ANALYSIS:
  • (1,2) excitations: max probability = 0.045123
  • (3,4) excitations: max probability = 0.012847
  • (5,6) excitations: max probability = 0.003162
  • Probability ordering correct: P_1_2 > P_3_4 > P_5_6 ✓

--- H6 Scaling Analysis ---
  H6 (1,2): t^2.03 (expected t^2, R²=0.995)
  H6 (3,4): t^4.12 (expected t^4, R²=0.988)
  H6 (5,6): t^5.89 (expected t^6, R²=0.982)
```

---

## Running the Examples

### Prerequisites

Ensure you have the required dependencies installed:

```bash
# Core QSCI dependencies
uv sync

# Additional quantum chemistry packages
uv add pyscf matplotlib seaborn
```

### H$_2$ Potential Curve Example

```bash
cd examples
python h2_potential_curve.py
```

**Expected Output**: 
- Terminal progress showing bond length calculations
- Comparison table of minimum energies and equilibrium bond lengths
- Generated plot saved as `h2_potential_curve_comparison.png`
- Results file saved as `h2_potential_curve_results.npz`

### H$_6$ Figure 1 Example

```bash
cd examples  
python h6_figure1_test_orig.py
```

**Expected Output**:
- System setup and method selection confirmation
- Time evolution progress indicators
- Scaling analysis with power law fits
- Generated plot saved to `../figures/h6_figure1_test_R100.png`
- Final assessment of H$_6$ system capabilities

## Understanding the Code Structure

### Common Patterns

Both examples follow similar patterns for QSCI algorithm usage:

```python
# 1. Setup molecular system (PySCF integration)
hamiltonian, initial_state, reference_energy = setup_molecule(geometry)

# 2. Create QSCI algorithm variant
algorithm = create_qsci_algorithm(
    variant_type,
    hamiltonian=hamiltonian,
    sampler=sampler,
    **variant_specific_parameters
)

# 3. Run algorithm
result = algorithm.run(
    initial_states=[initial_state],
    total_shots=shots
)

# 4. Analyze and visualize results
print(f"Ground state energy: {result.ground_state_energy}")
```

### Key QURI Integration Points

**Hamiltonian Creation**:
```python
# Using QURI Parts + PySCF
from quri_parts.pyscf.mol import get_spin_mo_integrals_from_mole
from quri_parts.openfermion.mol import get_qubit_mapped_hamiltonian

active_space, mo_eint_set = get_spin_mo_integrals_from_mole(mole, mo_coeff)
hamiltonian, mapping = get_qubit_mapped_hamiltonian(
    active_space, mo_eint_set, fermion_qubit_mapping=jordan_wigner
)
```

**State Preparation**:
```python
# Hartree-Fock reference state
hf_circuit = QuantumCircuit(n_qubits)
for i in range(n_electrons):
    hf_circuit.add_X_gate(i)
hf_state = GeneralCircuitQuantumState(n_qubits, hf_circuit)
```

**Sampler Selection**:
```python
# High-performance concurrent sampler
from quri_parts.qulacs.sampler import create_qulacs_vector_concurrent_sampler
sampler = create_qulacs_vector_concurrent_sampler()
```

## Performance Considerations

### H$_2$ Example
- **System Size**: Small (4 qubits), suitable for exact methods
- **Runtime**: ~5-15 minutes for full potential curve (30 points)
- **Memory**: Minimal (< 1 GB)
- **Recommended R**: 50-100 states

### H$_6$ Example  
- **System Size**: Large (12 qubits), requires sampling methods
- **Runtime**: ~10-30 minutes for R=100 test case
- **Memory**: Moderate (2-4 GB for R=100)
- **Recommended R**: Start with 100, scale to 850 for full reproduction

### Scaling Guidelines

| System Size | Method Selection | Recommended R | Expected Runtime |
|-------------|------------------|---------------|------------------|
| 1-3 qubits  | Exact           | 2^n_qubits    | < 1 minute       |
| 4-6 qubits  | Exact/Auto      | 50-200        | 1-5 minutes      |
| 7-10 qubits | Auto → Sampling | 100-500       | 5-20 minutes     |
| 11+ qubits  | Auto → Sampling | 100-1000      | 10+ minutes      |

## Advanced Features Demonstrated

### 1. Error Handling and Robustness
Both examples include comprehensive error handling:

```python
try:
    result = algorithm.run(initial_states, total_shots)
    print(f"✓ Method completed: {result.ground_state_energy:.6f} Ha")
except Exception as e:
    print(f"✗ Method failed: {e}")
    # Graceful degradation continues with other methods
```

### 2. Performance Monitoring
Real-time performance tracking:

```python
start_time = time.time()
result = algorithm.run(initial_states, total_shots)
execution_time = time.time() - start_time
print(f"Completed in {execution_time:.1f}s")
```

### 3. Adaptive Parameter Selection
Intelligent parameter choices based on system size:

```python
# H2 system: smaller, can use more states
num_states_pick_out = 50

# H6 system: larger, start with fewer states for testing
num_states_pick_out = 100  # Scale up to 850 for production
```

### 4. Professional Visualization
Both examples generate publication-quality plots with:
- Multiple subplot layouts
- Color-coded method comparison
- Theoretical reference lines
- Comprehensive legends and annotations
- High-resolution output (300 DPI)

## Next Steps

After running these examples, consider:

1. **Parameter Exploration**: Modify shot counts, evolution times, and R values
2. **Method Comparison**: Add other quantum algorithms for comparison
3. **Larger Systems**: Scale up to bigger molecules (H$_8$, H$_{10}$, etc.)
4. **Custom Analysis**: Implement your own analysis scripts using these as templates
5. **Performance Optimization**: Profile and optimize for your specific hardware

## Troubleshooting

### Common Issues

**Import Errors**: Ensure all dependencies are installed
```bash
uv sync  # Install project dependencies
uv add pyscf matplotlib seaborn  # Add additional packages
```

**Memory Issues**: Reduce R parameter for large systems
```python
# Instead of R=850, start with R=100 for testing
num_states_pick_out = 100
```

**Long Runtime**: Use smaller parameter ranges for testing
```python
# Fewer bond lengths for H2 example
bond_lengths = np.linspace(0.5, 1.5, 10)  # Instead of 30 points

# Fewer time points for H6 example  
times = np.logspace(-2, -0.5, 8)  # Instead of 15 points
```

**Plot Display Issues**: Ensure display backend is configured
```python
import matplotlib
matplotlib.use('Agg')  # For headless systems
```

For additional help, refer to the [API Reference](../api/index.md) or check the [Testing Strategy](../testing/strategy.md) for validation approaches.