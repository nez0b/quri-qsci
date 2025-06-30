"""
Comprehensive tests for qsci_vm_analysis.py - Core VM Components

Tests VMFactory, VMSampler, and QURIVMInterface functionality including:
- VM creation with proper configurations
- Error rate validation and boundary conditions  
- Sampling delegation and error handling
- Circuit analysis dispatch based on lowering levels
- Default VM creation and initialization
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from quri_parts.backend.units import TimeValue, TimeUnit

# Module under test
from src.qsci_vm_analysis import (
    VMFactory,
    VMSampler,
    QURIVMInterface,
)
from src.qsci_algo_interface import LoweringLevel


@pytest.fixture
def mock_circuit():
    """Provides a generic mock QuantumCircuit."""
    circuit = Mock(name="QuantumCircuit")
    circuit.gates = []
    circuit.depth = 5
    circuit.qubit_count = 4
    return circuit


class TestVMFactory:
    """Tests for the VMFactory class."""

    @patch('src.qsci_vm_analysis.VM')
    def test_create_ideal_vm(self, mock_vm_class):
        """
        Verify that create_ideal_vm creates an abstract VM instance.
        This corresponds to calling VM() with no arguments.
        """
        # Arrange
        mock_vm_instance = Mock(name="IdealVM")
        mock_vm_class.return_value = mock_vm_instance

        # Act
        vm = VMFactory.create_ideal_vm()

        # Assert
        mock_vm_class.assert_called_once_with()
        assert vm is mock_vm_instance

    @patch('src.qsci_vm_analysis.star_device')
    @patch('src.qsci_vm_analysis.VM')
    @patch('src.qsci_vm_analysis.TimeValue')
    @patch('src.qsci_vm_analysis.TimeUnit')
    def test_create_star_vm_success(self, mock_time_unit, mock_time_value, mock_vm_class, mock_star_device):
        """
        Verify that create_star_vm correctly configures a STAR VM.
        It should generate a device property with the given error rate
        and create a VM from that property.
        """
        # Arrange
        error_rate = 1e-3
        mock_device_prop = Mock(name="MockDeviceProperty")
        mock_star_vm_instance = Mock(name="StarVM")
        mock_qec_cycle = Mock(name="QECCycle")

        mock_time_value.return_value = mock_qec_cycle
        mock_time_unit.MICROSECOND = "MICROSECOND"
        mock_star_device.generate_device_property.return_value = mock_device_prop
        mock_vm_class.from_device_prop.return_value = mock_star_vm_instance

        # Act
        vm = VMFactory.create_star_vm(error_rate=error_rate)

        # Assert
        # Check that TimeValue was created correctly
        mock_time_value.assert_called_once_with(value=1.0, unit="MICROSECOND")
        
        # Check that the device property was generated with the correct error rate
        mock_star_device.generate_device_property.assert_called_once()
        call_args, call_kwargs = mock_star_device.generate_device_property.call_args
        assert call_kwargs.get("physical_error_rate") == error_rate
        assert call_kwargs.get("qubit_count") == 16
        assert call_kwargs.get("code_distance") == 7
        assert call_kwargs.get("qec_cycle") == mock_qec_cycle

        # Check that the VM was created from this property
        mock_vm_class.from_device_prop.assert_called_once_with(mock_device_prop)
        assert vm is mock_star_vm_instance

    @patch('src.qsci_vm_analysis.star_device')
    @patch('src.qsci_vm_analysis.VM')
    @patch('src.qsci_vm_analysis.TimeValue')
    @patch('src.qsci_vm_analysis.TimeUnit')
    def test_create_star_vm_default_error_rate(self, mock_time_unit, mock_time_value, mock_vm_class, mock_star_device):
        """Test that create_star_vm works with default error rate of 0.0."""
        # Arrange
        mock_device_prop = Mock()
        mock_vm_instance = Mock()
        mock_qec_cycle = Mock()

        mock_time_value.return_value = mock_qec_cycle
        mock_time_unit.MICROSECOND = "MICROSECOND"
        mock_star_device.generate_device_property.return_value = mock_device_prop
        mock_vm_class.from_device_prop.return_value = mock_vm_instance

        # Act
        vm = VMFactory.create_star_vm()  # No error_rate specified

        # Assert
        call_args, call_kwargs = mock_star_device.generate_device_property.call_args
        assert call_kwargs.get("physical_error_rate") == 0.0
        assert vm is mock_vm_instance

    @patch('src.qsci_vm_analysis.VMFactory.create_star_vm')
    def test_create_star_vm_for_h6(self, mock_create_star):
        """Test that create_star_vm_for_h6 delegates to create_star_vm."""
        # Arrange
        error_rate = 1e-4
        mock_vm = Mock()
        mock_create_star.return_value = mock_vm

        # Act
        vm = VMFactory.create_star_vm_for_h6(error_rate)

        # Assert
        mock_create_star.assert_called_once_with(error_rate)
        assert vm is mock_vm


class TestVMSampler:
    """Tests for the VMSampler class."""

    def test_init(self):
        """Test VMSampler initialization."""
        # Arrange
        mock_vm = Mock(name="MockVM")

        # Act
        sampler = VMSampler(vm_instance=mock_vm)

        # Assert
        assert sampler.vm is mock_vm

    def test_call_success(self):
        """
        Verify that the sampler correctly delegates to the vm.sample method
        for each circuit-shot pair.
        """
        # Arrange
        mock_vm = Mock()
        mock_vm.sample.side_effect = [
            "result1",
            "result2",
        ]  # Each call returns a different value
        sampler = VMSampler(vm_instance=mock_vm)

        circuit1, shots1 = Mock(name="circuit1"), 100
        circuit2, shots2 = Mock(name="circuit2"), 200
        circuit_shot_pairs = [(circuit1, shots1), (circuit2, shots2)]

        # Act
        results = sampler(circuit_shot_pairs)

        # Assert
        assert mock_vm.sample.call_count == 2
        mock_vm.sample.assert_any_call(circuit1, shots1)
        mock_vm.sample.assert_any_call(circuit2, shots2)
        assert results == ["result1", "result2"]

    def test_call_vm_raises_error(self):
        """
        Verify that if the underlying vm.sample method fails, the exception
        is wrapped in a RuntimeError.
        """
        # Arrange
        mock_vm = Mock()
        original_exception = ValueError("Underlying VM error")
        mock_vm.sample.side_effect = original_exception
        sampler = VMSampler(vm_instance=mock_vm)
        circuit_shot_pairs = [(Mock(), 100)]

        # Act & Assert
        with pytest.raises(RuntimeError, match="VM sampling failed: Underlying VM error"):
            sampler(circuit_shot_pairs)

    def test_call_with_empty_list(self):
        """
        Verify that calling the sampler with an empty list of jobs
        returns an empty list without calling the VM.
        """
        # Arrange
        mock_vm = Mock()
        sampler = VMSampler(vm_instance=mock_vm)

        # Act
        results = sampler([])

        # Assert
        assert results == []
        mock_vm.sample.assert_not_called()

    def test_call_with_zero_shots(self):
        """Test behavior with zero shots."""
        # Arrange
        mock_vm = Mock()
        mock_vm.sample.return_value = {}
        sampler = VMSampler(vm_instance=mock_vm)
        circuit_shot_pairs = [(Mock(), 0)]

        # Act
        results = sampler(circuit_shot_pairs)

        # Assert
        mock_vm.sample.assert_called_once()
        assert len(results) == 1

    def test_call_with_large_shots(self):
        """Test behavior with large shot counts."""
        # Arrange
        mock_vm = Mock()
        mock_vm.sample.return_value = {"00": 500000, "11": 500000}
        sampler = VMSampler(vm_instance=mock_vm)
        circuit_shot_pairs = [(Mock(), 1000000)]

        # Act
        results = sampler(circuit_shot_pairs)

        # Assert
        mock_vm.sample.assert_called_once()
        assert len(results) == 1


class TestQURIVMInterface:
    """Tests for the QURIVMInterface class."""

    @patch("src.qsci_vm_analysis.VMFactory")
    def test_init_with_defaults(self, mock_vm_factory):
        """
        Verify that if no VMs are provided, the interface creates default
        ideal and STAR VMs via the factory.
        """
        # Arrange
        mock_ideal_vm = Mock(name="DefaultIdealVM")
        mock_star_vm = Mock(name="DefaultStarVM")
        mock_vm_factory.create_ideal_vm.return_value = mock_ideal_vm
        mock_vm_factory.create_star_vm.return_value = mock_star_vm

        # Act
        interface = QURIVMInterface()

        # Assert
        mock_vm_factory.create_ideal_vm.assert_called_once()
        mock_vm_factory.create_star_vm.assert_called_once_with(0.0)
        assert interface.logical_vm is mock_ideal_vm
        assert interface.arch_vm is mock_star_vm

    def test_init_with_provided_vms(self):
        """
        Verify that the interface uses the VM instances provided to it.
        """
        # Arrange
        logical_vm = Mock(name="ProvidedLogicalVM")
        arch_vm = Mock(name="ProvidedArchVM")

        # Act
        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)

        # Assert
        assert interface.logical_vm is logical_vm
        assert interface.arch_vm is arch_vm

    def test_init_partial_vms(self):
        """Test initialization with only one VM provided."""
        # Arrange
        provided_logical_vm = Mock(name="ProvidedLogical")
        
        with patch("src.qsci_vm_analysis.VMFactory") as mock_factory:
            mock_factory.create_star_vm.return_value = Mock(name="DefaultStar")
            
            # Act
            interface = QURIVMInterface(logical_vm=provided_logical_vm)
            
            # Assert
            assert interface.logical_vm is provided_logical_vm
            mock_factory.create_star_vm.assert_called_once_with(0.0)

    @pytest.mark.parametrize(
        "level, expected_vm_attr",
        [
            (LoweringLevel.LogicalCircuit, "logical_vm"),
            (LoweringLevel.ArchLogicalCircuit, "arch_vm"),
        ],
    )
    def test_analyze_circuit_at_level(self, level, expected_vm_attr, mock_circuit):
        """
        Verify that circuit analysis is dispatched to the correct VM
        based on the specified lowering level.
        """
        # Arrange
        logical_vm = Mock(name="LogicalVM")
        logical_vm.analyze.return_value = {"type": "logical_analysis"}
        arch_vm = Mock(name="ArchVM")
        arch_vm.analyze.return_value = {"type": "arch_analysis"}

        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)
        
        # Act
        result = interface.analyze_circuit_at_level(mock_circuit, level)

        # Assert
        if level == LoweringLevel.LogicalCircuit:
            logical_vm.analyze.assert_called_once_with(mock_circuit)
            arch_vm.analyze.assert_not_called()
            assert result == {"type": "logical_analysis"}
        else:
            arch_vm.analyze.assert_called_once_with(mock_circuit)
            logical_vm.analyze.assert_not_called()
            assert result == {"type": "arch_analysis"}

    def test_analyze_circuit_unknown_level(self, mock_circuit):
        """Test that unknown lowering levels default to arch analysis."""
        # Arrange
        logical_vm = Mock()
        arch_vm = Mock()
        arch_vm.analyze.return_value = {"type": "default_arch"}
        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)
        
        # Act - Use a mock level that doesn't match known levels
        unknown_level = Mock()
        result = interface.analyze_circuit_at_level(mock_circuit, unknown_level)
        
        # Assert
        arch_vm.analyze.assert_called_once_with(mock_circuit)
        logical_vm.analyze.assert_not_called()
        assert result == {"type": "default_arch"}

    @pytest.mark.parametrize(
        "level, expected_vm_attr",
        [
            (LoweringLevel.LogicalCircuit, "logical_vm"),
            (LoweringLevel.ArchLogicalCircuit, "arch_vm"),
        ],
    )
    def test_create_sampler(self, level, expected_vm_attr):
        """
        Verify that sampler creation uses the correct VM instance
        based on the specified lowering level.
        """
        # Arrange
        vms = {
            "logical_vm": Mock(name="LogicalVM"),
            "arch_vm": Mock(name="ArchVM"),
        }
        interface = QURIVMInterface(logical_vm=vms["logical_vm"], arch_vm=vms["arch_vm"])

        # Act
        sampler = interface.create_sampler(level)

        # Assert
        assert isinstance(sampler, VMSampler)
        assert sampler.vm is vms[expected_vm_attr]

    def test_create_logical_sampler(self):
        """Test direct creation of logical sampler."""
        # Arrange
        logical_vm = Mock(name="LogicalVM")
        arch_vm = Mock(name="ArchVM")
        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)

        # Act
        sampler = interface.create_logical_sampler()

        # Assert
        assert isinstance(sampler, VMSampler)
        assert sampler.vm is logical_vm

    def test_create_arch_sampler(self):
        """Test direct creation of arch sampler."""
        # Arrange
        logical_vm = Mock(name="LogicalVM")
        arch_vm = Mock(name="ArchVM")
        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)

        # Act
        sampler = interface.create_arch_sampler()

        # Assert
        assert isinstance(sampler, VMSampler)
        assert sampler.vm is arch_vm

    def test_vm_analyze_error_propagation(self, mock_circuit):
        """Test that VM analyze errors are properly propagated."""
        # Arrange
        logical_vm = Mock()
        logical_vm.analyze.side_effect = RuntimeError("VM analyze failed")
        arch_vm = Mock()
        interface = QURIVMInterface(logical_vm=logical_vm, arch_vm=arch_vm)

        # Act & Assert
        with pytest.raises(RuntimeError, match="VM analyze failed"):
            interface.analyze_circuit_at_level(mock_circuit, LoweringLevel.LogicalCircuit)


class TestCreateStarVMForH6:
    """Tests for the create_star_vm_for_h6 function."""

    @patch("src.qsci_vm_analysis.VMFactory")
    def test_create_star_vm_for_h6_function(self, mock_vm_factory):
        """Test the standalone create_star_vm_for_h6 function."""
        # Arrange
        from src.qsci_vm_analysis import create_star_vm_for_h6
        error_rate = 1e-3
        mock_vm = Mock()
        mock_interface = Mock()
        mock_vm_factory.create_ideal_vm.return_value = Mock()
        mock_vm_factory.create_star_vm_for_h6.return_value = mock_vm

        # Act
        interface = create_star_vm_for_h6(error_rate)

        # Assert
        assert isinstance(interface, QURIVMInterface)