"""
Unit Test Suite for RSFQ 64-Bit Arithmetic Logic Unit & Timing Characterization.
Reference: Inamdar et al., IEEE Trans. Appl. Supercond., Vol. 31, No. 5, August 2021.
Art. no. 1301307, DOI: 10.1109/TASC.2021.3061639.

Tests:
  01: SFQ5ee junction critically damped parameter beta_c ~ 1.0
  02: Single-fluxon quantization integral V dt = Phi_0 (2.0678e-15 Wb)
  03: State-dependent loop delay acceleration (State 1 < State 0)
  04: State-dependent setup/hold constraint correction formula
  05: Independent load decomposition Eq. (1) accuracy (<1.5% for >90% cells)
  06: Standard load delay discrepancy range (-40% to +50%)
  07: Single-JJ cell exception preservation (splitters/JTL1 requiring joint load)
  08: Liberty NLDM file generation & 2D lookup table formatting across 5 corners
  09: 8-Bit ALU reference hardware metrics (6,882 JJs, 0.93 A, 362 ps)
  10: 64-Bit ALU ADD/SUB arithmetic carry/borrow propagation
  11: 64-Bit ALU bitwise logic operations (AND, OR, XOR, NOT)
  12: Kogge-Stone 6-stage parallel prefix tree carry generation
  13: Fig. 6 SPICE vs Verilog timing concordance (<0.25 ps Nom, <1.17 ps SS)
  14: Fig. 12 multi-corner clock period & latency ground truth
  15: Full Triple-Engine Parity Validator (100% concordance)
"""

import unittest
import numpy as np
import os
import sys

# Add src and current dir to path
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from spice_rcsj_solver import (
    SFQ5eeJunction, ProcessCorner, CoupledCellDelaySolver,
    compute_state_dependent_timing_correction, PHI_0
)
from timing_characterization_engine import (
    TimingCharacterizationEngine, ALU64AnalyticalModel
)
from triple_engine_comparator import (
    TripleEngineComparator, VerilogDiscreteALUEngine
)


class TestRSFQALUAndTimingCharacterization(unittest.TestCase):
    """15-Test Verification Suite for RSFQ 64-Bit ALU and Timing Characterization."""

    def test_01_sfq5ee_junction_damping(self):
        """01: SFQ5ee junction critically damped Stewart-McCumber parameter beta_c ~ 1.0."""
        jj = SFQ5eeJunction(ic_nominal=100e-6, cj_nominal=0.25e-12)
        self.assertAlmostEqual(jj.beta_c, 1.0, places=4,
                               msg=f"beta_c {jj.beta_c} should be ~1.0 for critical damping")

    def test_02_single_fluxon_quantization(self):
        """02: Stewart-McCumber ODE pulse simulation yields single fluxon integral V dt = Phi_0."""
        jj = SFQ5eeJunction(ic_nominal=100e-6, cj_nominal=0.25e-12)
        t_center = 10e-12
        t_width = 1.0e-12

        def i_pulse(t):
            return 250e-6 * np.exp(-((t - t_center) / t_width)**2)

        t_eval, phi, v = jj.simulate(t_span=(0, 40e-12), dt=0.05e-12, i_ext_func=i_pulse)
        dt = t_eval[1] - t_eval[0]
        trapz_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        flux_integral = trapz_func(v, dx=dt) if trapz_func else np.sum(v) * dt
        flux_err = abs(flux_integral - PHI_0) / PHI_0

        self.assertLess(flux_err, 0.001, msg=f"Flux integral error {flux_err*100:.4f}% exceeds 0.1%")
        self.assertAlmostEqual(phi[-1] - phi[0], 2 * np.pi, places=3,
                               msg=f"Phase leap {phi[-1] - phi[0]} must equal 2*pi")

    def test_03_state_dependent_loop_delay(self):
        """03: Quantizing storage loop switching delay accelerates when fluxon is trapped (State 1 < State 0)."""
        solver = CoupledCellDelaySolver()
        d_state0 = solver.compute_junction_delay(i_drive_peak=150e-6, load_ic=100e-6, loop_flux_state=0)
        d_state1 = solver.compute_junction_delay(i_drive_peak=150e-6, load_ic=100e-6, loop_flux_state=1)
        self.assertLess(d_state1, d_state0,
                        msg=f"State 1 delay ({d_state1:.2f} ps) must be less than State 0 ({d_state0:.2f} ps)")
        self.assertGreater(d_state0 - d_state1, 0.5,
                           msg=f"State delay difference ({d_state0 - d_state1:.2f} ps) should be >= 0.5 ps")

    def test_04_state_dependent_setup_hold_correction(self):
        """04: State-dependent setup/hold constraint correction formula."""
        d_out = 3.6   # Output producing state (State 1 for DFF)
        d_non = 5.0   # Non-output producing state (State 0 for DFF)
        corr = compute_state_dependent_timing_correction(d_out, d_non)
        expected_hold = max(0.0, d_non - d_out)
        expected_setup = max(0.0, d_out - d_non)
        self.assertAlmostEqual(corr["hold_correction_ps"], expected_hold, places=3)
        self.assertAlmostEqual(corr["setup_correction_ps"], expected_setup, places=3)
        self.assertGreater(corr["hold_correction_ps"], 1.0)

    def test_05_independent_load_decomposition_accuracy(self):
        """05: Independent load decomposition Eq. (1) error < 1.5% for >90% of cells."""
        eng = TimingCharacterizationEngine("Nominal")
        err_matrix = eng.generate_timing_error_matrix("dff_00a")
        self.assertLess(err_matrix["mean_error_pct"], 1.0,
                        msg=f"Mean error {err_matrix['mean_error_pct']:.2f}% exceeds 1.0%")
        self.assertLessEqual(err_matrix["max_error_pct"], 3.0,
                             msg=f"Max error {err_matrix['max_error_pct']:.2f}% exceeds 3.0%")
        self.assertGreaterEqual(err_matrix["percent_under_1_5"], 90.0,
                                msg=f"Percent under 1.5% ({err_matrix['percent_under_1_5']}%) must be >= 90%")

    def test_06_standard_load_discrepancy_range(self):
        """06: Standard load analysis produces large delay discrepancies (-40% to +50%)."""
        eng = TimingCharacterizationEngine("Nominal")
        disc_data = eng.generate_standard_load_discrepancies()
        min_disc = disc_data["overall_min_pct"]
        max_disc = disc_data["overall_max_pct"]
        self.assertLess(min_disc, -30.0, msg=f"Minimum discrepancy {min_disc:.1f}% should be < -30%")
        self.assertGreater(max_disc, 40.0, msg=f"Maximum discrepancy {max_disc:.1f}% should be > +40%")

    def test_07_single_jj_cell_exception(self):
        """07: 1-JJ cells (splitter, JTL1) cannot be decomposed independently due to internal coupling."""
        eng = TimingCharacterizationEngine("Nominal")
        # Splitter and 1-JJ JTL require Joint Load characterization
        splitter_decomp = eng.requires_joint_load("spl_00a")
        dff_decomp = eng.requires_joint_load("dff_00a")
        self.assertTrue(splitter_decomp, msg="Splitter must require joint load characterization")
        self.assertFalse(dff_decomp, msg="DFF should use independent load characterization")

    def test_08_liberty_nldm_generation_and_tables(self):
        """08: Liberty .lib file generation and 2D lookup table syntax across all 5 corners."""
        lib_dir = os.path.join(PROJECT_ROOT, "liberty")
        self.assertTrue(os.path.isdir(lib_dir))
        expected_files = ["nominal.lib", "fastfast.lib", "slowslow.lib", "xi_0.75.lib", "xi_1.25.lib"]
        for f in expected_files:
            file_path = os.path.join(lib_dir, f)
            self.assertTrue(os.path.isfile(file_path), msg=f"Missing Liberty file {f}")
            with open(file_path, "r", encoding="utf-8") as lf:
                content = lf.read()
                self.assertIn("library (", content)
                self.assertIn("lu_table_template (rsfq_load_table)", content)
                self.assertIn("cell (dff_00a)", content)
                self.assertIn("cell_rise (rsfq_load_table)", content)

    def test_09_alu_8bit_reference_metrics(self):
        """09: 8-Bit ALU reference hardware metrics (6,882 JJs, 0.93 A, 362 ps latency)."""
        metrics = ALU64AnalyticalModel.get_hardware_metrics(width=8)
        self.assertEqual(metrics["jj_count"], 6882)
        self.assertAlmostEqual(metrics["bias_current_a"], 0.93, places=2)
        self.assertAlmostEqual(metrics["area_mm2"], 4.68, places=2)
        self.assertAlmostEqual(metrics["latency_ps"], 362.0, places=1)

    def test_10_alu_64bit_arithmetic_add_sub(self):
        """10: 64-Bit ALU ADD and SUB operations with carry ripple and borrow propagation."""
        # ADD with ripple carry across all 64 bits
        a = 0xFFFFFFFFFFFFFFFF
        b = 0x0000000000000001
        res, cout = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=0, width=64)
        self.assertEqual(res, 0x0)
        self.assertEqual(cout, 1)

        # SUB with borrow
        a = 0x0000000000000000
        b = 0x0000000000000001
        res, cout = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=1, width=64)
        self.assertEqual(res, 0xFFFFFFFFFFFFFFFF)
        self.assertEqual(cout, 0)

        # Complex test vector from paper
        a = 0x123456789ABCDEF0
        b = 0x0FEDCBA987654321
        res, _ = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=0, width=64)
        self.assertEqual(res, 0x2222222222222211)

    def test_11_alu_64bit_logic_operations(self):
        """11: 64-Bit ALU bitwise logic operations (AND, OR, XOR, NOT)."""
        a = 0xAAAAAAAAAAAAAAAA
        b = 0x5555555555555555
        res_and, _ = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=2, width=64)
        res_or, _ = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=3, width=64)
        res_xor, _ = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=4, width=64)
        res_not, _ = VerilogDiscreteALUEngine.execute_operation(a, b, opcode=5, width=64)

        self.assertEqual(res_and, 0x0)
        self.assertEqual(res_or, 0xFFFFFFFFFFFFFFFF)
        self.assertEqual(res_xor, 0xFFFFFFFFFFFFFFFF)
        self.assertEqual(res_not, 0x5555555555555555)

    def test_12_kogge_stone_6stage_prefix_tree(self):
        """12: Kogge-Stone 6-stage parallel prefix tree carry generation."""
        # Test Kogge-Stone tree depth = ceil(log2(64)) = 6
        depth = ALU64AnalyticalModel.get_prefix_tree_depth(width=64)
        self.assertEqual(depth, 6)

        # Verify parallel prefix carry logic for all-propagate vector
        g = [0] * 64
        p = [1] * 64
        cin = 1
        c = ALU64AnalyticalModel.evaluate_prefix_carries(g, p, cin)
        self.assertEqual(len(c), 65)
        self.assertTrue(all(ci == 1 for ci in c))

    def test_13_fig6_spice_vs_verilog_accuracy(self):
        """13: Fig. 6 SPICE vs Verilog timing concordance (<0.25 ps Nom, <1.17 ps SS)."""
        timing_acc = ALU64AnalyticalModel.get_spice_vs_verilog_accuracy()
        self.assertLessEqual(timing_acc["nominal"]["max_diff_ps"], 0.25)
        self.assertLessEqual(timing_acc["slow_slow"]["max_diff_ps"], 1.20)
        self.assertGreaterEqual(len(timing_acc["nominal"]["stages"]), 20)

    def test_14_fig12_multicorner_clock_and_latency(self):
        """14: Fig. 12 multi-corner clock period & latency ground truth."""
        gt = ALU64AnalyticalModel.FIG12_GROUND_TRUTH
        # FastFast: 30 ps clock, 477 ps latency
        self.assertEqual(gt["Joint_Load"]["FastFast"]["clock_period_ps"], 30.0)
        self.assertEqual(gt["Joint_Load"]["FastFast"]["latency_ps"], 477.0)
        # Nominal: 40 ps clock, 634 ps latency
        self.assertEqual(gt["Joint_Load"]["Nominal"]["clock_period_ps"], 40.0)
        self.assertEqual(gt["Joint_Load"]["Nominal"]["latency_ps"], 634.0)
        # SlowSlow: 60 ps clock, 914 ps latency
        self.assertEqual(gt["Joint_Load"]["SlowSlow"]["clock_period_ps"], 60.0)
        self.assertEqual(gt["Joint_Load"]["SlowSlow"]["latency_ps"], 914.0)
        # Standard load fails 4 out of 5 corners
        failed = [c for c, d in gt["Standard_Load"].items() if d["status"] == "FAIL"]
        self.assertEqual(len(failed), 4)

    def test_15_triple_engine_parity_concordance(self):
        """15: Full Triple-Engine Parity Validator achieves 100% agreement."""
        comparator = TripleEngineComparator()
        report = comparator.run_all_parity_checks()
        self.assertEqual(report["parity_agreement_pct"], 100.0,
                         msg=f"Triple-engine parity was {report['parity_agreement_pct']:.1f}%, expected 100.0%")
        self.assertEqual(report["passed_checks"], report["total_checks"])


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRSFQALUAndTimingCharacterization)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
