"""Black-box failure contracts for the Paper I build driver (no scientific dependencies).

Run from NonStationaryMetrics/: python3 tests/test_build_paper1.py
Generators and latexmk are replaced in temporary directories; no real paper is modified.
"""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


DRIVER = Path(__file__).resolve().parents[1] / "paper1/build_paper1.sh"
FIGURES = (
    "fig_indicatrices", "fig_flrw_kinematics", "fig_flrw_variational",
    "fig_flrw_curvature", "fig_kodama_conservation", "fig_vaidya_trajectories",
    "fig_vaidya_plunge_t_tau", "fig_vaidya_no_inversion_evaporation",
    "fig_vaidya_bounce", "fig_vaidya_timing", "fig_vaidya_offshell",
    "fig_vaidya_kerr_a0", "fig_brachistochrone_minimum_check",
)


class BuildDriverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.paper = self.root / "paper1"
        self.paper.mkdir()
        source = DRIVER.read_text()
        (self.paper / DRIVER.name).write_text(source)
        self.generators = []
        for name in re.findall(r'^\s+"([^"\n]+\.py)"', source, re.M):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("raise SystemExit(0)\n")
            self.generators.append(path)
        self.assertTrue(self.generators, "fixture must exercise the listed generators")
        self.assets = self.root / "PaperFigures"
        self.assets.mkdir()
        for name in FIGURES:
            (self.assets / (name + ".pdf")).write_bytes(b"figure fixture\n")
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.executable("latexmk", """#!/bin/sh
printf invoked > latexmk-called
if [ "${PAPER1_TEST_LATEX_FAIL:-0}" = 1 ]; then exit 7; fi
if [ "${PAPER1_TEST_NO_PDF:-0}" != 1 ]; then
  printf '%s\n' 'new PDF fixture' > paper1_JMP.pdf
fi
""")
        self.env = dict(os.environ)
        self.env["PATH"] = os.pathsep.join((str(self.bin), str(Path(sys.executable).parent), self.env["PATH"]))

    def executable(self, name, source):
        path = self.bin / name
        path.write_text(source)
        path.chmod(0o755)

    def run_driver(self, *args):
        return subprocess.run(
            ["bash", DRIVER.name, *args], cwd=self.paper, env=self.env,
            capture_output=True, text=True, timeout=30,
        )

    def assert_stopped_before_latex(self, result):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((self.paper / "latexmk-called").exists())

    def test_full_build_success(self):
        result = self.run_driver()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("13/13 figures present", result.stdout)
        self.assertTrue((self.paper / "paper1_JMP.pdf").exists())
        self.assertTrue((self.paper / "latexmk-called").exists())

    def test_figures_mode_does_not_compile(self):
        result = self.run_driver("--figures")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("13/13 figures present", result.stdout)
        self.assertFalse((self.paper / "latexmk-called").exists())

    def test_missing_generator_fails_even_with_existing_figures(self):
        self.generators[0].unlink()
        self.assert_stopped_before_latex(self.run_driver("--figures"))

    def test_failed_generator_stops_full_build(self):
        self.generators[0].write_text("raise SystemExit(7)\n")
        self.assert_stopped_before_latex(self.run_driver())

    def test_missing_figure_stops_paper_mode(self):
        (self.assets / (FIGURES[0] + ".pdf")).unlink()
        self.assert_stopped_before_latex(self.run_driver("--paper"))

    def test_copy_failure_stops_paper_mode(self):
        self.executable("cp", "#!/bin/sh\nexit 7\n")
        self.assert_stopped_before_latex(self.run_driver("--paper"))

    def test_destination_creation_failure_stops_paper_mode(self):
        directory = self.root / "paper"
        directory.mkdir()
        (directory / "Immagini").write_text("not a directory\n")
        self.assert_stopped_before_latex(self.run_driver("--paper"))

    def test_latex_failure_preserves_old_pdf_and_fails(self):
        pdf = self.paper / "paper1_JMP.pdf"
        pdf.write_bytes(b"previous PDF fixture\n")
        self.env["PAPER1_TEST_LATEX_FAIL"] = "1"
        result = self.run_driver("--paper")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.paper / "latexmk-called").exists())
        self.assertNotIn("built paper1_JMP.pdf", result.stdout)
        self.assertEqual(pdf.read_bytes(), b"previous PDF fixture\n")

    def test_successful_latex_without_pdf_fails(self):
        self.env["PAPER1_TEST_NO_PDF"] = "1"
        result = self.run_driver("--paper")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.paper / "latexmk-called").exists())

    def test_paper_mode_success(self):
        result = self.run_driver("--paper")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.paper / "latexmk-called").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
