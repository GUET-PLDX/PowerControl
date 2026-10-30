"""Build and exercise the production RLS implementation on the host."""

import pathlib
import subprocess
import tempfile
import unittest


class RlsNumericalSafetyTest(unittest.TestCase):
    def test_covariance_safety(self):
        module = pathlib.Path(__file__).resolve().parents[1]
        workspace = next(
            parent for parent in module.parents
            if (parent / "Middlewares/Third_Party/LibXR/lib/Eigen/Eigen/Core").is_file()
        )
        eigen = workspace / "Middlewares/Third_Party/LibXR/lib/Eigen"
        with tempfile.TemporaryDirectory() as directory:
            binary = pathlib.Path(directory) / "rls_numerical_safety_test"
            subprocess.run(
                ["c++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-Werror",
                 "-ffp-contract=off", "-I", str(module), "-I", str(eigen),
                 str(module / "tests/rls_numerical_safety_test.cpp"), "-o", str(binary)],
                check=True,
            )
            subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    unittest.main()
