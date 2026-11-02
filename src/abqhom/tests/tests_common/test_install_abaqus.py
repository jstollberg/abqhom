# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Tests for the isolated Abaqus runtime installer."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import abqhom.install_abaqus as installer
from abqhom import __version__


TEST_REQUIREMENTS = (
    "gmsh>=4.9.4",
    "paraqus>=1.0.2",
    "backports.weakref>=1.0.post1",
)


@unittest.skipIf(sys.version_info < (3, 9),
                 "the packaging installer requires Python 3.9 or newer")
class _InstallerTestCase(unittest.TestCase):
    """Provide isolated runtime fixtures for installer tests."""

    def setUp(self):
        """Create an isolated directory for each installer test."""
        self.directory = tempfile.mkdtemp(prefix="abqhom_installer_")

    def tearDown(self):
        """Remove the isolated installer directory."""
        shutil.rmtree(self.directory)

    def _create_runtime(self, prefix, distribution="site-packages",
                        include_weakref=True, include_library=True,
                        version=__version__):
        """Create the minimum layout expected after package installation."""
        package_directory = os.path.join(
            prefix, "Lib", "python", distribution)
        os.makedirs(os.path.join(package_directory, "abqhom"))
        os.makedirs(os.path.join(package_directory, "paraqus"))
        os.makedirs(os.path.join(package_directory, "backports"))
        abqhom_init = os.path.join(
            package_directory, "abqhom", "__init__.py")
        with open(abqhom_init, "w") as output:
            output.write('__version__ = "{}"\n'.format(version))
        for path in (
                os.path.join(package_directory, "paraqus", "__init__.py"),
                os.path.join(package_directory, "gmsh.py")):
            with open(path, "w"):
                pass
        if include_weakref:
            with open(os.path.join(
                    package_directory, "backports", "weakref.py"), "w"):
                pass
        if include_library:
            with open(os.path.join(prefix, "gmsh-test.dll"), "w"):
                pass
        return package_directory

    def _write_valid_marker(self, prefix, package_directory):
        """Write a valid ownership marker for a test runtime."""
        installer._write_marker(
            prefix, prefix, package_directory, TEST_REQUIREMENTS)


class TestRequirementsAndCommand(_InstallerTestCase):
    """Test dependency selection and pip command construction."""

    def test_default_prefix_is_versioned(self):
        """The default prefix separates abqhom and Python versions."""
        self.assertIn(__version__, installer._default_prefix())

    def test_runtime_requirements_follow_package_metadata(self):
        """Only Abaqus-facing dependencies enter the isolated runtime."""
        from importlib import metadata

        original_requires = metadata.requires
        metadata.requires = lambda name: [
            "numpy>=1.24",
            "paraqus>=1.0.2",
            "gmsh>=4.9.4",
            "backports.weakref>=1.0.post1",
        ]
        try:
            requirements = installer._runtime_requirements()
        finally:
            metadata.requires = original_requires

        self.assertEqual(requirements, TEST_REQUIREMENTS)

    def test_pip_command_installs_only_runtime_packages(self):
        """The pip command disables dependency installation into the prefix."""
        command = installer._pip_command(
            self.directory, "local.whl", TEST_REQUIREMENTS)
        self.assertIn("--ignore-installed", command)
        self.assertIn("--no-deps", command)
        self.assertIn("local.whl", command)
        self.assertNotIn("numpy", " ".join(command).lower())
        for requirement in TEST_REQUIREMENTS:
            self.assertIn(requirement, command)


class TestManagedRuntime(_InstallerTestCase):
    """Test managed-prefix ownership, layout, and runtime validation."""

    def test_nonempty_unmanaged_prefix_is_rejected(self):
        """An unrelated directory cannot be reused as a runtime prefix."""
        path = os.path.join(self.directory, "unmanaged")
        os.mkdir(path)
        with open(os.path.join(path, "user_file.txt"), "w"):
            pass
        self.assertRaises(RuntimeError, installer._prepare_prefix, path)

    def test_empty_marker_is_rejected(self):
        """An empty marker does not establish ownership of a prefix."""
        with open(os.path.join(self.directory, installer.MARKER_FILE), "w"):
            pass
        self.assertRaises(
            RuntimeError, installer._prepare_prefix, self.directory)

    def test_marker_for_another_prefix_is_rejected(self):
        """A marker cannot authorize replacement of another directory."""
        package_directory = self._create_runtime(self.directory)
        other = os.path.join(self.directory, "other-runtime")
        installer._write_marker(
            self.directory, other, package_directory, TEST_REQUIREMENTS)
        self.assertRaises(
            RuntimeError, installer._prepare_prefix, self.directory)

    def test_marker_package_directory_must_be_inside_prefix(self):
        """A managed package directory cannot escape its runtime prefix."""
        outside = os.path.join(os.path.dirname(self.directory), "outside")
        installer._write_marker(
            self.directory, self.directory, outside, TEST_REQUIREMENTS)
        self.assertRaises(
            RuntimeError, installer._prepare_prefix, self.directory)

    def test_marker_records_and_validates_runtime_configuration(self):
        """The ownership marker records a validated runtime configuration."""
        package_directory = self._create_runtime(self.directory)
        self._write_valid_marker(self.directory, package_directory)
        data = installer._read_marker(self.directory)
        self.assertEqual(data["abqhom_version"], __version__)
        self.assertEqual(data["managed_by"], "abqhom")
        self.assertEqual(data["prefix"], self.directory)
        self.assertEqual(data["site_packages"], package_directory)
        self.assertEqual(data["requirements"], list(TEST_REQUIREMENTS))

    def test_valid_managed_prefix_can_be_reused(self):
        """A valid managed prefix may be replaced transactionally."""
        package_directory = self._create_runtime(self.directory)
        self._write_valid_marker(self.directory, package_directory)
        self.assertEqual(
            installer._prepare_prefix(self.directory), self.directory)

    def test_site_packages_layout_is_found(self):
        """A conventional site-packages directory is discovered."""
        expected = self._create_runtime(self.directory)
        self.assertEqual(
            installer._find_site_packages(self.directory), expected)

    def test_dist_packages_layout_is_found(self):
        """A distribution-specific dist-packages directory is discovered."""
        expected = self._create_runtime(
            self.directory, distribution="dist-packages")
        self.assertEqual(
            installer._find_site_packages(self.directory), expected)

    def test_missing_weakref_is_rejected(self):
        """A runtime without the weak-reference backport is incomplete."""
        package_directory = self._create_runtime(
            self.directory, include_weakref=False)
        self.assertRaises(
            RuntimeError, installer._check_runtime_contents,
            self.directory, package_directory)

    def test_missing_gmsh_library_is_rejected(self):
        """A runtime without a native Gmsh library is incomplete."""
        package_directory = self._create_runtime(
            self.directory, include_library=False)
        self.assertRaises(
            RuntimeError, installer._check_runtime_contents,
            self.directory, package_directory)

    def test_external_numpy_directory_is_rejected(self):
        """The runtime must not shadow the NumPy bundled with Abaqus."""
        package_directory = self._create_runtime(self.directory)
        os.mkdir(os.path.join(package_directory, "numpy"))
        self.assertRaises(
            RuntimeError, installer._check_runtime_contents,
            self.directory, package_directory)

    def test_external_numpy_module_is_rejected(self):
        """A standalone NumPy module is also excluded from the runtime."""
        package_directory = self._create_runtime(self.directory)
        with open(os.path.join(package_directory, "numpy.py"), "w"):
            pass
        self.assertRaises(
            RuntimeError, installer._check_runtime_contents,
            self.directory, package_directory)

    def test_mismatched_abqhom_version_is_rejected(self):
        """The runtime package must match the running installer version."""
        package_directory = self._create_runtime(
            self.directory, version="9.9")
        self.assertRaises(
            RuntimeError, installer._check_runtime_imports,
            self.directory, package_directory)


class TestActivationAndCommandLine(_InstallerTestCase):
    """Test activation scripts and command-line entry points."""

    def test_activation_script_sets_pythonpath(self):
        """The generated shell activation changes PYTHONPATH when executed."""
        package_directory = self._create_runtime(self.directory)
        scripts = installer._write_activation_scripts(
            self.directory, package_directory)
        environment = os.environ.copy()
        environment["PYTHONPATH"] = ""
        if os.name == "nt":
            batch = [path for path in scripts if path.endswith(".bat")][0]
            command = os.path.join(self.directory, "check_activation.bat")
            with open(command, "w") as output_file:
                output_file.write(
                    '@call "{}"\n@echo %PYTHONPATH%\n'.format(batch))
            output = subprocess.check_output(
                ["cmd.exe", "/d", "/c", command],
                env=environment,
                universal_newlines=True)
            self.assertEqual(output.strip(), package_directory)

            powershell = [path for path in scripts
                          if path.endswith(".ps1")][0]
            quoted = powershell.replace("'", "''")
            command = ". '{}'; [Console]::Out.Write($env:PYTHONPATH)".format(
                quoted)
            powershell_executable = os.path.join(
                os.environ.get("SystemRoot", r"C:\Windows"),
                "System32", "WindowsPowerShell", "v1.0", "powershell.exe")
            output = subprocess.check_output(
                [powershell_executable, "-NoProfile", "-ExecutionPolicy",
                 "Bypass", "-Command", command],
                env=environment,
                universal_newlines=True)
        else:
            shell = [path for path in scripts if path.endswith(".sh")][0]
            output = subprocess.check_output(
                ["/bin/sh", "-c", '. "$1"; printf %s "$PYTHONPATH"',
                 "sh", shell],
                env=environment,
                universal_newlines=True)
        self.assertEqual(output.strip(), package_directory)

    def test_main_reports_installed_runtime(self):
        """The command-line entry point reports paths returned by install."""
        original_install = installer.install
        original_stdout = sys.stdout
        script = os.path.join(self.directory, "activate_abqhom.sh")
        output = StringIO()
        installer.install = lambda prefix, source: (self.directory, [script])
        sys.stdout = output
        try:
            result = installer.main(["--prefix", self.directory])
        finally:
            installer.install = original_install
            sys.stdout = original_stdout

        self.assertEqual(result, 0)
        self.assertIn(self.directory, output.getvalue())
        self.assertIn(script, output.getvalue())

    def test_console_bootstrap_ignores_shadowing_checkout(self):
        """The entry point imports abqhom from its own distribution."""
        source_root = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", ".."))
        shadow_root = os.path.join(self.directory, "shadow")
        shadow_package = os.path.join(shadow_root, "abqhom")
        os.makedirs(shadow_package)
        with open(os.path.join(shadow_package, "__init__.py"), "w") as output:
            output.write('__version__ = "0.0-shadow"\n')

        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join(
            (shadow_root, source_root))
        output = subprocess.check_output(
            [sys.executable, "-m", "abqhom_installer", "--help"],
            cwd=self.directory,
            env=environment,
            universal_newlines=True)

        self.assertIn("usage: abqhom-install-abaqus", output)


class TestTransactionalInstallation(_InstallerTestCase):
    """Test staged installation, rollback, and managed replacement."""

    def _create_existing_runtime(self):
        """Create a managed runtime containing a recognizable stale file."""
        prefix = os.path.join(self.directory, "runtime")
        os.mkdir(prefix)
        package_directory = self._create_runtime(prefix)
        self._write_valid_marker(prefix, package_directory)
        old_file = os.path.join(prefix, "old-runtime.txt")
        with open(old_file, "w"):
            pass
        return prefix, old_file

    def test_failed_installation_preserves_previous_runtime(self):
        """A failed staged installation leaves the previous runtime intact."""
        prefix, old_file = self._create_existing_runtime()
        original_install = installer._install_packages
        original_requirements = installer._runtime_requirements

        def fail_install(stage, source, requirements):
            """Simulate failure after partially populating staging."""
            with open(os.path.join(stage, "partial.txt"), "w"):
                pass
            raise RuntimeError("simulated installation failure")

        installer._install_packages = fail_install
        installer._runtime_requirements = lambda: TEST_REQUIREMENTS
        try:
            self.assertRaises(RuntimeError, installer.install, prefix)
        finally:
            installer._install_packages = original_install
            installer._runtime_requirements = original_requirements

        self.assertTrue(os.path.isfile(old_file))
        self.assertEqual(installer._read_marker(prefix)["managed_by"],
                         "abqhom")

    def test_failed_validation_preserves_previous_runtime(self):
        """A failed staged validation also preserves the previous runtime."""
        prefix, old_file = self._create_existing_runtime()
        original_install = installer._install_packages
        original_requirements = installer._runtime_requirements

        def incomplete_install(stage, source, requirements):
            """Create a staged runtime that fails content validation."""
            self._create_runtime(stage, include_library=False)

        installer._install_packages = incomplete_install
        installer._runtime_requirements = lambda: TEST_REQUIREMENTS
        try:
            self.assertRaises(RuntimeError, installer.install, prefix)
        finally:
            installer._install_packages = original_install
            installer._runtime_requirements = original_requirements

        self.assertTrue(os.path.isfile(old_file))
        self.assertEqual(installer._read_marker(prefix)["managed_by"],
                         "abqhom")

    def test_successful_update_replaces_previous_runtime(self):
        """A staged update removes stale files from the managed runtime."""
        prefix, old_file = self._create_existing_runtime()
        original_install = installer._install_packages
        original_requirements = installer._runtime_requirements

        def fake_install(stage, source, requirements):
            """Create a complete staged runtime without invoking pip."""
            self._create_runtime(stage)

        installer._install_packages = fake_install
        installer._runtime_requirements = lambda: TEST_REQUIREMENTS
        try:
            package_directory, scripts = installer.install(prefix)
        finally:
            installer._install_packages = original_install
            installer._runtime_requirements = original_requirements

        self.assertFalse(os.path.exists(old_file))
        self.assertTrue(os.path.isdir(package_directory))
        self.assertTrue(all(os.path.isfile(path) for path in scripts))
        self.assertEqual(installer._read_marker(prefix)["site_packages"],
                         package_directory)


if __name__ == "__main__":
    unittest.main()
