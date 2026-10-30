"""
Unit tests for abqhom.

The tests are split into two groups:

- ``abqhom.tests.tests_common`` contains tests that only need Gmsh (the
  boundary detection, the utility functions, ...). They run in any Python
  interpreter with Gmsh available, including the one shipped with Abaqus.
- ``abqhom.tests.tests_abaqus`` contains tests for AbaqusRVE, i.e. the
  full homogenization workflow. They must be run in the Abaqus Python
  interpreter (``abaqus cae noGUI=...``).

Use ``run_abaqus_tests(gmsh_path)`` (in Abaqus Python) to run everything,
or ``run_python_tests(gmsh_path)`` (in any compatible Python interpreter)
to run the Gmsh-only subset. ``gmsh_path`` must be the directory containing
the prescribed ``gmsh.py`` installation. Both runners isolate generated
files and write a report to the current working directory.

"""
import os
import shutil
import sys
import tempfile
import time
import unittest


def _canonical_path(path):
    """Return a case-normalized physical path for reliable comparison."""
    return os.path.normcase(os.path.realpath(os.path.abspath(path)))


def _load_gmsh(gmsh_path):
    """Load Gmsh exclusively from the prescribed library directory."""

    if gmsh_path is None:
        raise ValueError(
            "A path to the prescribed Gmsh library must be provided.")
    try:
        gmsh_path = os.fspath(gmsh_path)
    except AttributeError:  # Legacy Abaqus releases use Python 2.7
        pass
    gmsh_path = _canonical_path(gmsh_path)
    gmsh_file = os.path.join(gmsh_path, "gmsh.py")
    if not os.path.isfile(gmsh_file):
        raise IOError("The prescribed Gmsh module does not exist: {}"
                      .format(gmsh_file))

    loaded = sys.modules.get("gmsh")
    if loaded is not None:
        loaded_directory = _canonical_path(os.path.dirname(loaded.__file__))
        if loaded_directory != gmsh_path:
            raise RuntimeError(
                "Gmsh is already loaded from '{}', not from the prescribed "
                "directory '{}'. Start a fresh interpreter and provide the "
                "Gmsh path before importing gmsh."
                .format(loaded_directory, gmsh_path))
        return loaded

    sys.path.insert(0, gmsh_path)
    try:
        import gmsh
    finally:
        if sys.path[0] == gmsh_path:
            del sys.path[0]

    loaded_directory = _canonical_path(os.path.dirname(gmsh.__file__))
    if loaded_directory != gmsh_path:
        raise RuntimeError(
            "Gmsh was loaded from '{}', not from the prescribed directory "
            "'{}'.".format(loaded_directory, gmsh_path))
    return gmsh


def run_abaqus_tests(gmsh_path):
    """
    Run all unit tests (Gmsh-only and Abaqus) using Abaqus Python.

    Generated Abaqus files are isolated in a temporary subdirectory. A report
    is written to ``abqhom_test_report_abaqus.txt`` in the current working
    directory. The temporary directory is retained when a test fails.

    Parameters
    ----------
    gmsh_path : str or os.PathLike
        Directory containing the prescribed ``gmsh.py`` installation.

    Raises
    ------
    RuntimeError
        If Abaqus cannot be imported (i.e. not running in Abaqus Python).

    """
    try:
        import abaqus  # noqa: F401
    except ImportError:
        raise RuntimeError("run_abaqus_tests() must be run in Abaqus Python.")

    _load_gmsh(gmsh_path)

    from abqhom.tests import tests_common, tests_abaqus
    suites = [_find_tests(m) for m in (tests_common, tests_abaqus)]
    _run_isolated(
        unittest.TestSuite(suites), "abqhom_test_report_abaqus.txt")


def run_python_tests(gmsh_path):
    """
    Run the Gmsh-only subset of the unit tests.

    Generated files are isolated in a temporary subdirectory. A report is
    written to ``abqhom_test_report_python.txt`` in the current working
    directory. The temporary directory is retained when a test fails.

    Parameters
    ----------
    gmsh_path : str or os.PathLike
        Directory containing the prescribed ``gmsh.py`` installation.

    """
    _load_gmsh(gmsh_path)

    from abqhom.tests import tests_common
    _run_isolated(
        _find_tests(tests_common), "abqhom_test_report_python.txt")


def _find_tests(module):
    """Build a test suite of all tests in a package."""
    loader = unittest.TestLoader()
    return loader.discover(module.__path__[0])


def _run_isolated(suite, report_name):
    """Run a suite in an isolated directory and retain failures."""

    original_directory = os.getcwd()
    report_path = os.path.join(original_directory, report_name)
    work_directory = tempfile.mkdtemp(
        prefix=".abqhom_tests_", dir=original_directory)
    temporary_directory = os.path.join(work_directory, "temp")
    os.mkdir(temporary_directory)
    original_tempdir = tempfile.tempdir

    try:
        tempfile.tempdir = temporary_directory
        os.chdir(work_directory)
        _run_tests(suite, report_path)
    except Exception:
        print("Test artifacts retained in {}".format(work_directory))
        raise
    else:
        os.chdir(original_directory)
        tempfile.tempdir = original_tempdir
        _remove_test_directory(work_directory)
    finally:
        os.chdir(original_directory)
        tempfile.tempdir = original_tempdir


def _remove_test_directory(path, attempts=8, delay=0.25):
    """Remove test artifacts, allowing Abaqus workers time to release them."""
    for attempt in range(attempts):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if attempt + 1 == attempts:
                print("Test artifacts could not be removed from {}"
                      .format(path))
                return
            time.sleep(delay)


def _run_tests(suite, file_name):
    """Run a test suite, write a report and print a short summary."""
    print("Starting test run...")
    logfile = open(file_name, "w")
    try:
        result = unittest.TextTestRunner(logfile, verbosity=3).run(suite)
    finally:
        logfile.close()

    print("Test run finished.\n   tests run: %d\n   failures: %d\n"
          "   errors: %d" % (result.testsRun, len(result.failures),
                             len(result.errors)))

    if not result.wasSuccessful():
        report = open(file_name, "r")
        try:
            for line in report.readlines():
                print(line.rstrip())
        finally:
            report.close()
        raise RuntimeError("Not all tests passed.")
