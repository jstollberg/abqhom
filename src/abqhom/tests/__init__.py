"""
Unit tests for abqhom.

The tests are split into two groups:

- ``abqhom.tests.tests_common`` contains tests that only need Gmsh (the
  boundary detection, the utility functions, ...). They run in any Python
  interpreter with Gmsh available, including the one shipped with Abaqus.
- ``abqhom.tests.tests_abaqus`` contains tests for AbaqusRVE, i.e. the
  full homogenization workflow. They must be run in the Abaqus Python
  interpreter (``abaqus cae noGUI=...``).

Use ``run_abaqus_tests()`` (in Abaqus Python) to run everything, or
``run_python_tests()`` (in any Python with Gmsh) to run the Gmsh-only
subset. Both write a report file to the current working directory.

If Gmsh is not importable, set the environment variable ``GMSH_LIB`` to
the directory containing ``gmsh.py`` before running the tests.

"""
import os
import sys
import unittest


def _ensure_gmsh():
    """Make sure Gmsh is importable, falling back to $GMSH_LIB."""
    try:
        import gmsh  # noqa: F401
    except ImportError:
        gmsh_lib = os.environ.get("GMSH_LIB")
        if gmsh_lib and gmsh_lib not in sys.path:
            sys.path.insert(0, gmsh_lib)
        import gmsh  # noqa: F401


def run_abaqus_tests():
    """
    Run all unit tests (Gmsh-only and Abaqus) using Abaqus Python.

    A report is written to ``abqhom_test_report_abaqus.txt`` in the
    current working directory.

    Raises
    ------
    RuntimeError
        If Abaqus cannot be imported (i.e. not running in Abaqus Python).

    """
    try:
        import abaqus  # noqa: F401
    except ImportError:
        raise RuntimeError("run_abaqus_tests() must be run in Abaqus Python.")

    _ensure_gmsh()

    from abqhom.tests import tests_common, tests_abaqus
    suites = [_find_tests(m) for m in (tests_common, tests_abaqus)]
    _run_tests(unittest.TestSuite(suites), "abqhom_test_report_abaqus.txt")


def run_python_tests():
    """
    Run the Gmsh-only subset of the unit tests.

    A report is written to ``abqhom_test_report_python.txt`` in the
    current working directory.

    """
    _ensure_gmsh()

    from abqhom.tests import tests_common
    _run_tests(_find_tests(tests_common), "abqhom_test_report_python.txt")


def _find_tests(module):
    """Build a test suite of all tests in a package."""
    loader = unittest.TestLoader()
    return loader.discover(module.__path__[0])


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
