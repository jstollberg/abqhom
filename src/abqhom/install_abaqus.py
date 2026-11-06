# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Install the packages required by abqhom into an isolated Abaqus runtime."""

import argparse
import fnmatch
from importlib import metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from abqhom import __version__


RUNTIME_PROJECTS = ("gmsh", "paraqus", "backports.weakref")
MARKER_FILE = ".abqhom-runtime.json"
MARKER_SCHEMA = 1


def _canonical_name(name):
    """Return the normalized project name defined by Python packaging."""

    return re.sub(r"[-_.]+", "-", name).lower()


def _requirement_name(requirement):
    """Extract and normalize a project name from a requirement string."""

    match = re.match(r"\s*([A-Za-z0-9][A-Za-z0-9._-]*)", requirement)
    if match is None:
        raise RuntimeError("Invalid package requirement: {}".format(
            requirement))
    return _canonical_name(match.group(1))


def _runtime_requirements():
    """Read the Abaqus-facing requirements from installed package metadata."""

    try:
        requirements = metadata.requires("abqhom") or []
    except metadata.PackageNotFoundError:
        raise RuntimeError(
            "abqhom package metadata is unavailable. Install abqhom before "
            "running the Abaqus installer.")

    # The regular package metadata remains the single source of version
    # constraints. Select only packages that must be visible inside Abaqus.
    requested = {_canonical_name(name): name for name in RUNTIME_PROJECTS}
    selected = {}
    for requirement in requirements:
        name = _requirement_name(requirement)
        if name in requested:
            selected[name] = requirement

    missing = [name for name in requested if name not in selected]
    if missing:
        raise RuntimeError(
            "Package metadata does not define the required Abaqus runtime "
            "project(s): {}.".format(", ".join(sorted(missing))))
    return tuple(selected[_canonical_name(name)] for name in RUNTIME_PROJECTS)


def _default_prefix():
    """Return the versioned default directory for the Abaqus runtime."""

    # Separate runtimes by abqhom and installer-Python version so an update
    # cannot silently replace an environment created for another release.
    python_tag = "py{}.{}".format(sys.version_info[0], sys.version_info[1])
    return os.path.join(os.path.expanduser("~"), ".abqhom",
                        "abaqus-runtime", __version__ + "-" + python_tag)


def _pip_command(prefix, abqhom_source, requirements=None):
    """Build the pip command for an isolated, dependency-free prefix."""

    source = abqhom_source
    if source is None:
        source = "abqhom=={}".format(__version__)
    if requirements is None:
        requirements = _runtime_requirements()

    # Install the selected runtime projects explicitly. --no-deps prevents
    # pip from adding packages such as NumPy that must come from Abaqus.
    return [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--ignore-installed",
        "--no-deps",
        "--prefix",
        prefix,
        source,
    ] + list(requirements)


def _canonical_path(path):
    """Return an absolute, case-normalized physical path."""

    return os.path.normcase(os.path.realpath(os.path.abspath(path)))


def _is_within(path, directory):
    """Return whether path lies in directory or is the directory itself."""

    path = _canonical_path(path)
    directory = _canonical_path(directory)
    return path == directory or path.startswith(directory + os.sep)


def _read_marker(prefix):
    """Read and validate the ownership marker of a managed prefix."""

    marker = os.path.join(prefix, MARKER_FILE)
    try:
        with open(marker, "r") as source:
            data = json.load(source)
    except (IOError, OSError, ValueError) as error:
        raise RuntimeError(
            "The runtime prefix has an invalid ownership marker: {}"
            .format(error))

    # The marker is both an ownership record and a guard against replacing an
    # unrelated directory. Do not trust paths read from it without checking.
    if not isinstance(data, dict):
        raise RuntimeError("The runtime ownership marker must contain a "
                           "JSON object.")
    if data.get("managed_by") != "abqhom":
        raise RuntimeError("The runtime prefix is not marked as managed by "
                           "abqhom.")
    if data.get("schema_version") != MARKER_SCHEMA:
        raise RuntimeError("The runtime prefix uses an unsupported marker "
                           "schema.")
    marker_prefix = data.get("prefix")
    package_directory = data.get("site_packages")
    if not isinstance(marker_prefix, str):
        raise RuntimeError("The runtime marker contains an invalid prefix.")
    if not isinstance(package_directory, str):
        raise RuntimeError("The runtime marker contains an invalid package "
                           "directory.")
    if _canonical_path(marker_prefix) != _canonical_path(prefix):
        raise RuntimeError("The runtime marker belongs to a different "
                           "installation prefix.")
    if not _is_within(package_directory, prefix):
        raise RuntimeError("The runtime marker contains an invalid package "
                           "directory.")
    return data


def _prepare_prefix(prefix):
    """Validate the destination and create its parent directory."""

    prefix = os.path.abspath(os.path.expanduser(prefix))

    # Refuse broad or ambiguous replacement targets. A nonempty destination
    # is reusable only when a valid abqhom ownership marker is present.
    if prefix == os.path.dirname(prefix):
        raise RuntimeError("A filesystem root cannot be a runtime prefix.")
    if os.path.islink(prefix):
        raise RuntimeError("A symbolic link cannot be a runtime prefix: {}"
                           .format(prefix))
    if os.path.isdir(prefix):
        if os.listdir(prefix):
            _read_marker(prefix)
    elif os.path.exists(prefix):
        raise RuntimeError(
            "The runtime prefix exists but is not a directory: {}"
            .format(prefix))

    parent = os.path.dirname(prefix)
    if not os.path.isdir(parent):
        os.makedirs(parent)
    return prefix


def _find_site_packages(prefix):
    """Locate the installed package directory below the prefix."""

    # A pip --prefix layout varies between platforms and Python versions.
    # Locate it through the installed packages instead of assuming a path.
    matches = []
    for root, directories, files in os.walk(prefix):
        if os.path.basename(root).lower() != "abqhom":
            continue
        if "__init__.py" not in files:
            continue

        package_directory = os.path.abspath(os.path.dirname(root))
        required = (
            os.path.join(package_directory, "gmsh.py"),
            os.path.join(package_directory, "paraqus", "__init__.py"),
        )
        if all(os.path.isfile(path) for path in required):
            matches.append(package_directory)
        directories[:] = []

    matches = sorted(set(matches))
    if len(matches) != 1:
        raise RuntimeError(
            "Expected one complete runtime package directory below {}, "
            "found {}.".format(prefix, len(matches)))
    return matches[0]


def _gmsh_libraries(prefix):
    """Return native Gmsh libraries installed below the runtime prefix."""

    patterns = ("gmsh*.dll", "libgmsh*.so*", "libgmsh*.dylib")
    libraries = []
    for root, _, files in os.walk(prefix):
        for name in files:
            if any(fnmatch.fnmatch(name.lower(), pattern)
                   for pattern in patterns):
                libraries.append(os.path.join(root, name))
    return libraries


def _check_runtime_contents(prefix, package_directory):
    """Validate required packages and exclude an external NumPy."""

    weakref_candidates = (
        os.path.join(package_directory, "backports", "weakref.py"),
        os.path.join(package_directory, "backports", "weakref", "__init__.py"),
    )
    if not any(os.path.isfile(path) for path in weakref_candidates):
        raise RuntimeError("The runtime prefix does not contain "
                           "backports.weakref.")
    if not _gmsh_libraries(prefix):
        raise RuntimeError("The runtime prefix does not contain the native "
                           "Gmsh shared library.")

    # An external NumPy can shadow the version bundled with Abaqus and cause
    # binary incompatibilities when the runtime is added to PYTHONPATH.
    for name in os.listdir(package_directory):
        normalized = name.lower()
        if (normalized in ("numpy", "numpy.py")
                or normalized.startswith(("numpy-", "numpy_", "numpy."))):
            raise RuntimeError(
                "The runtime prefix contains an external NumPy package. "
                "Abaqus must use its bundled NumPy installation.")


def _check_runtime_imports(prefix, package_directory):
    """Import isolated runtime packages and validate the abqhom version."""

    # Run the check in a clean subprocess with only the staged package
    # directory on PYTHONPATH. Each required module must resolve inside it.
    script = """
import os
import sys
import abqhom
import gmsh
import paraqus
import backports.weakref as weakref_backport

root = os.path.normcase(os.path.realpath(os.path.abspath(sys.argv[1])))
for module in (abqhom, gmsh, paraqus, weakref_backport):
    path = os.path.normcase(os.path.realpath(os.path.abspath(module.__file__)))
    if path != root and not path.startswith(root + os.sep):
        raise RuntimeError("Module loaded outside the runtime: " + path)
print("ABQHOM_VERSION=" + abqhom.__version__)
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = package_directory
    output = subprocess.check_output(
        [sys.executable, "-c", script, package_directory],
        cwd=prefix,
        env=environment,
        universal_newlines=True)
    version_lines = [line for line in output.splitlines()
                     if line.startswith("ABQHOM_VERSION=")]
    if version_lines != ["ABQHOM_VERSION=" + __version__]:
        raise RuntimeError(
            "The installed abqhom version does not match the installer "
            "version {}.".format(__version__))


def _powershell_quote(value):
    """Quote a literal path for a PowerShell single-quoted string."""

    return "'{}'".format(value.replace("'", "''"))


def _shell_quote(value):
    """Quote a literal path for a POSIX single-quoted string."""

    return "'{}'".format(value.replace("'", "'\"'\"'"))


def _write_activation_scripts(directory, package_directory):
    """Write platform-specific scripts that prepend the runtime path."""

    # Add only the dedicated runtime package directory. Exposing the complete
    # installer environment could shadow modules bundled with Abaqus.
    paths = []
    if os.name == "nt":
        powershell = os.path.join(directory, "activate_abqhom.ps1")
        with open(powershell, "w") as output:
            output.write(
                "$abqhomSitePackages = {}\n"
                "if ([string]::IsNullOrEmpty($env:PYTHONPATH)) {{\n"
                "    $env:PYTHONPATH = $abqhomSitePackages\n"
                "}} else {{\n"
                "    $env:PYTHONPATH = $abqhomSitePackages + ';' + "
                "$env:PYTHONPATH\n"
                "}}\n".format(_powershell_quote(package_directory)))
        paths.append(powershell)

        batch = os.path.join(directory, "activate_abqhom.bat")
        with open(batch, "w") as output:
            output.write(
                "@echo off\n"
                "if defined PYTHONPATH (\n"
                "    set \"PYTHONPATH={};%PYTHONPATH%\"\n"
                ") else (\n"
                "    set \"PYTHONPATH={}\"\n"
                ")\n".format(package_directory, package_directory))
        paths.append(batch)
    else:
        shell = os.path.join(directory, "activate_abqhom.sh")
        with open(shell, "w") as output:
            output.write(
                "#!/bin/sh\n"
                "abqhom_site_packages={}\n"
                "if [ -n \"${{PYTHONPATH:-}}\" ]; then\n"
                "    export PYTHONPATH=\"$abqhom_site_packages:$PYTHONPATH\"\n"
                "else\n"
                "    export PYTHONPATH=\"$abqhom_site_packages\"\n"
                "fi\n".format(_shell_quote(package_directory)))
        os.chmod(shell, 0o755)
        paths.append(shell)
    return paths


def _write_marker(directory, prefix, package_directory, requirements):
    """Record the managed runtime layout for safe repeated installation."""

    data = {
        "abqhom_version": __version__,
        "managed_by": "abqhom",
        "prefix": prefix,
        "requirements": list(requirements),
        "schema_version": MARKER_SCHEMA,
        "site_packages": package_directory,
    }
    marker = os.path.join(directory, MARKER_FILE)
    with open(marker, "w") as output:
        json.dump(data, output, indent=2, sort_keys=True)
        output.write("\n")


def _install_packages(prefix, abqhom_source, requirements):
    """Install the selected distributions into a staging prefix."""

    subprocess.check_call(
        _pip_command(prefix, abqhom_source, requirements=requirements))


def _promote_stage(stage, prefix):
    """Replace a managed runtime with a validated staging directory."""

    backup = stage + ".previous"
    had_previous = os.path.lexists(prefix)

    # Rename the previous runtime rather than deleting it first. It can then
    # be restored if promoting the staged installation fails.
    if had_previous:
        os.rename(prefix, backup)
    try:
        os.rename(stage, prefix)
    except Exception:
        if had_previous and not os.path.lexists(prefix):
            os.rename(backup, prefix)
        raise
    if had_previous:
        try:
            shutil.rmtree(backup)
        except OSError as error:
            print("Warning: the previous runtime could not be removed from "
                  "{}: {}".format(backup, error), file=sys.stderr)


def install(prefix, abqhom_source=None):
    """
    Install and configure an isolated Python runtime for Abaqus.

    Parameters
    ----------
    prefix : str or os.PathLike
        Destination directory of the managed runtime.
    abqhom_source : str or os.PathLike, optional
        Wheel, source directory, or package requirement used to install
        abqhom. By default, the running abqhom version is installed from the
        configured package index.

    Returns
    -------
    package_directory : str
        Directory to prepend to ``PYTHONPATH`` before launching Abaqus.
    scripts : list of str
        Paths to the generated platform-specific activation scripts.

    Raises
    ------
    OSError
        If the runtime directories or activation files cannot be managed.
    RuntimeError
        If the destination is unmanaged or the staged runtime is invalid.
    subprocess.CalledProcessError
        If pip cannot install the selected packages.

    """

    prefix = _prepare_prefix(prefix)
    requirements = _runtime_requirements()
    parent = os.path.dirname(prefix)

    # Build beside the final prefix so the final directory replacement stays
    # on one filesystem and can use atomic rename operations where supported.
    stage = tempfile.mkdtemp(
        prefix=".{}-staging-".format(os.path.basename(prefix)), dir=parent)

    try:
        # Populate and validate staging before modifying the active runtime.
        _install_packages(stage, abqhom_source, requirements)
        stage_packages = _find_site_packages(stage)
        _check_runtime_contents(stage, stage_packages)
        _check_runtime_imports(stage, stage_packages)

        # Write activation data for the final location, not the staging path.
        relative_packages = os.path.relpath(stage_packages, stage)
        final_packages = os.path.join(prefix, relative_packages)
        stage_scripts = _write_activation_scripts(stage, final_packages)
        final_scripts = [
            os.path.join(prefix, os.path.relpath(path, stage))
            for path in stage_scripts
        ]
        _write_marker(
            stage, prefix, final_packages, requirements=requirements)

        # Replace an existing managed runtime only after validation succeeds.
        _promote_stage(stage, prefix)
        stage = None
    finally:
        # Failed installations never modify the active runtime. Remove only
        # the temporary directory still owned by this invocation.
        if stage is not None and os.path.isdir(stage):
            shutil.rmtree(stage)

    return final_packages, final_scripts


def _parser():
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="abqhom-install-abaqus",
        description=(
            "Install abqhom and its required Abaqus-facing packages into "
            "an isolated prefix without installing an external NumPy."))
    parser.add_argument(
        "--prefix",
        default=_default_prefix(),
        help="runtime installation prefix (default: %(default)s)")
    parser.add_argument(
        "--abqhom-source",
        help=(
            "abqhom wheel, source directory or requirement to install; "
            "the running version is used by default"))
    return parser


def main(argv=None):
    """Install the Abaqus runtime and print its activation instructions."""

    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        package_directory, scripts = install(
            arguments.prefix, arguments.abqhom_source)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "abqhom installation failed: {}\n".format(error))

    print("Abaqus runtime installed successfully.")
    print("PYTHONPATH entry: {}".format(package_directory))
    print("Activation script{}:".format("s" if len(scripts) != 1 else ""))
    for path in scripts:
        print("  {}".format(path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
