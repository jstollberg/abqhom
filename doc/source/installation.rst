Installation
============

Install ``abqhom`` in a regular Python environment and use its installation
helper to create an isolated Abaqus runtime. The Abaqus installation remains
unchanged, which makes the workflow suitable for workstations and compute
clusters without administrative access. Abaqus 2018 or newer must be installed
separately, and the regular installation environment requires Python 3.9 or
newer.

Install the package
-------------------

Until the first PyPI release, install from a source checkout:

.. code-block:: console

   git clone https://github.com/jstollberg/abqhom.git
   cd abqhom
   python -m pip install .
   abqhom-install-abaqus --abqhom-source .

The helper installs only ``abqhom``, Gmsh, Paraqus, and
``backports.weakref`` into a versioned prefix below
``~/.abqhom/abaqus-runtime``. Dependencies are deliberately disabled during
this step. Abaqus therefore continues to use its bundled NumPy installation.
Use ``--prefix PATH`` to select another runtime location. Running
``abqhom-install-abaqus --help`` lists the available options. Re-running the
helper safely replaces an existing runtime managed at the same location.

Activate the runtime
--------------------

The command reports the installed ``site-packages`` path and creates an
activation script in the runtime prefix. Run that script in the shell from
which Abaqus will be launched.

On Windows PowerShell:

.. code-block:: powershell

   . "C:/Users/name/.abqhom/abaqus-runtime/.../activate_abqhom.ps1"

On Linux:

.. code-block:: bash

   . ~/.abqhom/abaqus-runtime/.../activate_abqhom.sh

Successful activation prints ``abqhom activated:`` followed by the managed
runtime package directory.

Only this dedicated runtime directory is added to ``PYTHONPATH``. Other
packages installed in the regular Python environment are not exposed to
Abaqus. For batch systems, the reported path can equivalently be added by the
job scheduler environment.

Verify the setup
----------------

The geometry layer can be checked in a regular Python interpreter:

.. code-block:: console

   python -c "from abqhom.gmshmodel import GmshRVE"

Then run one of the bundled examples from its directory to verify the Abaqus
side:

.. code-block:: console

   cd examples/beam_lattice_simple_cubic_3d
   abaqus cae noGUI=main.py

The activated runtime supplies abqhom, Gmsh, Paraqus, and the legacy
weak-reference backport to the examples.

Linux Gmsh compatibility mode
-----------------------------

Some Linux Abaqus installations terminate at ``gmsh.initialize()`` with a
native-library error such as ``free(): invalid pointer``. If the same managed
runtime works with regular Python but fails only inside Abaqus, reinstall it
with the optional compatibility mode:

.. code-block:: console

   abqhom-install-abaqus --abqhom-source . --gmsh-deepbind

This option changes only how the managed Gmsh shared library is loaded by
Abaqus on Linux. It uses ``RTLD_DEEPBIND`` to isolate Gmsh from conflicting
symbols already loaded by the Abaqus executable. Regular Python, other
platforms, and the Gmsh NumPy interface retain their standard behavior. Use
the option only when the native-library conflict occurs.

The selected mode is recorded as ``gmsh_deepbind`` in the runtime marker
``.abqhom-runtime.json``. Reinstalling the same managed prefix without
``--gmsh-deepbind`` restores the unmodified Gmsh wrapper.

.. tip::

   Do not copy packages into the Abaqus ``site-packages`` directory. The
   isolated runtime and its activation script keep the solver installation
   untouched and make the same analysis scripts portable between machines.
