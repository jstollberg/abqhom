Installation
============

``abqhom`` is currently used directly from a source checkout. Clone the
repository and add its ``src`` directory, the Gmsh SDK, and optionally
Paraqus to ``PYTHONPATH``. This approach works without write access to the
Abaqus installation and is suitable for local workstations and compute
clusters.

Obtain the source
-----------------

.. code-block:: console

   git clone https://github.com/jstollberg/abqhom.git

Configure the Python path
-------------------------

The Gmsh path below is the directory containing ``gmsh.py`` and the Gmsh
shared library. The Paraqus entry is optional.

On Windows PowerShell:

.. code-block:: powershell

   $env:PYTHONPATH = (
       "C:/path/to/abqhom/src;" +
       "C:/path/to/gmsh/lib;" +
       "C:/path/to/Paraqus/src;" +
       $env:PYTHONPATH
   )

On Linux or macOS:

.. code-block:: bash

   export PYTHONPATH="/path/to/abqhom/src:/path/to/gmsh/lib:/path/to/Paraqus/src:${PYTHONPATH}"

Set the variable in the shell from which Abaqus is launched. For a persistent
installation, define the same entries in the environment setup used by your
workstation or job scheduler.

Verify the setup
----------------

The geometry layer can be checked in a regular Python interpreter:

.. code-block:: console

   python -c "import gmsh; from abqhom.gmshmodel import GmshRVE; print('OK')"

Then run one of the bundled examples from its directory to verify the Abaqus
side:

.. code-block:: console

   cd examples/beam_lattice_simple_cubic_3d
   abaqus cae noGUI=main.py

The bundled examples locate the repository's ``src`` directory relative to
their own files, so only Gmsh and optional Paraqus must be exposed to them
through ``PYTHONPATH``.

.. tip::

   Avoid adding imports to every analysis script or copying packages into the
   Abaqus ``site-packages`` directory. A configured ``PYTHONPATH`` keeps the
   solver installation untouched and makes the same script portable between
   machines.
