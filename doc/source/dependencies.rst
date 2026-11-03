Dependencies
============

Required software
-----------------

**Abaqus 2018 or newer**
   Model construction, finite element analysis, and result recovery use the
   Abaqus scripting interface. Scripts that instantiate
   :class:`~abqhom.abqmodel.AbaqusRVE` must therefore run in an Abaqus CAE
   Python kernel, for example with ``abaqus cae noGUI=main.py``. Both legacy
   Python 2.7-based and current Python 3-based Abaqus releases are supported.

**Python 3.9 or newer**
   A regular Python environment is required to install ``abqhom`` and create
   the isolated Abaqus runtime. The Abaqus interpreter itself may still use
   Python 2.7.

**Gmsh**
   Geometry handling and meshing use Gmsh 4.9.4 or newer. The installation
   helper obtains its Python API and shared library from the published wheel.

**NumPy**
   Array and tensor operations use NumPy. A compatible NumPy distribution is
   included with current Abaqus installations.

**Paraqus**
   :meth:`~abqhom.abqmodel.AbaqusRVE.export_to_vtk` uses
   `Paraqus <https://github.com/InstituteOfMechanics/Paraqus>`_ to write VTK
   files. Version 1.0.2 or newer provides the required API.

**backports.weakref**
   The runtime includes this pure-Python backport so current Gmsh releases can
   also be imported by legacy Abaqus versions using Python 2.7. It is harmless
   when the standard-library implementation is available.

.. note::

   Abaqus uses its own Python interpreter. Installing a package in a regular
   Python environment does not automatically make it available to Abaqus.
   :doc:`installation` therefore creates a dedicated runtime prefix and adds
   only that prefix to ``PYTHONPATH``.
