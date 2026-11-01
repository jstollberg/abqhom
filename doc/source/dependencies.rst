Dependencies
============

Required software
-----------------

**Abaqus**
   Model construction, finite element analysis, and result recovery use the
   Abaqus scripting interface. Scripts that instantiate
   :class:`~abqhom.abqmodel.AbaqusRVE` must therefore run in an Abaqus CAE
   Python kernel, for example with ``abaqus cae noGUI=main.py``. The source is
   compatible with both the Python 2.7 interface of legacy Abaqus releases and
   the Python 3 interface of current releases.

**Gmsh**
   Geometry handling and meshing use the Gmsh Python API. Use the Python
   module supplied with the official `Gmsh SDK <https://gmsh.info/>`_ and
   keep its accompanying shared libraries in place.

**NumPy**
   Array and tensor operations use NumPy. A compatible NumPy distribution is
   included with current Abaqus installations.

Optional software
-----------------

**Paraqus**
   :meth:`~abqhom.abqmodel.AbaqusRVE.export_to_vtk` uses
   `Paraqus <https://github.com/InstituteOfMechanics/Paraqus>`_ to write VTK
   files. All meshing, analysis, and homogenization functionality remains
   available without it.

.. note::

   Abaqus uses its own Python interpreter. Installing a package in a regular
   Python environment does not automatically make it available to Abaqus.
   :doc:`installation` therefore configures source directories through
   ``PYTHONPATH`` rather than modifying the Abaqus installation.
