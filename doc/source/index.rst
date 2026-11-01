abqhom
======

``abqhom`` is a Python toolbox for the computational homogenization of
periodic two- and three-dimensional unit cells with
`Gmsh <https://gmsh.info/>`_ and
`Abaqus <https://www.3ds.com/products/simulia/abaqus>`_. It supports
continuum, truss, beam, and shell constituents, including mixed-dimensional
cells, finite kinematics, and path-dependent material behavior available in
Abaqus.

The package connects geometry and periodic meshing in Gmsh with model
construction, analysis, and effective-response recovery in Abaqus. It also
removes redundant copies of structural elements on periodic cell boundaries
before the boundary-value problem is solved.

Start with :doc:`installation`, then run the complete example in
:doc:`usage`. The :doc:`examples` page provides larger tutorial problems, and
:doc:`api` documents the public Python interface.

Getting started
---------------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   introduction
   dependencies
   installation
   usage

User guide
----------

.. toctree::
   :maxdepth: 2
   :caption: User guide

   examples
   limitations

Reference
---------

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   tests
   development

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
