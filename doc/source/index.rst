abqhom
======

``abqhom`` is an open-source Python framework for first-order computational
homogenization of two- and three-dimensional unit cells and representative
volume elements with `Gmsh <https://gmsh.info/>`_ and
`Abaqus <https://www.3ds.com/products/simulia/abaqus>`_. It supports continuum,
truss, beam, and shell constituents, including mixed-dimensional cells,
finite kinematics, and path-dependent material behavior available in Abaqus.

The package connects geometry and meshing in Gmsh with model
construction, analysis, and effective-response recovery in Abaqus. It also
removes redundant copies of structural elements on periodic cell boundaries
before the boundary-value problem is solved.

Start with :doc:`installation`, then run the complete example in
:doc:`usage`. The :doc:`examples` page presents the tutorial and
manuscript-reproduction cases, while :doc:`api` documents the public Python
interface. Questions and bug reports can be submitted through the
`GitHub issue tracker <https://github.com/jstollberg/abqhom/issues>`_.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   introduction
   dependencies
   installation
   usage

.. toctree::
   :maxdepth: 2
   :caption: User guide

   examples
   boundary
   limitations

.. toctree::
   :maxdepth: 2
   :caption: Reference

   citation
   api
   tests
   development

Indices
-------

* :ref:`genindex`
* :ref:`search`
