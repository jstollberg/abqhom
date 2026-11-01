API reference
=============

Only the public user interface is listed here. The geometry classes can be
used in a regular Python interpreter with Gmsh. The analysis classes require
the Abaqus CAE Python kernel.

Geometry and mesh
-----------------

.. currentmodule:: abqhom.gmshmodel

.. autoclass:: GmshRVE
   :members:

.. autoclass:: Boundary
   :members:

Abaqus model and analysis
-------------------------

.. currentmodule:: abqhom.abqmodel

.. autoclass:: AbaqusRVE
   :members:

.. autoclass:: LoadHistory
   :members:
   :special-members: __call__, __getitem__, __iter__

Utilities
---------

.. automodule:: abqhom.utils
   :members:
