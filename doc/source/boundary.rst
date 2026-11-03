.. _manual-boundary:

Manual boundary definition
==========================

Normally, :class:`~abqhom.gmshmodel.GmshRVE` constructs its
:class:`~abqhom.gmshmodel.Boundary` from the model geometry. If unusual or
trimmed CAD entities prevent automatic classification, a boundary map can be
supplied explicitly. Its keys are boundary-group names and its values are
lists of Gmsh dimension tags ``(dimension, tag)``, where dimensions 0, 1, and
2 denote points, curves, and surfaces. A group may contain several entities,
including lower-dimensional constituent entities whose complete support lies
on that group, such as a beam curve on a cell face.

A manual map bypasses automatic classification but does not extend the
supported cell geometries. A two-dimensional cell must remain rectangular in
the global :math:`X_1`--:math:`X_2` plane with :math:`X_3=0`. A
three-dimensional cell must remain axis-aligned and cuboidal. The geometric
assignment of the supplied entities is the user's responsibility.

Boundary-group convention
-------------------------

For bounds :math:`(x_{\min},y_{\min},z_{\min})` and
:math:`(x_{\max},y_{\max},z_{\max})`, the groups are:

* ``POINT_1``--``POINT_8`` at, respectively,
  ``(xmin,ymin,zmin)``, ``(xmax,ymin,zmin)``, ``(xmax,ymax,zmin)``,
  ``(xmin,ymax,zmin)``, ``(xmin,ymin,zmax)``, ``(xmax,ymin,zmax)``,
  ``(xmax,ymax,zmax)``, and ``(xmin,ymax,zmax)``.
* ``EDGE_1``--``EDGE_12`` connecting, respectively, points ``1-2``, ``2-3``,
  ``3-4``, ``4-1``, ``5-6``, ``6-7``, ``7-8``, ``8-5``, ``1-5``, ``2-6``,
  ``3-7``, and ``4-8``.
* ``FACE_1`` and ``FACE_2`` at ``zmin`` and ``zmax``, ``FACE_3`` and
  ``FACE_4`` at ``ymin`` and ``ymax``, and ``FACE_5`` and ``FACE_6`` at
  ``xmin`` and ``xmax``.

A two-dimensional cell uses ``POINT_1``--``POINT_4`` and
``EDGE_1``--``EDGE_4`` from the same convention. The authoritative coordinate
signatures are returned by
:meth:`~abqhom.gmshmodel.Boundary.group_coordinate_identifiers`:

.. code-block:: python

   import numpy as np
   import gmsh

   from abqhom.gmshmodel import Boundary, GmshRVE


   bounds_min = np.array([xmin, ymin, zmin])
   bounds_max = np.array([xmax, ymax, zmax])
   identifiers = Boundary.group_coordinate_identifiers(
       bounds_min, bounds_max,
       tolerance=gmsh.option.get_number("Geometry.Tolerance"))

   for name, coordinates in identifiers.items():
       print(name, coordinates)

An identifier component equal to ``numpy.nan`` is free along the corresponding
coordinate direction. For example, ``EDGE_1`` has the signature
``(nan, ymin, zmin)`` and ``FACE_1`` has ``(nan, nan, zmin)``.

Constructing the map
--------------------

Use the coordinate signatures, Gmsh GUI, or existing physical groups to
associate entity tags with the boundary groups. Include every group key in the
mapping and use an empty list when no entity belongs to a group:

.. code-block:: python

   boundary_data = {name: [] for name in identifiers}

   boundary_data["POINT_1"] = [(0, point_1_tag)]
   boundary_data["EDGE_1"] = [(1, edge_1_tag)]
   boundary_data["FACE_1"] = [(2, face_1_tag)]
   # Assign the remaining boundary entities in the same way.

   boundary = Boundary("CELL", initial_data=boundary_data)
   gmsh_cell = GmshRVE("CELL", boundary=boundary)

With the default ``check_periodic=True``,
:class:`~abqhom.gmshmodel.Boundary` checks the supplied geometry and orders
corresponding primary and secondary entities automatically. Set
``check_periodic=False`` and ``periodic_default=True`` only when the automatic
check fails and corresponding entities have already been placed at matching
list indices in every periodic primary--secondary group pair.
