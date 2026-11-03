Quickstart
==========

This example creates a square, homogeneous plane-strain cell, applies a
periodic uniaxial strain, and evaluates its effective stress. Save it as
``quickstart.py`` and run it from a writable working directory with
``abaqus cae noGUI=quickstart.py``.

.. code-block:: python

   import numpy as np
   import gmsh

   from abqhom.abqmodel import AbaqusRVE, LoadHistory
   from abqhom.gmshmodel import GmshRVE
   from abqhom.utils import tensor_to_voigt, voigt_to_tensor


   gmsh.initialize()
   try:
       # Create a 2D cell and identify its material region.
       gmsh.model.add("CELL")
       surface = gmsh.model.occ.add_rectangle(-0.5, -0.5, 0.0, 1.0, 1.0)
       gmsh.model.occ.synchronize()
       group = gmsh.model.add_physical_group(2, [surface])
       gmsh.model.set_physical_name(2, group, "MATERIAL")

       # Generate a matching periodic mesh.
       gmsh_cell = GmshRVE("CELL")
       gmsh_cell.add_periodic_mesh_constraints()
       gmsh_cell.generate_mesh(element_size=0.1)

       # Import the mesh into Abaqus and define the constituent material.
       cell = AbaqusRVE(gmsh_cell)
       material = cell.model.Material(name="MATERIAL")
       material.Elastic(table=[[100.0, 0.3]])
       cell.assign_solid_section("MATERIAL", "MATERIAL")
       cell.model.StaticStep(name="LOAD", previous="Initial")

       # Ramp the macroscopic strain eps_11 to 3 percent.
       strain_magnitude = 0.03

       def displacement_gradient(time):
           strain = np.array([strain_magnitude, 0.0, 0.0])
           return time*voigt_to_tensor(strain, "STRAIN")

       history = LoadHistory([0.0, 1.0], displacement_gradient)
       cell.apply_periodic_bc("LOAD", history)
       cell.run_analysis(field_outputs=["S"])

       # Recover the effective first Piola--Kirchhoff stress.
       stress = cell.average_stress("LOAD")
       print(tensor_to_voigt(stress, "STRESS"))
   finally:
       if gmsh.is_initialized():
           gmsh.finalize()

What happens in the quickstart
------------------------------

#. A Gmsh physical group named ``MATERIAL`` identifies the elements that
   receive the Abaqus section with the same name.
#. :class:`~abqhom.gmshmodel.GmshRVE` detects the axis-aligned cell boundary.
   :meth:`~abqhom.gmshmodel.GmshRVE.add_periodic_mesh_constraints` makes the
   meshes on opposite boundaries match before mesh generation.
#. :class:`~abqhom.abqmodel.AbaqusRVE` writes and imports the mesh, transfers
   all generated sets, and creates the Abaqus part and assembly instance.
#. :class:`~abqhom.abqmodel.LoadHistory` stores the total-time history of the
   macroscopic displacement gradient
   :math:`\bar{\boldsymbol{H}}=\bar{\boldsymbol{F}}-\boldsymbol{I}`. In this
   geometrically linear example it equals the prescribed infinitesimal strain.
#. :meth:`~abqhom.abqmodel.AbaqusRVE.apply_periodic_bc` creates the periodic
   equations and the macroscopic loading reference points. A rigid-body anchor
   is selected automatically.
#. :meth:`~abqhom.abqmodel.AbaqusRVE.average_stress` evaluates the effective
   first Piola--Kirchhoff stress from the element nodal-force virial.

.. _manual-boundary:

Supplying a boundary manually
-----------------------------

Normally, :class:`~abqhom.gmshmodel.GmshRVE` constructs its
:class:`~abqhom.gmshmodel.Boundary` from the model geometry. If unusual or
trimmed CAD entities prevent automatic classification, supply a boundary map
explicitly. Its keys are the boundary-group names and its values are lists of
Gmsh dimension tags ``(dimension, tag)``, where dimensions 0, 1, and 2 denote
points, curves, and surfaces. A group may contain several entities. It may
also contain lower-dimensional constituent entities whose complete support
lies on that group, such as a beam curve on a cell face.

A manual map bypasses the automatic geometric classification but does not
extend the supported reference-cell geometries. The cell must remain
rectangular in the global :math:`X_1`--:math:`X_2` plane with :math:`X_3=0`
in two dimensions or axis-aligned and cuboidal in three dimensions. The
mapping's geometric assignment is the user's responsibility. The constituent
topology inside and on that boundary may have arbitrary orientations.

The numbering follows the axis-aligned cell coordinates. For bounds
:math:`(x_{\min},y_{\min},z_{\min})` and
:math:`(x_{\max},y_{\max},z_{\max})`, the groups are:

* ``POINT_1``--``POINT_8`` at, respectively,
  ``(xmin,ymin,zmin)``, ``(xmax,ymin,zmin)``, ``(xmax,ymax,zmin)``,
  ``(xmin,ymax,zmin)``, ``(xmin,ymin,zmax)``, ``(xmax,ymin,zmax)``,
  ``(xmax,ymax,zmax)``, and ``(xmin,ymax,zmax)``;
* ``EDGE_1``--``EDGE_12`` connecting, respectively, points ``1-2``, ``2-3``,
  ``3-4``, ``4-1``, ``5-6``, ``6-7``, ``7-8``, ``8-5``, ``1-5``, ``2-6``,
  ``3-7``, and ``4-8``;
* ``FACE_1`` and ``FACE_2`` at ``zmin`` and ``zmax``, ``FACE_3`` and
  ``FACE_4`` at ``ymin`` and ``ymax``, and ``FACE_5`` and ``FACE_6`` at
  ``xmin`` and ``xmax``.

A two-dimensional cell in the global x-y plane uses ``POINT_1``--``POINT_4``
and ``EDGE_1``--``EDGE_4`` from the same convention. The static method
:meth:`~abqhom.gmshmodel.Boundary.group_coordinate_identifiers` is the
authoritative way to obtain these coordinate signatures:

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
``(nan, ymin, zmin)`` and ``FACE_1`` has ``(nan, nan, zmin)``. Use these
signatures, the Gmsh GUI, or existing physical groups to associate the model's
entity tags with the boundary groups. ``gmsh.model.get_entities(dimension)``
returns the candidate tags of a given dimension. Include every group key in
the mapping, using an empty list when a group has no associated entity:

.. code-block:: python

   identifiers = Boundary.group_coordinate_identifiers(
       bounds_min, bounds_max,
       tolerance=gmsh.option.get_number("Geometry.Tolerance"))
   boundary_data = {name: [] for name in identifiers}

   boundary_data["POINT_1"] = [(0, point_1_tag)]
   boundary_data["EDGE_1"] = [(1, edge_1_tag)]
   boundary_data["FACE_1"] = [(2, face_1_tag)]
   # Assign the remaining boundary entities in the same way.

   boundary = Boundary("CELL", initial_data=boundary_data)
   gmsh_cell = GmshRVE("CELL", boundary=boundary)

With the default ``check_periodic=True``, :class:`~abqhom.gmshmodel.Boundary`
checks the supplied geometry and orders corresponding primary and secondary
entities automatically. If that check also fails, set ``check_periodic=False``
and ``periodic_default=True`` only after arranging corresponding entities at
the same list indices in every periodic primary--secondary group pair.

Effective stiffness
-------------------

For a linear three-dimensional problem, repeat the analysis for the six
independent strain vectors in the Voigt ordering
:math:`(11,22,33,23,13,12)`. A two-dimensional problem requires the three
components :math:`(11,22,12)`. Each effective stress vector divided by the
applied strain magnitude forms one column of the effective stiffness matrix.
The tutorial examples implement these loops and include numerical or
analytical checks where available.

Nonlinear analyses
------------------

For finite kinematics, let the load history return
:math:`\bar{\boldsymbol{H}}(\tau)=\bar{\boldsymbol{F}}(\tau)-\boldsymbol{I}`
and create an Abaqus step with geometric nonlinearity enabled. Constitutive
nonlinearity is defined through the usual Abaqus material and section APIs;
``abqhom`` does not replace or restrict that constitutive model.

The following methods support nonlinear post-processing:

* :meth:`~abqhom.abqmodel.AbaqusRVE.average_stress` returns the effective
  first Piola--Kirchhoff stress at a requested output frame.
* :meth:`~abqhom.abqmodel.AbaqusRVE.average_power_increment` evaluates the
  microscopic incremental-work density between consecutive frames from nodal
  translations, rotations, forces, and moments.
* :func:`~abqhom.utils.hill_mandel_residual` compares independently evaluated
  macroscopic and microscopic incremental-work densities.

See :doc:`examples` for complete finite-strain and path-dependent workflows.
