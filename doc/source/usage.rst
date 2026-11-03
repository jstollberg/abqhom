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

Automatic boundary detection is sufficient for the usual axis-aligned cells.
For imported or unusually partitioned geometry, :doc:`boundary` explains how
to provide the boundary entities explicitly.

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
nonlinearity is defined through the usual Abaqus material and section APIs.
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
