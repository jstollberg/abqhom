Quickstart
==========

This example creates a square, homogeneous plane-strain cell, applies a
periodic uniaxial strain, and evaluates its effective stress. Save it as
``quickstart.py`` and run it from a writable working directory with
``abaqus cae noGUI=quickstart.py``.

.. code-block:: python

   import sys

   import numpy as np
   import gmsh

   from abqhom.abqmodel import AbaqusRVE, LoadHistory
   from abqhom.gmshmodel import GmshRVE
   from abqhom.utils import tensor_to_voigt, voigt_to_tensor


   # Abaqus replaces stdout. Restore output to the launching shell.
   if sys.__stdout__ is not None:
       sys.stdout = sys.__stdout__

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
       cell.apply_periodic_displacement_bc("LOAD", history)
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
#. :meth:`~abqhom.abqmodel.AbaqusRVE.apply_periodic_displacement_bc` creates
   the periodic equations and the macroscopic loading reference points. A
   rigid-body anchor is selected automatically.
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

Stress control
--------------

Periodic cells can alternatively be loaded by a prescribed effective stress
through periodic traction boundary conditions. The method introduces shared
macroscopic kinematic controls and applies their work-conjugate generalized
forces. For a control :math:`q_\alpha` with tensor basis
:math:`\boldsymbol{B}_\alpha`, this force is

.. math::

   Q_\alpha(\tau)
   =V\,\bar{\boldsymbol{T}}(\tau)
   :\boldsymbol{B}_\alpha.

In a geometrically linear step, :math:`\bar{\boldsymbol{T}}` is the Cauchy
stress and the controls are engineering strains. Their shear loads are
therefore :math:`V\bar{\sigma}_{ij}` without an additional factor of two. In
a geometrically nonlinear step, :math:`\bar{\boldsymbol{T}}` is the first
Piola--Kirchhoff stress. An upper-triangular deformation gradient removes the
arbitrary macroscopic rigid rotation while retaining six independent
deformation components in three dimensions and three in two dimensions.
Components with identically zero generalized loads remain free and enforce
zero effective stress in their respective directions.

.. code-block:: python

   def stress_history(time):
       stress = np.zeros((3, 3))
       stress[0, 0] = time*1.0
       return stress

   history = LoadHistory([0.0, 1.0], stress_history)
   cell.apply_periodic_traction_bc("LOAD", history)
   cell.run_analysis(field_outputs=["U"])

   stress = cell.average_stress("LOAD")
   deformation_gradient = cell.average_deformation_gradient("LOAD")

The deformation gradient is assembled from the solved macroscopic control
DOFs. In the geometrically linear example above, the infinitesimal strain
follows from the symmetric part of its displacement gradient. For a
geometrically nonlinear step, the history function instead returns the
effective first Piola--Kirchhoff stress. Its upper-triangular components are
work-conjugate to the independent deformation-gradient controls. The
remaining stress components arise as reactions to the rotation constraint
and are recovered by :meth:`~abqhom.abqmodel.AbaqusRVE.average_stress`.

The effective stress obtained independently from
:meth:`~abqhom.abqmodel.AbaqusRVE.average_stress` can be compared with the
prescribed work-conjugate components as an equilibrium check. Mixed
stress--deformation control is not currently implemented.

The simple-cubic beam-lattice tutorial demonstrates linear displacement and
stress control. The spherical-inclusion tutorial recovers a
displacement-controlled finite-strain state under the corresponding stress.

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

Efficient ODB post-processing
-----------------------------

Post-processing methods open and close the output database automatically for
individual calls. For repeated operations, open it once around the complete
loop. Methods then reuse the open database and leave it open until the caller
closes it:

.. code-block:: python

   cell.open_odb_file()
   try:
       stresses = [
           cell.average_stress("LOAD", frame=index)
           for index in range(cell.get_number_of_frames("LOAD"))
       ]
   finally:
       cell.close_odb_file()

See :doc:`examples` for complete finite-strain and path-dependent workflows.
