Overview
========

``abqhom`` automates periodic-cell analyses without prescribing a particular
constituent material model. The finite element response is computed by
Abaqus, while the package supplies the scale-transition machinery around that
analysis.

Core capabilities
-----------------

* **Mixed-dimensional unit cells.** Continuum, truss, beam, and shell
  constituents can occur separately or in the same model. Embedded-region
  constraints can couple nonconforming structural and continuum meshes.
* **Linear and nonlinear analyses.** The loading history is expressed through
  the macroscopic displacement gradient. Finite kinematics and nonlinear or
  path-dependent material models are handled by the corresponding Abaqus
  analysis step and constitutive definition.
* **Periodic and affine boundary conditions.** Matching periodic meshes are
  generated in Gmsh or verified for imported meshes. Periodic translations
  and, where present, nodal rotations are constrained automatically. Affine
  displacement conditions are available for continuum cells.
* **Consistent periodic-cell representation.** Redundant structural elements
  on secondary boundary entities are removed before analysis while their
  nodes are retained for the periodic constraints.
* **Effective-response recovery.** The effective first Piola--Kirchhoff stress
  is evaluated from element nodal forces. Incremental microscopic work can be
  recovered independently for a discrete Hill--Mandel check.
* **Visualization.** Abaqus results can be exported to VTK through
  `Paraqus <https://github.com/InstituteOfMechanics/Paraqus>`_. Stored beam
  and truss radii can be included for visualization.

Execution model
---------------

The workflow has two layers:

#. :class:`~abqhom.gmshmodel.GmshRVE` operates on a Gmsh model. It identifies
   the cell boundary, imposes periodic mesh relations, generates or reads the
   mesh, and creates the node and element sets needed by Abaqus.
#. :class:`~abqhom.abqmodel.AbaqusRVE` imports that representation into
   Abaqus and provides methods for section assignment, macroscopic boundary
   conditions, analysis execution, and effective-response recovery from the
   output database.

The geometry layer needs only Gmsh and NumPy. Importing or using the Abaqus
layer requires the Abaqus Python kernel. This separation is useful when a
mesh is prepared outside the solver, although the bundled examples execute
the complete workflow in one Abaqus CAE session.

Cell convention
---------------

The reference cell must be rectangular in the global :math:`X_1`--:math:`X_2`
plane with :math:`X_3=0` in two dimensions or cuboidal and aligned with the
global coordinate axes in three dimensions. This restriction also applies
when the boundary is supplied manually, in which case its geometric
assignment is the user's responsibility. Physical groups in the Gmsh model
identify the constituents to which Abaqus sections are later assigned. See
:doc:`limitations` for the precise modeling restrictions.
