Testing
========

The tests under ``src/abqhom/tests`` are split into two groups:

``tests_common``
   Geometry, meshing, boundary classification, utilities, and other tests
   that do not require the Abaqus kernel.

``tests_abaqus``
   Model construction, constraints, analysis, and homogenization tests that
   require Abaqus.

The test runners use the Gmsh dependency installed in the active Python
environment. Before launching Abaqus, activate the isolated runtime generated
by ``abqhom-install-abaqus`` so that Abaqus imports the matching abqhom, Gmsh
and Paraqus packages.

Run the common tests
--------------------

From a Python interpreter that can import ``abqhom``:

.. code-block:: python

   import abqhom.tests as tests

   tests.run_python_tests()

Run the complete suite
----------------------

Execute the following through an Abaqus CAE Python kernel:

.. code-block:: python

   import abqhom.tests as tests

   tests.run_abaqus_tests()

Both runners isolate generated files in a temporary directory and write a
test report to the original working directory. Artifacts from failing tests
are retained for diagnosis.

Documentation checks
--------------------

Build the documentation with warnings treated as errors:

.. code-block:: console

   python -m sphinx -W --keep-going -b html doc/source doc/build/html
