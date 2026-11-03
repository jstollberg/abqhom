Development
===========

The `contribution guide <https://github.com/jstollberg/abqhom/blob/main/CONTRIBUTING.md>`_
describes issue reporting, development setup, code conventions, testing, and
pull-request expectations.

Repository layout
-----------------

.. code-block:: text

   src/abqhom/
       gmshmodel.py    geometry, boundary detection, and periodic meshing
       abqmodel.py     Abaqus model, boundary conditions, and averaging
       utils.py        tensor conversion, bounds, residuals, and file helpers
       constants.py    boundary topology and element-type maps
       tests/          common and Abaqus-dependent tests
   examples/           concise tutorial analyses
   doc/                Sphinx and Read the Docs sources

Building the documentation
---------------------------

Install the documentation dependencies in a regular Python environment:

.. code-block:: console

   python -m pip install -r doc/requirements.txt

Then build from the repository root:

.. code-block:: console

   python -m sphinx -W --keep-going -b html doc/source doc/build/html

Alternatively, run ``make html`` or ``.\make.bat html`` from ``doc``. The
HTML output is written to ``doc/build/html``. Autodoc mocks the external
solver and meshing modules, so the documentation build does not require an
Abaqus, Gmsh, or Paraqus installation.

Read the Docs
-------------

The repository-level ``.readthedocs.yaml`` selects the Python version,
installs ``doc/requirements.txt``, and points Sphinx to
``doc/source/conf.py``. Keep the local warning-clean build above passing
before publishing changes.

Code and documentation conventions
----------------------------------

Public classes and functions use NumPy-style docstrings. New public behavior
should include its assumptions, input and output conventions, and underlying
mathematical quantity where that prevents ambiguity. Examples should remain
self-contained and use the same terminology as the API: *primary* and
*secondary* periodic entities, *displacement gradient*, *effective stress*,
and *incremental work*.

License
-------

``abqhom`` is distributed under the GNU General Public License v3.0 or later.
It uses Gmsh for meshing and Paraqus for VTK export.
