# Contributing to abqhom

Contributions are welcome. Bug reports, documentation corrections, focused
feature proposals, tests, and example improvements all help the project.

## Before making a change

Search the [issue tracker](https://github.com/jstollberg/abqhom/issues) for an
existing discussion. Open an issue before starting a substantial feature,
public-API change, or change to the supported modeling assumptions. Small bug
fixes and editorial corrections may be submitted directly.

Bug reports should include:

- the abqhom version or commit
- the versions of Abaqus, Gmsh, Paraqus, Python, and the operating system
- the command or script that produced the problem
- the complete error message or traceback
- a minimal reproducible model when possible

Do not include confidential or proprietary input files in a public issue.

## Development setup

Fork and clone the repository, then create a branch for the change. Install
the checkout in editable mode in a Python 3.9 or newer environment:

```console
python -m pip install --editable .
python -m pip install -r doc/requirements.txt
```

The common tests run in this regular Python environment. Changes that require
Abaqus should also be tested through the isolated runtime described in the
[installation guide](https://abqhom.readthedocs.io/en/latest/installation.html).
Create or refresh that runtime from the checkout with:

```console
abqhom-install-abaqus --abqhom-source .
```

## Code conventions

- Follow PEP 8 and the existing module structure.
- Keep package modules, the test suite, and example scripts compatible with
  both legacy Python 2.7-based and current Python 3-based Abaqus releases.
- Use NumPy-style docstrings for public classes, methods, and functions. Keep
  private helper docstrings concise.
- Document assumptions, units, tensor conventions, and mathematical
  quantities when they are not evident from the interface.
- Use *primary* and *secondary* for periodic entities. Follow the established
  terms *displacement gradient*, *effective stress*, and *incremental work*.
- Prefer focused functions and comments that explain mechanical or numerical
  choices rather than restating the code.
- Add the project SPDX copyright and GPL-3.0-or-later identifiers to new source
  and configuration files where appropriate.

Avoid semicolons that join independent clauses in prose. Documentation and
comments should use direct sentences and remain understandable without the
accompanying research publication.

## Tests

Behavior changes should include a focused regression test. Run the common
test suite from a Python interpreter that can import the editable checkout:

```python
import abqhom.tests as tests

tests.run_python_tests()
```

Run the complete suite through the Abaqus CAE Python kernel when a change
affects model construction, constraints, analysis, output databases, or the
homogenization workflow:

```python
import abqhom.tests as tests

tests.run_abaqus_tests()
```

The [testing documentation](https://abqhom.readthedocs.io/en/latest/tests.html)
explains the test groups and generated reports. If an Abaqus test cannot be
run because the required solver or license is unavailable, state that clearly
in the pull request.

## Documentation and examples

Update docstrings and user documentation whenever public behavior,
requirements, or modeling assumptions change. Build the documentation with
warnings treated as errors:

```console
python -m sphinx -W --keep-going -b html doc/source doc/build/html
```

Tutorial examples should remain self-contained, readable, and limited to the
behavior they demonstrate. Add automated regression coverage to the test
suite rather than using examples as tests. Do not commit Abaqus jobs, output
databases, VTK output, generated figures, or other analysis artifacts unless
they are intentional source assets required by an example.

See the
[development documentation](https://abqhom.readthedocs.io/en/latest/development.html)
for the repository layout and documentation build details.

## Pull requests

Keep each pull request focused. Its description should explain:

- the problem and the chosen solution
- any changes to the public interface or modeling assumptions
- the tests and documentation checks that were run
- the Abaqus release used for solver-dependent verification
- any remaining limitations or checks that could not be performed

Before submitting, check the diff for whitespace errors:

```console
git diff --check
```

Review the diff for unrelated changes and generated files. Use clear commit
messages and keep the branch up to date with `main`.

## License

By submitting a contribution, you agree that it may be distributed under the
GNU General Public License v3.0 or later. Only contribute material that you
have the right to submit, and identify any incorporated third-party material
and its license.
