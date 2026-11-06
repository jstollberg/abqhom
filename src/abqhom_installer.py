# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Bootstrap the Abaqus installer from the adjacent abqhom package."""

import os
import sys


def main():
    """Run the installer from the same distribution as this entry point."""

    # Console-script wrappers can inherit a working-directory checkout ahead
    # of site-packages. Prefer the package adjacent to this installed wrapper
    # so the installer and reported version always belong to one distribution.
    package_root = os.path.dirname(os.path.abspath(__file__))
    if not sys.path or sys.path[0] != package_root:
        sys.path.insert(0, package_root)

    from abqhom.install_abaqus import main as install_main
    return install_main()


if __name__ == "__main__":
    sys.exit(main())
