# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Bootstrap the Abaqus installer from the adjacent abqhom package."""

import os
import sys


def main():
    """Run the installer from the same distribution as this entry point."""

    package_root = os.path.dirname(os.path.abspath(__file__))
    if not sys.path or sys.path[0] != package_root:
        sys.path.insert(0, package_root)

    from abqhom.install_abaqus import main as install_main
    return install_main()


if __name__ == "__main__":
    sys.exit(main())
