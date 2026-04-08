"""Entry point for SAMANTHEX application.

Can be run directly or imported by run.py / PyInstaller executable.
"""

import sys
import os


def main():
    """Launch the SAMANTHEX GUI application."""
    # Ensure the package is importable when run from various locations
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    from samanthex.gui.app import SamanthexApp

    app = SamanthexApp()
    app.mainloop()


if __name__ == "__main__":
    main()
