#!/usr/bin/env python3
"""Build script for creating a standalone SAMANTHEX executable.

This script uses PyInstaller to create a single-file Windows executable
that includes Python, all dependencies, and the SAMANTHEX application.

Usage:
    python build_exe.py

Requirements:
    pip install pyinstaller

The executable will be created in the 'dist/' directory.
"""

import os
import subprocess
import sys


def check_pyinstaller():
    """Check if PyInstaller is installed."""
    try:
        import PyInstaller
        print(f"[OK] PyInstaller {PyInstaller.__version__} found")
        return True
    except ImportError:
        print("[!] PyInstaller not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
        print("[OK] PyInstaller installed")
        return True


def build():
    """Build the standalone executable."""
    print("=" * 60)
    print("  SAMANTHEX — Building Standalone Executable")
    print("=" * 60)
    print()

    # Check dependencies
    check_pyinstaller()

    # Determine paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_script = os.path.join(script_dir, "run.py")

    if not os.path.exists(run_script):
        print(f"[ERROR] run.py not found at {run_script}")
        sys.exit(1)

    # Build PyInstaller command
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onefile",
        "--windowed",
        "--name", "Samanthex",
        # Hidden imports that PyInstaller might miss
        "--hidden-import", "Bio",
        "--hidden-import", "Bio.SeqIO",
        "--hidden-import", "Bio.Seq",
        "--hidden-import", "Bio.Align",
        "--hidden-import", "Bio.Blast",
        "--hidden-import", "Bio.Blast.NCBIWWW",
        "--hidden-import", "Bio.Blast.NCBIXML",
        "--hidden-import", "numpy",
        # Collect all BioPython data
        "--collect-all", "Bio",
        # Working directory
        "--distpath", os.path.join(script_dir, "dist"),
        "--workpath", os.path.join(script_dir, "build"),
        "--specpath", script_dir,
        # Clean build
        "--clean",
        # Entry point
        run_script,
    ]

    # Add icon if available
    icon_path = os.path.join(script_dir, "icon.ico")
    if os.path.exists(icon_path):
        cmd.extend(["--icon", icon_path])
        print(f"[OK] Using icon: {icon_path}")

    print(f"\n[...] Running PyInstaller...\n")
    print(f"  Command: {' '.join(cmd)}\n")

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Build failed with exit code {e.returncode}")
        sys.exit(1)

    # Verify output
    exe_path = os.path.join(script_dir, "dist", "Samanthex.exe")
    if os.path.exists(exe_path):
        size_mb = os.path.getsize(exe_path) / (1024 * 1024)
        print(f"\n{'=' * 60}")
        print(f"  BUILD SUCCESSFUL!")
        print(f"  Executable: {exe_path}")
        print(f"  Size: {size_mb:.1f} MB")
        print(f"{'=' * 60}")
    else:
        print(f"\n[WARNING] Expected executable not found at {exe_path}")
        print("  Check the dist/ directory for the output.")


if __name__ == "__main__":
    build()
