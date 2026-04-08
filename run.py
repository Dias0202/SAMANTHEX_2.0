#!/usr/bin/env python3
"""SAMANTHEX - Sanger Consensus Assembler.

Launch script for the SAMANTHEX application.
Run this file to start the GUI:

    python run.py
"""

import sys
import os

# Ensure project root is in path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from samanthex.main import main

if __name__ == "__main__":
    main()
