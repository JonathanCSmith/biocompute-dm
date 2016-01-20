import os
import sys

sys.path.insert(0, os.path.realpath(__file__))

from biocomputedm.biocomputedm import create_app

__author__ = "jon"

# Build the app last
application = create_app()
