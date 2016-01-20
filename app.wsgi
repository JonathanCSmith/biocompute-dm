import os
import sys

from biocomputedm.biocomputedm import create_app

__author__ = "jon"

sys.path.insert(0, os.path.realpath(__file__))

# Build the app last
application = create_app()
