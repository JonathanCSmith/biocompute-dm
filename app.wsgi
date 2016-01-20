import sys
import os

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from biocomputedm import biocomputedm

__author__ = "jon"

# Build the app last
application = biocomputedm.create_app()
