import sys

from biocomputedm import biocomputedm

sys.path.insert(0, "./")

__author__ = "jon"

# Build the app last
application = biocomputedm.create_app()
