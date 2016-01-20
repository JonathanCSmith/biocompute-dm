import sys

sys.path.insert(0, "./")

from biocomputedm.biocomputedm import create_app

__author__ = "jon"

# Build the app last
application = create_app()
