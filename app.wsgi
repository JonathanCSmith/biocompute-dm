import sys

sys.path.insert(0, "/var/www/biocomputedm")

from biocomputedm import biocomputedm

__author__ = "jon"

# Build the app last
application = biocomputedm.create_app()
