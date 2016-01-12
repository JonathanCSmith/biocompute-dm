import sys
import config

from biocomputedm.biocomputedm import create_app

__author__ = "jon"

sys.path.insert(0, config.ROOT_WEBSERVER_DIRECTORY)

# Build the app last
application = create_app()
