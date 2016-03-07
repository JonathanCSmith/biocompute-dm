import os

__author__ = "jon"

ROOT_WEBSERVER_DIRECTORY = ""
SECRET_KEY = ""

ERROR_LOG_LOCATION = "error.log"

MAIL_USERNAME = ""
MAIL_PASSWORD = ""
MAIL_SERVER = ""
MAIL_SOURCE_ADDRESS = ""

DATABASE_USERNAME = ""
DATABASE_PASSWORD = ""
DATABASE_LOCATION = ""
DATABASE_NAME = ""

basedir = os.path.abspath(os.path.dirname(__file__))
SQLALCHEMY_DATABASE_URI = "mysql+pymysql://" + DATABASE_USERNAME + ":" + DATABASE_PASSWORD + "@" + DATABASE_LOCATION + "/" + DATABASE_NAME
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, "db_repository")

SITE_ADMIN_USERNAME = ""
SITE_ADMIN_PASSWORD = ""
SITE_ADMIN_EMAIL = ""

# Currently used to assess whether or not we need to fake job submission (i.e. no hpc is present)
HPC_DEBUG = "False"
LOCAL_WEBSERVER_PORT = "5000"  # Flask default
NETWORK_PATH_TO_WEBSERVER_FROM_HPC = ""
NETWORK_PATH_TO_HPC_FROM_WEBSERVER = ""
HPC_USERNAME = ""  # Must have full permissions, i.e. ssh credentials to access from webserver

# ======================== NOTE Everything below here should be OS independent paths (i.e. //) =========================

# Allows different relative path to be provided per environment - however, the contents should be identical (i.e. remote
# mount, sshfs, samba etc - they can also be exactly the same path if the webserver is executed from the HPC environ)
HPC_ROOT_PATH = ""
WEBSERVER_ROOT_PATH = ""
SFTP_USER_ROOT_PATH = ""

# THIS CONFIGURATION IS PROVIDED FOR FLEXIBILITY BUT IT IS ASSUMED THAT THE FILE NAMES & STRUCTURE WILL BE MAINTAINED!
MANAGEMENT_SCRIPTS_PATH_AFTER_RELATIVE_ROOT = "scripts"

# Paths relative to the above
PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT = ""
SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT = ""
SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT = ""
PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT = ""
REFERENCE_DATA_PATH_AFTER_RELATIVE_ROOT = ""
