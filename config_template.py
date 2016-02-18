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

LOCAL_WEBSERVER_PORT = "5000"  # Flask default
NETWORK_PATH_TO_WEBSERVER_FROM_HPC = ""

# NOTE Everything below here should be OS independent paths (i.e. //)
HPC_JOB_SUBMISSION_FILE = ""
SCRIPTS_PATH = ""

# SFTP root folder
SFTP_USER_ROOT_PATH = ""

# Allows different relative path to be provided per environment - however, the contents should be identical (i.e. remote
# mount, sshfs, samba etc - they can also be exactly the same path if the webserver is executed from the HPC environ)
HPC_ROOT_PATH = ""
WEBSERVER_ROOT_PATH = ""

# Paths relative to the above
PIPELINE_SCRIPTS_PATH_AFTER_RELATIVE_ROOT = ""
PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT = ""
SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT = ""
SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT = ""
PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT = ""
