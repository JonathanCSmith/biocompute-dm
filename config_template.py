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

HPC_JOB_SUBMISSION_FILE = ""
SFTP_SCRIPTS_PATH = ""

PIPELINES_PATH_ON_WEBSERVER = ""
PIPELINES_PATH_ON_HPC = ""

RAW_DATA_PATH_ON_WEBSERVER = ""
RAW_DATA_PATH_ON_HPC = ""

SUBMISSIONS_PATH_ON_WEBSERVER = ""
SUBMISSIONS_PATH_ON_HPC = ""

SAMPLES_PATH_ON_WEBSERVER = ""
SAMPLES_PATH_ON_HPC = ""

INVESTIGATIONS_PATH_ON_WEBSERVER = ""
INVESTIGATIONS_PATH_ON_HPC = ""
