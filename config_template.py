__author__ = 'jon'

SECRET_KEY = ""

import os
basedir = os.path.abspath(os.path.dirname(__file__))

USERNAME = ""
PASSWORD = ""
LOCATION = ""
DATABASE = ""

SQLALCHEMY_DATABASE_URI = "mysql+pymysql://" + USERNAME + ":" + PASSWORD + "@" + LOCATION + "/" + DATABASE
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, "db_repository")

SITE_ADMIN_USERNAME = ""
SITE_ADMIN_PASSWORD = ""
SITE_ADMIN_EMAIL = ""

REMOTE_PIPELINES_PATH = ""
REMOTE_RAW_PATH = ""
REMOTE_SUBMISSIONS_PATH = ""
REMOTE_SAMPLES_PATH = ""
REMOTE_INVESTIGATIONS_PATH = ""
