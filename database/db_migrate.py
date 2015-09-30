#!flask/bin/python
__author__ = 'jon'

import imp

from migrate.versioning import api
from app import db
from config import SQLALCHEMY_DATABASE_URI
from config import SQLALCHEMY_MIGRATE_REPO

# Create model
print("Creating db model")
v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
migration = SQLALCHEMY_MIGRATE_REPO + ("/versions/%03d_migration.py" % (v+1))
tmp_module = imp.new_module("old_model")
old_model = api.create_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
exec(old_model, tmp_module.__dict__)

# Create migration script
print("Creating migration script")
script = api.make_update_script_for_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, tmp_module.meta, db.metadata)
file = open(migration, "wt")
chomp = script.split("\n", 1)
file.write(chomp[0] + "\n")
file.write("from sqlalchemy.dialects.mysql import *\n")
file.write(chomp[1])
file.close()

# Migrating
print("Migrating")
api.upgrade(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
print("New migration saved as " + migration)
print("Current database version: " + str(v))
