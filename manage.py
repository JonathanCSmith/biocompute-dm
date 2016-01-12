import os

from biocomputedm.biocomputedm import create_app, load_defaults
from coverage.backward import imp
from flask.ext.script import Manager
from migrate.versioning import api

__author__ = "jon"

app = create_app()
manager = Manager(app)

@manager.command
def debug():
    app.run(debug=True)


@manager.command
def run():
    app.run()


@manager.command
def initdb():
    from biocomputedm.extensions import db
    db.drop_all()
    db.create_all()

    # Build a version repository
    if not os.path.exists(app.config["SQLALCHEMY_MIGRATE_REPO"]):
        api.create(app.config["SQLALCHEMY_MIGRATE_REPO"], "database repository")
        api.version_control(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    else:
        api.version_control(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"], api.version(app.config["SQLALCHEMY_MIGRATE_REPO"]))

    # Load the minimum information into the app
    load_defaults(app)


@manager.command
def migrate_db():
    from biocomputedm.extensions import db
    # Create model
    print("Creating db model")
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    migration = app.config["SQLALCHEMY_MIGRATE_REPO"] + ("/versions/%03d_migration.py" % (v + 1))
    tmp_module = imp.new_module("old_model")
    old_model = api.create_model(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    exec(old_model, tmp_module.__dict__)

    # Create migration script
    print("Creating migration script")
    script = api.make_update_script_for_model(
            app.config["SQLALCHEMY_DATABASE_URI"],
            app.config["SQLALCHEMY_MIGRATE_REPO"],
            tmp_module.meta,
            db.metadata
    )
    file = open(migration, "wt")
    chomp = script.split("\n", 1)
    file.write(chomp[0] + "\n")
    file.write("from sqlalchemy.dialects.mysql import *\n")
    file.write(chomp[1])
    file.close()

    # Migrating
    print("Migrating")
    api.upgrade(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    print("New migration saved as " + migration)
    print("Current database version: " + str(v))


@manager.command
def upgrade_db():
    api.upgrade(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    print("Current database version: " + str(v))


@manager.command
def downgrade_db():
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    api.downgrade(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"], v - 1)
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    print("Current database version is: " + str(v))


if __name__ == "__main__":
    manager.run()
