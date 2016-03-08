import os
import subprocess

from biocomputedm import utils
from biocomputedm.biocomputedm import create_app, load_defaults, get_registered_users
from biocomputedm.extensions import db
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
    # Create model
    print("Creating db model")
    v = api.db_version(app.config["SQLALCHEMY_DATABASE_URI"], app.config["SQLALCHEMY_MIGRATE_REPO"])
    migration = app.config["SQLALCHEMY_MIGRATE_REPO"] + ("/versions/%03d_migration.py" % (v + 1))
    import imp
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


@manager.command
def clean(key):
    if key != "FORCE":
        return

    users = get_registered_users(app)

    # Cleanup folder
    cleanup_directory = utils.get_path("scripts", "webserver")
    cleanup_directory = os.path.join(cleanup_directory, "cleanup")

    # Dump sftp users
    for user in users:
        subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_user.sh"),
                "-u=" + user
            ]
        ).wait()

    # Dump database
    db.session.close()
    db.session.bind.dispose()
    subprocess.Popen(
        [
            "sudo",
            os.path.join(cleanup_directory, "wipe_database.sh"),
            "-u=" + app.config["DATABASE_USERNAME"],
            "-p=" + app.config["DATABASE_PASSWORD"],
            "-l=" + app.config["DATABASE_LOCATION"],
            "-n=" + app.config["DATABASE_NAME"]
        ]
    ).wait()

    # Clean directories
    subprocess.Popen(
        [
            "sudo",
            os.path.join(cleanup_directory, "wipe_directory.sh"),
            "-p=" + app.config["SFTP_USER_ROOT_PATH"]
        ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Serve the directories
    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    link_path = os.path.dirname(os.path.realpath(__file__))
    subprocess.Popen(
        [
            "sudo",
            serve_script,
            "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"]),
            "-t=" + link_path + "/biocomputedm/static/serve/" + app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
        ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" + app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" + app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"], app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" +link_path + "/biocomputedm/static/serve/" + app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    initdb()

if __name__ == "__main__":
    manager.run()
