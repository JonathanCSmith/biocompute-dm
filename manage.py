import os
import subprocess

from biocomputedm import utils
from biocomputedm.biocomputedm import create_app, load_defaults, get_registered_users
from biocomputedm.extensions import db
from biocomputedm.extensions import migrate as migrator
from flask.ext.migrate import MigrateCommand, init, migrate, upgrade
from flask.ext.script import Manager
from sqlalchemy import MetaData

__author__ = "jon"

app = create_app()
migrator.init_app(app, db)

manager = Manager(app)
manager.add_command("db", MigrateCommand)


@manager.command
def debug():
    app.run(debug=True)


@manager.command
def run():
    app.run()


@manager.command
def initialize_master_database_and_clean_versioning():
    # Cleanup folder
    cleanup_directory = utils.get_path("scripts", "webserver")
    cleanup_directory = os.path.join(cleanup_directory, "cleanup")

    # Wipe any existing
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

    # Migrate setup
    init()
    migrate()
    upgrade()

    # Load the minimum information into the app
    load_defaults(app)


@manager.command
def generate_migration():
    migrate()
    upgrade()


@manager.command
def update():
    upgrade()


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

    con = db.engine.connect()
    trans = con.begin()
    meta = MetaData(bind=db.engine, reflect=True)
    for table in reversed(meta.sorted_tables):
        if table.fullname == "alembic_version":
            continue
        con.execute(table.delete())
    trans.commit()

    # Refresh database
    db.session.close()
    db.session.bind.dispose()

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
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Clean directories
    subprocess.Popen(
            [
                "sudo",
                os.path.join(cleanup_directory, "wipe_directory.sh"),
                "-p=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"])
            ]
    ).wait()

    # Serve the directories
    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    link_path = os.path.dirname(os.path.realpath(__file__))
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" +
                app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" +
                app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" +
                app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    serve_script = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "serve.sh")
    subprocess.Popen(
            [
                "sudo",
                serve_script,
                "-s=" + os.path.join(app.config["WEBSERVER_ROOT_PATH"],
                                     app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"]),
                "-t=" + link_path + "/biocomputedm/static/serve/" +
                app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"].split("/")[-1]
            ]
    ).wait()

    # Load the minimum information into the app
    load_defaults(app)


if __name__ == "__main__":
    manager.run()
