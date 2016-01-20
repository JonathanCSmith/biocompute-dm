from biocomputedm.admin.models import User, create_group, Group
from biocomputedm.decorators import login_required
from biocomputedm.pipelines.models import Pipeline, create_module, create_option
from biocomputedm.pipelines.models import create_pipeline
from flask import Flask
from flask.ext.bootstrap import Bootstrap
from .extensions import db, mail, login_manager

__author__ = "jon"


# TODO: Mark out dead users/groups and rm from folders


def create_app():
    app = Flask(__name__)
    app.config.from_object("config")
    app.secret_key = app.config["SECRET_KEY"]

    setup_logging(app)
    setup_extensions(app)
    setup_hooks(app)
    setup_modules(app)
    setup_default_routes(app)

    return app


def setup_logging(app):
    # Don't create a file or send us mail when we are testing
    if app.debug:
        return

    import logging
    from logging.handlers import RotatingFileHandler, SMTPHandler

    # Set the overall logging level
    app.logger.setLevel(logging.INFO)

    # Create a file specific log handler
    log = RotatingFileHandler(app.config["LOG_LOCATION"], maxBytes=10000, backupCount=10)
    log.setLevel(logging.INFO)
    log.setFormatter(
            logging.Formatter(
                    "%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]"
            )
    )
    app.logger.addHandler(log)

    # Create a mail specific log handler
    mail_log = SMTPHandler(
            app.config["MAIL_SERVER"],
            app.config["MAIL_USERNAME"],
            app.config["MAIL_SOURCE_ADDRESS"],
            '0_ops... Biocompute-DM failed!',
            (app.config["MAIL_USERNAME"], app.config["MAIL_PASSWORD"])
    )
    mail_log.setLevel(logging.ERROR)
    mail_log.setFormatter(
            logging.Formatter(
                    "%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]"
            )
    )
    app.logger.addHandler(mail_log)


def setup_extensions(app):
    # flask-sqlalchemy
    db.init_app(app)

    # flask-mail
    mail.init_app(app)

    # flask-login
    login_manager.init_app(app)
    login_manager.login_view = "admin.login"
    login_manager.refresh_view = "admin.reauth"

    # flask-bootstrap
    Bootstrap(app)


def setup_hooks(app):
    from flask.ext.login import current_user

    @app.before_request
    def before_request():
        from flask import g
        g.user = current_user

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.filter_by(display_key=user_id).first()


def setup_modules(app):
    from biocomputedm.content.views import content
    from biocomputedm.admin.views import admin
    from biocomputedm.manage.views import manage
    from biocomputedm.pipelines.views import pipelines

    app.register_blueprint(content)
    app.register_blueprint(admin)
    app.register_blueprint(manage)
    app.register_blueprint(pipelines)


def setup_default_routes(app):
    from flask import render_template, redirect, url_for, g

    @app.errorhandler(403)
    def forbidden_page(error):
        return render_template("forbidden_page.html"), 403

    @app.errorhandler(404)
    def page_not_found(error):
        return render_template("page_not_found.html"), 404

    @app.errorhandler(500)
    def server_error(error):
        return render_template("server_error.html"), 500

    @app.route("/")
    @app.route("/index")
    def index():
        if g.user is not None and g.user.is_authenticated:
            return redirect(url_for("content.activity"))

        else:
            return redirect(url_for("content.welcome"))

    @app.route("/empty")
    def empty():
        return render_template("empty.html", title="Alice has fallen down the rabbit hole!")


def load_defaults(app):
    import os
    import json

    # Build our default user
    admin = User.query.filter_by(username=app.config["SITE_ADMIN_USERNAME"]).first()

    if admin is None:
        group = Group.create(name="Site Admins")
        user = User.create(username=app.config["SITE_ADMIN_USERNAME"], password=app.config["SITE_ADMIN_PASSWORD"],
                           email=app.config["SITE_ADMIN_EMAIL"], group=group)
        group.set_administrator(user)
        group.save()
        user.set_role("Site Admin")
        user.save()

    # Build our default pipelines
    path = app.config["PIPELINES_PATH_ON_WEBSERVER"]
    directories = os.listdir(path)
    for directory in directories:
        directory_path = os.path.join(path, directory)
        if not os.path.isdir(directory_path):
            continue

        file = os.path.join(directory_path, directory + ".json")
        if not os.path.isfile(file):
            continue

        from biocomputedm.pipelines.helpers import pipeline_mappings_template as template_helper

        if not template_helper.validate(file):
            continue

        json_instance = json.loads(open(file).read())
        name = json_instance.get("name")
        description = json_instance.get("description")
        author = json_instance.get("author")
        version = json_instance.get("version")
        type = json_instance.get("pipeline_type")

        pipeline = Pipeline.query.filter_by(name=name, description=description, author=author, version=version,
                                            type=type).first()
        if pipeline is not None:
            continue

        pipeline = create_pipeline(name, description, author, version, type)

        for module in json_instance.get("modules"):
            mod = create_module(module.get("name"), module.get("description"), module.get("executor"),
                                module.get("index_in_execution_order"), pipeline)

            for option in module.get("options"):
                opt = create_option(option.get("display_name"), option.get("parameter_name"),
                                    option.get("default_value"), option.get("user_interaction_type"), mod)
