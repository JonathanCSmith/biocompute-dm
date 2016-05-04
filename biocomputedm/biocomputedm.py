from biocomputedm.admin.models import User, Person, UserGroup
from biocomputedm.admin.models import refresh_reference_data_library
from biocomputedm.pipelines.models import refresh_pipelines
from flask import Flask
from flask.ext.bootstrap import Bootstrap
from flask.ext.login import current_user
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
    login_manager.login_message_category = "warning"

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
        return Person.query.filter_by(display_key=user_id).first()


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
    from flask import render_template, redirect, url_for

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
        if current_user is not None and current_user.is_authenticated:
            return redirect(url_for("content.activity"))

        else:
            return redirect(url_for("content.welcome"))

    @app.route("/empty")
    def empty():
        return render_template("empty.html", title="Alice has fallen down the rabbit hole!")


def load_defaults(app):
    # Build our default user
    admin = User.query.filter_by(username=app.config["SITE_ADMIN_USERNAME"]).first()

    if admin is None:
        group = UserGroup.create(
            group_name="Site Admins",
            admin_name=app.config["SITE_ADMIN_USERNAME"],
            admin_password=app.config["SITE_ADMIN_PASSWORD"],
            admin_email=app.config["SITE_ADMIN_EMAIL"]
        )
        admin = group.members.first()
        admin.set_role("Site Admin")
        admin.save()

    # Build our default pipelines
    refresh_pipelines()

    # Build our reference data library
    refresh_reference_data_library()


def get_registered_users(app):
    persons = User.query.all()
    names = []
    for person in persons:
        names.append(person.username)

    return names
