from flask.ext.migrate import Migrate
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.mail import Mail
from flask.ext.login import LoginManager

__author__ = "jon"

db = SQLAlchemy()
migrate = Migrate()

mail = Mail()
login_manager = LoginManager()
