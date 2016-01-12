from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.mail import Mail
from flask.ext.login import LoginManager

__author__ = "jon"

db = SQLAlchemy()
mail = Mail()
login_manager = LoginManager()
