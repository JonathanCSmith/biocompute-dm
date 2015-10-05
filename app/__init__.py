__author__ = 'jon'

# Application setup
from flask import Flask, g

app = Flask(__name__)
app.config.from_object("config")

from config import SECRET_KEY

app.secret_key = SECRET_KEY

# Database setup
from flask.ext.sqlalchemy import SQLAlchemy

db = SQLAlchemy(app)

# Flask-login setup
from flask.ext.login import LoginManager, current_user

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"
login_manager.login_message_category = "warning"


@app.before_request
def before_request():
    g.user = current_user


@login_manager.user_loader
def load_user(user_id):
    from app.models import User
    return User.query.filter_by(id=int(user_id)).first()


# Function to allow role management on pages
def login_required(role="ANY"):
    def wrapper(fn):
        from functools import wraps

        @wraps(fn)
        def decorated_view(*args, **kwargs):
            if current_user is None:
                return login_manager.unauthorized()

            if not current_user.is_authenticated:
                return login_manager.unauthorized()

            if (current_user.get_user_role() != role) and (role != "ANY"):
                return login_manager.unauthorized()

            return fn(*args, **kwargs)

        return decorated_view

    return wrapper


# Setup the actual website
from app import views, models
