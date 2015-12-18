__author__ = 'jon'

# Application setup
from flask import Flask, g

app = Flask(__name__)
app.config.from_object("config")

import config

app.secret_key = config.SECRET_KEY

# Database setup
from flask.ext.sqlalchemy import SQLAlchemy

db = SQLAlchemy(app)

# Flask-login setup
from flask.ext.login import LoginManager, current_user

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"
login_manager.login_message_category = "warning"

# Flask bootstrap setup
from flask_bootstrap import Bootstrap

Bootstrap(app)


@app.before_request
def before_request():
    g.user = current_user


@login_manager.user_loader
def load_user(user_id):
    from app.models import User
    return User.query.filter_by(id=int(user_id)).first()


# Function to allow role management on pages
def login_required(*role):
    def wrapper(fn):
        from functools import wraps

        @wraps(fn)
        def decorated_view(*args, **kwargs):
            if current_user is not None and current_user.is_authenticated:
                if current_user.get_role() == "Site Admin":
                    return fn(*args, **kwargs)

                elif role == "ANY":
                    return fn(*args, **kwargs)

                else:
                    for item in role:
                        if item == "ANY":
                            return fn(*args, **kwargs)

                        if current_user.get_role() == item:
                            return fn(*args, **kwargs)

            return login_manager.unauthorized()

        return decorated_view

    return wrapper


# Setup the actual website
from app import views, models
