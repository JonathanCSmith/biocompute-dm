__author__ = 'jon'

from flask.ext.wtf import Form
from wtforms import StringField, PasswordField, SubmitField
from wtforms.validators import DataRequired, Email
from app.models import User


class RegisterForm(Form):
    username = StringField("Username", validators=[DataRequired()])
    password = PasswordField("Password", validators=[DataRequired()])
    email = StringField("Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            return False

        user = User.query.filter_by(email=str(self.email.data.lower())).first()
        if user is None:
            user = User.query.filter_by(username=str(self.username.data.lower())).first()

        if user:
            self.email.errors.append("That email is already taken.")
            return False

        else:
            return True


class LoginForm(Form):
    username = StringField("Username", validators=[DataRequired()])
    password = PasswordField("Password", validators=[DataRequired()])
    submit = SubmitField("Login")

    def validate(self):
        if not Form.validate(self):
            return False

        user = User.query.filter_by(username=str(self.username.data)).first()
        check = user.check_password(str(self.password.data))
        if user is not None and check:
            return True

        else:
            self.password.errors.append("The username or password supplied was incorrect")
            return False
