from biocomputedm.admin.models import User, Group, Person
from flask import flash
from flask.ext.wtf import Form
from wtforms import PasswordField, SubmitField
from wtforms import StringField
from wtforms.validators import DataRequired, Email


class LoginForm(Form):
    username = StringField("Username", validators=[DataRequired()])
    password = PasswordField("Password", validators=[DataRequired()])
    submit = SubmitField("Login")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        user = Person.query.filter_by(username=str(self.username.data)).first()
        if user:
            if user.check_password(str(self.password.data)):
                return True
            else:
                return False

        else:
            self.password.errors.append("The username or password supplied was incorrect")
            return False


class CreateGroupForm(Form):
    group_name = StringField("Group Name", validators=[DataRequired()])
    admin_login = StringField("Group Administrator Name", validators=[DataRequired()])
    admin_password = PasswordField("Group Administrator Password", validators=[DataRequired()])
    admin_email = StringField("Group Administrator Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            return False

        user = User.query.filter_by(email=str(self.admin_email.data)).first()
        if user is None:
            user = User.query.filter_by(username=str(self.admin_login.data)).first()

        if user:
            flash("This username already exists", "error")
            return False

        group = Group.query.filter_by(name=str(self.group_name.data)).first()
        if group:
            flash("This group name already exists", "error")
            return False

        return True


class CreatePerson(Form):
    login_name = StringField("Login Name", validators=[DataRequired()])
    login_password = PasswordField("Password", validators=[DataRequired()])
    login_email = StringField("Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        user = User.query.filter_by(email=str(self.login_email.data)).first()
        if user is None:
            user = User.query.filter_by(username=str(self.login_name.data)).first()

        if user:
            flash("This username already exists", "error")
            return False

        return True


class UpdateCustomerLinkForm(Form):
    submit = SubmitField("Submit")
