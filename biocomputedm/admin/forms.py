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


class ChangePasswordForm(Form):
    old_password = PasswordField("Original Password", validators=[DataRequired()])
    new_password = PasswordField("New Password", validators=[DataRequired()])
    new_password_2 = PasswordField("Repeat New Password", validators=[DataRequired()])
    submit = SubmitField("Submit")

    def validate(self):
        if not Form.validate(self):
            return False;

        if str(self.new_password.data) != str(self.new_password_2.data):
            flash("Your password submissions did not match", "warning")
            return False

        return True


class CreateGroupForm(Form):
    group_name = StringField("Group Name", validators=[DataRequired()])
    admin_login = StringField("Administrator Name", validators=[DataRequired()])
    #admin_password = PasswordField("Group Administrator Password", validators=[DataRequired()])
    admin_email = StringField("Administrator Email", validators=[DataRequired(), Email()])
    admin_email_2 = StringField("Repeat Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            return False

        trimmed_name = ''.join(str(self.group_name.data).split())
        if trimmed_name != str(self.group_name.data):
            flash("Group names cannot contain any whitespace characters", "warning")
            return False

        user = User.query.filter_by(email=str(self.admin_email.data)).first()
        if user is None:
            user = User.query.filter_by(username=str(self.admin_login.data)).first()

        if user:
            flash("This username already exists", "warning")
            return False

        group = Group.query.filter_by(name=str(self.group_name.data)).first()
        if group:
            flash("This group name already exists", "warning")
            return False

        if str(self.admin_email.data) != str(self.admin_email_2.data):
            flash("The emails provided were not the same", "warning")
            return False

        return True


class CreatePerson(Form):
    login_name = StringField("Login Name", validators=[DataRequired()])
    # login_password = PasswordField("Password", validators=[DataRequired()])
    login_email = StringField("Email", validators=[DataRequired(), Email()])
    login_email_2 = StringField("Repeat Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "warning")
            return False

        trimmed_name = ''.join(str(self.login_name.data).split())
        if trimmed_name != str(self.login_name.data):
            flash("Usernames cannot contain any whitespace characters", "warning")
            return False

        user = User.query.filter_by(email=str(self.login_email.data)).first()
        if user is None:
            user = User.query.filter_by(username=str(self.login_name.data)).first()

        if user:
            flash("This username already exists", "warning")
            return False

        if str(self.login_email.data) != str(self.login_email_2.data):
            flash("The emails provided were not the same", "warning")
            return False

        return True


class UpdateCustomerLinkForm(Form):
    submit = SubmitField("Submit")
