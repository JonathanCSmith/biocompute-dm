__author__ = 'jon'

from flask import flash
from flask.ext.wtf import Form
from flask.ext.wtf.file import FileField, FileRequired
from wtforms import StringField, PasswordField, SubmitField, DateField, IntegerField, BooleanField
from wtforms.validators import DataRequired, Email


class LoginForm(Form):
    username = StringField("Username", validators=[DataRequired()])
    password = PasswordField("Password", validators=[DataRequired()])
    submit = SubmitField("Login")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        from app.models import User
        user = User.query.filter_by(login_name=str(self.username.data)).first()
        if user:
            if user.check_password(str(self.password.data)):
                return True
            else:
                return False

        else:
            flash("There was an error in your submission", "error")
            self.password.errors.append("The username or password supplied was incorrect")
            return False


class NewInvestigationForm(Form):
    investigation_name = StringField("Investigation Name", validators=[DataRequired()])
    investigation_lead = StringField("Investigation Lead", validators=[DataRequired()])
    submit = SubmitField("Submit")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        return True


class CreateGroupForm(Form):
    group_name = StringField("Group Name", validators=[DataRequired()])
    admin_login = StringField("Group Administrator Name", validators=[DataRequired()])
    admin_password = PasswordField("Group Administrator Password", validators=[DataRequired()])
    admin_email = StringField("Group Administrator Email", validators=[DataRequired(), Email()])
    submit = SubmitField("Create Account")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        from app.models import User, Group
        user = User.query.filter_by(email=str(self.admin_email.data)).first()
        if user is None:
            user = User.query.filter_by(login_name=str(self.admin_login.data)).first()

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

        from app.models import User
        user = User.query.filter_by(email=str(self.login_email.data)).first()
        if user is None:
            user = User.query.filter_by(login_name=str(self.login_name.data)).first()

        if user:
            flash("This username already exists", "error")
            return False

        return True


class AddDocumentForm(Form):
    description = StringField("Document Desription", validators=[DataRequired()])
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")


class UploadSequencingRunForm(Form):
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")


class NewSequencingProjectForm(Form):
    sequence_run_name = StringField("Sequence Run Name", validators=[DataRequired()])
    flow_cell_id = StringField("Flow Cell ID", validators=[DataRequired()])
    from datetime import datetime
    start_date = DateField("Start Date", format="%Y-%m-%d", default=datetime.today, validators=[DataRequired()])
    completion_date = DateField("Completion Date", format="%Y-%m-%d", default=datetime.today,
                                validators=[DataRequired()])
    genomics_lead = StringField("Genomics Lead", validators=[DataRequired()])
    data_location = StringField("Data Location", validators=[DataRequired()])
    index_tag_cycles = IntegerField("Number of Cycles for Index Tag 1", validators=[DataRequired()])
    index_tag_cycles_2 = IntegerField("Number of Cycles for Index Tag 2", validators=[DataRequired()])
    read_cycles = IntegerField("Number of Cycles for Read 1", validators=[DataRequired()])
    read_cycles_2 = IntegerField("Number of Cycles for Read 2", validators=[DataRequired()])
    paired_end = BooleanField("Paired End")
    submit = SubmitField("Submit")


class UploadSequencingProjectForm(Form):
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")
