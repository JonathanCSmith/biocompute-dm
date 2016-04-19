from flask import flash
from flask.ext.wtf import Form
from wtforms import StringField, SubmitField, SelectField, FileField, TextAreaField
from wtforms.validators import DataRequired


class NewSubmissionForm(Form):
    submission_name = StringField("Submission Name", validators=[DataRequired()])
    submission_description = TextAreaField("Submission Description", validators=[DataRequired()])
    submit = SubmitField("Submit")


class NewSampleGroupForm(Form):
    name = StringField("Sample Group Name", validators=[DataRequired()])
    submit_field = SubmitField("Submit")


class UpdateSampleGroupForm(Form):
    source_type = SelectField("Data Source (Pipeline)")
    submit = SubmitField("Submit")

    def fill(self, source_pipelines):
        self.source_type.choices=[(s, s) for s in source_pipelines]
        return


class NewProjectForm(Form):
    investigation_name = StringField("Project Name", validators=[DataRequired()])
    investigation_description = StringField("Project Description", validators=[DataRequired()])
    submit = SubmitField("Submit")

    def validate(self):
        if not Form.validate(self):
            flash("There was an error in your submission", "error")
            return False

        return True


class AddDocumentForm(Form):
    description = StringField("Document Desription", validators=[DataRequired()])
    file_upload = FileField("Upload File")
    submit = SubmitField("Upload")


class UpdateProjectForm(Form):
    submit = SubmitField("Submit")
