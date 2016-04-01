from biocomputedm.pipelines.models import Pipeline
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

    choices = [("NO", "None")]
    pipelines = Pipeline.query.all()
    for pipeline in pipelines:
        choices.append((pipeline.display_key, pipeline.name + " (" + pipeline.version + ")"))

    pipeline = SelectField("Pipeline Source", default=choices[0], choices=choices, validators=[DataRequired()])
    submit_field = SubmitField("Submit")


class UpdateSampleGroupForm(Form):
    submit = SubmitField("Submit")


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
