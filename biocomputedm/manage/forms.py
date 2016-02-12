from flask import flash
from flask.ext.wtf import Form
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired


class NewSubmissionForm(Form):
    submission_name = StringField("Submission Name", validators=[DataRequired()])
    submission_description = StringField("Submission Description", validators=[DataRequired()])
    submit = SubmitField("Submit")

