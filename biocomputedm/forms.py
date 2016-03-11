__author__ = 'jon'

from flask import flash
from flask.ext.wtf import Form
from flask.ext.wtf.file import FileField, FileRequired
from wtforms import StringField, PasswordField, SubmitField, DateField, IntegerField, BooleanField
from wtforms.validators import DataRequired, Email












class AddDocumentForm(Form):
    description = StringField("Document Desription", validators=[DataRequired()])
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")


class UploadSequencingSubmissionForm(Form):
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")


class NewSequencingSubmissionForm(Form):
    sequence_submission_name = StringField("Sequence Submission Name", validators=[DataRequired()])
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


class UploadSequencingSampleMappingsForm(Form):
    file_upload = FileField("Upload File", validators=[FileRequired()])
    submit = SubmitField("Upload")
