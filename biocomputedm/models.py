import os
import datetime
import uuid

from biocomputedm.extensions import db
from flask import current_app
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash

__author__ = 'jon'


# class Investigation(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)
#
#     name = db.Column(db.String(40), nullable=False)
#     leader = db.Column(db.String(40), nullable=False)
#     directory = db.Column(db.String(500))
#     description = db.Column(db.TEXT)
#     open_date = db.Column(db.Date, nullable=False)
#     last_update = db.Column(db.Date, nullable=False)
#
#     submitter_id = db.Column(db.Integer, db.ForeignKey("Person.id"))
#     group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))
#
#     sample_group = db.RelationshipProperty("SampleGroup", backref="investigation", lazy="dynamic")
#     document = db.RelationshipProperty("Document", backref="investigation", lazy="dynamic")
#
#     __tablename__ = "Investigation"
#
#     def __init__(self, investigation_name, investigation_lead):
#         self.name = investigation_name
#         self.leader = investigation_lead
#
#         today = datetime.date.today()
#         self.open_date = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
#         self.last_update = self.open_date
#
#     def __repr__(self):
#         return "<Investigation %r %r>" % (self.investigation_name, self.investigation_lead)
#
#     def set_last_update(self):
#         today = datetime.date.today()
#         self.last_update = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
#
#     def validate_investigation_directory(self):
#         if self.directory is None:
#             # Local variant as this code is execute webserver side
#             import config
#             path = config.INVESTIGATIONS_PATH_ON_WEBSERVER
#             self.directory = os.path.join(path, str(self.id) + "_" + self.name)
#
#         # Create our investigation directory if necessary
#         try:
#             if not os.path.exists(self.directory):
#                 os.mkdir(self.directory)
#                 os.chmod(self.directory, 0o777)
#
#         except OSError as e:
#             print(e)
#             return None
#
#         # Create our documents directory if necessary
#         docs = os.path.join(self.directory, "documents")
#         try:
#             if not os.path.exists(docs):
#                 os.mkdir(docs)
#                 os.chmod(docs, 0o777)
#
#         except OSError as e:
#             print(e)
#             return None
#
#         return docs
#
#

