import os

from biocomputedm.database import SurrogatePK, Model, reference_col, relationship
from biocomputedm.extensions import db
from biocomputedm.admin.helpers import helper_functions
from flask import current_app
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


# Permissions wrapper & environment contextualiser
class Group(SurrogatePK, Model):
    name = db.Column(db.String(50), unique=True, nullable=False)
    directory = db.Column(db.String(120), nullable=False)

    member = relationship("Person", backref="group", lazy="dynamic")

    # TODO: CHANGE
    investigation = db.RelationshipProperty("Investigation", backref="group", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="group", lazy="dynamic")
    submission = db.RelationshipProperty("Submission", backref="group", lazy="dynamic")
    sample_group = db.RelationshipProperty("SampleGroup", backref="group", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="group", lazy="dynamic")

    __tablename__ = "Group"

    def __init__(self, name):
        db.Model.__init__(self, name=name)
        self.directory = os.path.join(current_app.config["RAW_DATA_PATH_ON_WEBSERVER"], self.display_key)
        helper_functions.create_group_directory(self)

    def __repr__(self):
        return "<Group %r>" % self.name

    def set_administrator(self, administrator):
        # Programming error!
        if administrator.get_role == "Site Admin":
            return

        administrator.set_role("Group Admin")
        self.member.append(administrator)


def create_group(group_name, admin_name, admin_password, admin_email):
    group = Group.create(name=group_name)
    user = User.create(name=admin_name, email=admin_email, password=admin_password, group=group)
    group.set_administrator(user)
    return group


# Person table, abstract parent for individuals interacting with the software
class Person(UserMixin, SurrogatePK, Model):
    username = db.Column(db.String(50), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
    role = db.Column(db.Enum("Member", "Group Admin", "Site Admin"), default="Member")
    type = db.Column(db.String(50), nullable=False)

    group_id = reference_col("Group")

    # TODO: CHANGE
    investigation = db.RelationshipProperty("Investigation", backref="submitter", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="submitter", lazy="dynamic")

    _password = db.Column(db.String(160))

    __tablename__ = "Person"
    __mapper_args__ = {"polymorphic_on": type}

    def __init__(self, username, email, password, group):
        db.Model.__init__(self, username=username, email=email)
        self.set_password(password)
        group.member.append(self)
        group.save()

    def __repr__(self):
        return "<Member %r %r>" % (self.username, self.email)

    def get_id(self):
        return self.display_key

    def set_password(self, password):
        self._password = generate_password_hash(password)

    def check_password(self, pwd):
        return check_password_hash(self._password, pwd)

    def get_role(self):
        return self.role

    def set_role(self, role):
        self.role = role


# User - someone who can submit jobs to the software
class User(Person):
    id = reference_col("Person", primary_key=True)

    # TODO: Change
    submission = db.RelationshipProperty("Submission", backref="submitter", lazy="dynamic")
    sample_group = db.RelationshipProperty("SampleGroup", backref="submitter", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="submitter", lazy="dynamic")

    __tablename__ = "User"
    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        Person.__init__(self, username=username, email=email, password=password, group=group)
        helper_functions.create_user_directory(self, password)

    def __repr__(self):
        return "<User %r %r>" % (self.login_name, self.email)


# Customer table
class Customer(Person):
    id = reference_col("Person", primary_key=True)

    __tablename__ = "Customer"
    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        Person.__init__(self, username=username, email=email, password=password, group=group)

    def __repr__(self):
        return "<Customer %r %r>" % (self.login_name, self.email)
