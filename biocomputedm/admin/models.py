import os

from biocomputedm import utils
from biocomputedm.admin.helpers import admin_helper
from biocomputedm.database import *
from biocomputedm.extensions import db
from flask.ext.login import UserMixin
from sqlalchemy import update
from werkzeug.security import generate_password_hash, check_password_hash

sample_customer_table = Table("SampleToCustomer",
                              Column("sample_id", Integer, ForeignKey("Sample.id")),
                              Column("customer_id", Integer, ForeignKey("Customer.id")))

sample_group_customer_table = Table("SampleGroupToCustomer",
                                    Column("sample_group_id", Integer, ForeignKey("SampleGroup.id")),
                                    Column("customer_id", Integer, ForeignKey("Customer.id")))

project_customer_table = Table("ProjectToCustomer",
                               Column("project_id", Integer, ForeignKey("Project.id")),
                               Column("customer_id", Integer, ForeignKey("Customer.id")))


# Permissions wrapper & environment contextualiser
class Group(SurrogatePK, Model):
    name = Column(String(50), unique=True, nullable=False)

    parent_id = reference_col("Group", nullable=True)

    members = relationship("Person", backref="group", lazy="dynamic")
    submissions = relationship("Submission", lazy="dynamic", backref="group")
    pipeline_instances = relationship("PipelineInstance", lazy="dynamic", backref="pipelines")
    samples = relationship("Sample", lazy="dynamic", backref="group")
    sample_groups = relationship("SampleGroup", lazy="dynamic", backref="group")
    projects = relationship("Project", lazy="dynamic", backref="group")

    __tablename__ = "Group"

    def __init__(self, name):
        db.Model.__init__(self, name=name)

        # Hack, we need the default values here
        db.session.add(self)
        db.session.commit()

        admin_helper.create_group_directory(self)

    def __repr__(self):
        return "<Group %r>" % self.name

    def set_administrator(self, administrator):
        # Programming error!
        if administrator.get_role == "Site Admin":
            return

        administrator.set_role("Group Admin")
        self.members.append(administrator)


def create_group(group_name, admin_name, admin_password, admin_email):
    group = Group.create(name=group_name)
    user = User.create(username=admin_name, email=admin_email, password=admin_password, group=group)
    group.set_administrator(user)
    return group


def create_customer_group(group_name, admin_name, admin_password, admin_email):
    group = Group.create(name=group_name)
    user = Customer.create(username=admin_name, email=admin_email, password=admin_password, group=group)
    group.set_administrator(user)
    return group


# Person table, abstract parent for individuals interacting with the software
class Person(UserMixin, SurrogatePK, Model):
    username = Column(String(50), unique=True, nullable=False)
    email = Column(String(120), unique=True, nullable=False)
    role = Column(db.Enum("Member", "Group Admin", "Site Admin"), default="Member")
    type = Column(String(50), nullable=False)

    group_id = reference_col("Group")

    _password = Column(String(160))

    __tablename__ = "Person"
    __mapper_args__ = {"polymorphic_on": type}

    def __init__(self, username, email, password, group):
        db.Model.__init__(self, username=username, email=email)
        self.set_password(password)
        group.members.append(self)
        group.save()

        # Hack, we need the default values here
        db.session.add(self)
        db.session.commit()

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
        self.update(role=role)


# User - someone who can submit jobs to the software
class User(Person):
    id = reference_col("Person", primary_key=True)

    submissions = relationship("Submission", backref="submitter", lazy="dynamic")

    samples = relationship("Sample", backref="creator", lazy="dynamic")
    sample_groups = relationship("SampleGroup", backref="creator", lazy="dynamic")
    projects = relationship("Project", backref="creator", lazy="dynamic")

    __tablename__ = "User"
    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        Person.__init__(self, username=username, email=email, password=password, group=group)
        admin_helper.create_user_directory(self, password)

    def __repr__(self):
        return "<User %r %r>" % (self.username, self.email)


# Customer table
class Customer(Person):
    id = reference_col("Person", primary_key=True)

    samples = relationship("Sample", secondary=sample_customer_table, lazy="dynamic")
    sample_groups = relationship("SampleGroup", secondary=sample_group_customer_table, lazy="dynamic")
    projects = relationship("Project", secondary=project_customer_table, lazy="dynamic")

    __tablename__ = "Customer"
    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        Person.__init__(self, username=username, email=email, password=password, group=group)

    def __repr__(self):
        return "<Customer %r %r>" % (self.login_name, self.email)


# Reference data set
class ReferenceData(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(db.String(500), nullable=False)
    version = Column(String(50), nullable=False)
    current = Column(Boolean(), default=False)

    __tablename__ = "ReferenceData"
    __table_args__ = (db.UniqueConstraint("name", "description", "version", name="_unique"),)

    def __init__(self, name, description, version):
        db.Model.__init__(self, name=name, description=description, version=version)

    def __repr__(self):
        return "<Reference Data with name: %s, description; %s and version: %s" % (
            self.name, self.description, self.version)


def refresh_reference_data_library():
    # Mark all as legacy on refresh then re-add
    db.session.execute(update(ReferenceData, values={ReferenceData.current: False}))
    db.session.commit()

    # HPC Side as we need the paths to be correct
    path = utils.get_path("reference_data", "webserver")
    directories = os.listdir(path)
    has_new = False
    for directory in directories:
        directory_path = os.path.join(path, directory)
        if not os.path.isdir(directory_path):
            continue

        file = os.path.join(directory_path, directory + ".json")
        if not os.path.isfile(file):
            continue

        from biocomputedm.admin.helpers import resource_helper as template_helper
        if not template_helper.validate(file):
            continue

        has_new |= template_helper.build(file)

    return has_new
