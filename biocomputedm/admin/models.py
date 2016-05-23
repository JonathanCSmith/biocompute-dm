import os

from biocomputedm import utils
from biocomputedm.admin.helpers import admin_helper
from biocomputedm.database import *
from biocomputedm.extensions import db, mail
from flask.ext.login import UserMixin
from flask.ext.mail import Message
from sqlalchemy import update
from werkzeug.security import generate_password_hash, check_password_hash

customer_to_sample_mapping = Table("CustomerSampleMapping", Column("customer_id", Integer, ForeignKey("CustomerGroup.id")), Column("sample_id", Integer, ForeignKey("Sample.id")))
customer_to_data_group_mapping = Table("CustomerDataGroupMapping", Column("customer_id", Integer, ForeignKey("CustomerGroup.id")), Column("data_group_id", Integer, ForeignKey("DataGroup.id")))
customer_to_data_item_mapping = Table("CustomerDataItemMapping", Column("customer_id", Integer, ForeignKey("CustomerGroup.id")), Column("data_item_id", Integer, ForeignKey("DataItem.id")))
customer_to_project_mapping = Table("CustomerProjectMapping", Column("customer_id", Integer, ForeignKey("CustomerGroup.id")), Column("project_id", Integer, ForeignKey("Project.id")))


# Permissions wrapper & environment contextualiser
class Group(SurrogatePK, Model):
    name = Column(String(50), unique=True, nullable=False)
    type = Column(String(50), nullable=False)

    members = relationship("Person", backref="group", lazy="dynamic")

    _password = Column(String(160))

    __tablename__ = "Group"
    __mapper_args__ = {"polymorphic_on": type}

    def __init__(self, name, password):
        db.Model.__init__(self, name=name)
        self.set_password(password)

        # Hack, we need the default values here
        db.session.add(self)
        db.session.commit()

        admin_helper.create_group_directory(self, password)

    def __repr__(self):
        return "<Group %r>" % self.name

    def set_password(self, password):
        if self._password is not None:
            admin_helper.change_password("biocompute-DM_group_" + self.name, password)

        self._password = generate_password_hash(password)

    def check_password(self, pwd):
        return check_password_hash(self._password, pwd)

    def set_administrator(self, administrator):
        # Programming error!
        if administrator.get_role == "Site Admin":
            return

        administrator.set_role("Group Admin")
        self.members.append(administrator)


class UserGroup(Group):
    id = reference_col("Group", primary_key=True)

    submissions = relationship("Submission", lazy="dynamic", backref="group")
    pipeline_instances = relationship("PipelineInstance", lazy="dynamic", backref="group")
    samples = relationship("Sample", lazy="dynamic", backref="group")
    data_groups = relationship("DataGroup", lazy="dynamic", backref="group")
    data_items = relationship("DataItem", lazy="dynamic", backref="group")
    projects = relationship("Project", lazy="dynamic", backref="group")

    __tablename__ = "UserGroup"
    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Group.id)}

    def __init__(self, group_name, password, admin_name, admin_password, admin_email):
        if password is None:
            password = uuid.uuid4().hex

        Group.__init__(self, name=group_name, password=password)
        user = User.create(username=admin_name, email=admin_email, password=admin_password, group=self)
        self.set_administrator(user)

        msg = Message(subject="Welcome to Biocompute-DM",
                      body="Your email address has been registered as a biocompute-DM group admin. Your current group password is set to: " + password + " please change it from your administration panel as soon as possible",
                      recipients=[admin_email])
        mail.send(msg)

    def __repr__(self):
        return "<UserGroup %r>" % self.name


class CustomerGroup(Group):
    id = reference_col("Group", primary_key=True)

    parent_id = reference_col("Group", nullable=True)

    parent_group = relationship("UserGroup", backref=backref("customer_groups", uselist=True, lazy="dynamic"), foreign_keys=parent_id, uselist=False)
    samples = relationship("Sample", secondary=customer_to_sample_mapping, lazy="dynamic")
    data_groups = relationship("DataGroup", secondary=customer_to_data_group_mapping, lazy="dynamic")
    data_items = relationship("DataItem", secondary=customer_to_data_item_mapping, lazy="dynamic")
    projects = relationship("Project", secondary=customer_to_project_mapping, lazy="dynamic", backref="customers")

    __tablename__ = "CustomerGroup"
    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Group.id)}

    def __init__(self, group_name, password, admin_name, admin_password, admin_email, parent_group):
        if password is None:
            password = uuid.uuid4().hex

        Group.__init__(self, name=group_name, password=password)
        self.parent_group = parent_group
        user = Customer.create(username=admin_name, email=admin_email, password=admin_password, group=self)
        self.set_administrator(user)

        msg = Message(subject="Welcome to Biocompute-DM",
                      body="Your email address has been registered as a biocompute-DM group admin. Your current group password is set to: " + password + " please change it from your administration panel as soon as possible",
                      recipients=[admin_email])
        mail.send(msg)

    def __repr__(self):
        return "<CustomerGroup %r>" % self.name


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

        msg = Message(subject="Welcome to Biocompute-DM",
                      body="Your email address has been registered as a biocompute-DM user. Your current password is set to: " + password + " please change it from your administration panel as soon as possible",
                      recipients=[email])
        mail.send(msg)

        # Hack, we need the default values here
        db.session.add(self)
        db.session.commit()
        admin_helper.create_user_directory(self, password)

    def __repr__(self):
        return "<Member %r %r>" % (self.username, self.email)

    def get_id(self):
        return self.display_key

    def set_password(self, password):
        if self._password is not None:
            admin_helper.change_password("biocompute-DM_user_" + self.username, password)

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

    submissions = relationship("Submission", backref="user", lazy="dynamic")
    data_groups = relationship("DataGroup", backref="user", lazy="dynamic")
    pipelines = relationship("PipelineInstance", backref="user", lazy="dynamic")
    samples = relationship("Sample", backref="user", lazy="dynamic")
    projects = relationship("Project", backref="user", lazy="dynamic")

    __tablename__ = "User"
    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        if password is None:
            password = uuid.uuid4().hex

        Person.__init__(self, username=username, email=email, password=password, group=group)

    def __repr__(self):
        return "<User %r %r>" % (self.username, self.email)


# Customer table
class Customer(Person):
    id = reference_col("Person", primary_key=True)

    __tablename__ = "Customer"
    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Person.id)}

    def __init__(self, username, email, password, group):
        if password is None:
            password = uuid.uuid4().hex

        Person.__init__(self, username=username, email=email, password=password, group=group)

    def __repr__(self):
        return "<Customer %r %r>" % (self.login_name, self.email)


# Reference data set
class ReferenceData(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    path = Column(String(50), nullable=False)
    description = Column(db.String(500), nullable=False)
    version = Column(String(50), nullable=False)
    current = Column(SmallInteger(), default=False)

    __tablename__ = "ReferenceData"
    __table_args__ = (db.UniqueConstraint("name", "path", "description", "version", name="_unique"),)

    def __init__(self, name, path, description, version):
        db.Model.__init__(self, name=name, path=path, description=description, version=version)

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
