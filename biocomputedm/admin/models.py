from biocomputedm.admin.helpers import helper_functions
from biocomputedm.database import SurrogatePK, Model, reference_col, relationship, String, Column
from biocomputedm.extensions import db
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


# Permissions wrapper & environment contextualiser
class Group(SurrogatePK, Model):
    name = Column(String(50), unique=True, nullable=False)

    members = relationship("Person", backref="group", lazy="dynamic")
    submissions = relationship("Submission", backref="group", lazy="dynamic")

    # # TODO: CHANGE
    # investigation = db.RelationshipProperty("Investigation", backref="group", lazy="dynamic")
    # document = db.RelationshipProperty("Document", backref="group", lazy="dynamic")
    # submission = db.RelationshipProperty("Submission", backref="group", lazy="dynamic")
    # sample_group = db.RelationshipProperty("SampleGroup", backref="group", lazy="dynamic")
    # sample = db.RelationshipProperty("Sample", backref="group", lazy="dynamic")

    __tablename__ = "Group"

    def __init__(self, name):
        db.Model.__init__(self, name=name)

        # Hack, we need the default values here
        db.session.add(self)
        db.session.commit()

        helper_functions.create_group_directory(self)

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
    user = User.create(name=admin_name, email=admin_email, password=admin_password, group=group)
    group.set_administrator(user)
    return group


# Person table, abstract parent for individuals interacting with the software
class Person(UserMixin, SurrogatePK, Model):
    username = Column(String(50), unique=True, nullable=False)
    email = Column(String(120), unique=True, nullable=False)
    role = Column(db.Enum("Member", "Group Admin", "Site Admin"), default="Member")
    type = Column(String(50), nullable=False)

    group_id = reference_col("Group")

    # TODO: CHANGE
    # investigation = db.RelationshipProperty("Investigation", backref="submitter", lazy="dynamic")
    # document = db.RelationshipProperty("Document", backref="submitter", lazy="dynamic")

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
        self.role = role


# User - someone who can submit jobs to the software
class User(Person):
    id = reference_col("Person", primary_key=True)

    submissions = relationship("Submission", backref="submitter", lazy="dynamic")

    # TODO: Change
    # submission = db.RelationshipProperty("Submission", backref="submitter", lazy="dynamic")
    # sample_group = db.RelationshipProperty("SampleGroup", backref="submitter", lazy="dynamic")
    # sample = db.RelationshipProperty("Sample", backref="submitter", lazy="dynamic")

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
