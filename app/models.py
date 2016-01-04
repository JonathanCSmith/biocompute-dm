import uuid

__author__ = 'jon'

import os
import datetime
from app import db
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


# Permissions wrapper & environment contextualiser
class Group(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), unique=True, nullable=False)

    member = db.RelationshipProperty("Person", backref="group", lazy="dynamic")
    investigation = db.RelationshipProperty("Investigation", backref="group", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="group", lazy="dynamic")
    submission = db.RelationshipProperty("Submission", backref="group", lazy="dynamic")
    sample_group = db.RelationshipProperty("SampleGroup", backref="group", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="group", lazy="dynamic")

    __tablename__ = "Group"

    def __repr__(self):
        return "<Group %r>" % self.name

    def add_administrator(self, person):
        if person.get_role(self) == "Site Admin":
            return

        person.set_user_role("Group Admin")
        self.member.append(person)

    def remove_administrator(self, person):
        person.set_user_role("Member")


# Person table, abstract parent for individuals interacting with the software
class Person(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    login_name = db.Column(db.String(50), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
    role = db.Column(db.Enum("Member", "Group Admin", "Site Admin"), default="Member")
    type = db.Column(db.String(50), nullable=False)

    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)

    investigation = db.RelationshipProperty("Investigation", backref="submitter", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="submitter", lazy="dynamic")

    _password = db.Column(db.String(160))

    __tablename__ = "Person"
    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Member %r %r>" % (self.username, self.email)

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
    id = db.Column(db.Integer, db.ForeignKey("Person.id"), primary_key=True)

    submission = db.RelationshipProperty("Submission", backref="submitter", lazy="dynamic")
    sample_group = db.RelationshipProperty("SampleGroup", backref="submitter", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="submitter", lazy="dynamic")

    __tablename__ = "User"
    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Person.id)}

    def __repr__(self):
        return "<User %r %r>" % (self.login_name, self.email)


# Customer table
class Customer(Person):
    id = db.Column(db.Integer, db.ForeignKey("Person.id"), primary_key=True)

    __tablename__ = "Customer"
    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Person.id)}

    def __repr__(self):
        return "<Customer %r %r>" % (self.login_name, self.email)


class Investigation(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(40), nullable=False)
    leader = db.Column(db.String(40), nullable=False)
    directory = db.Column(db.String(500))
    description = db.Column(db.TEXT)
    open_date = db.Column(db.Date, nullable=False)
    last_update = db.Column(db.Date, nullable=False)

    submitter_id = db.Column(db.Integer, db.ForeignKey("Person.id"))
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))

    sample_group = db.RelationshipProperty("SampleGroup", backref="investigation", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="investigation", lazy="dynamic")

    __tablename__ = "Investigation"

    def __init__(self, investigation_name, investigation_lead):
        self.name = investigation_name
        self.leader = investigation_lead

        today = datetime.date.today()
        self.open_date = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
        self.last_update = self.open_date

    def __repr__(self):
        return "<Investigation %r %r>" % (self.investigation_name, self.investigation_lead)

    def set_last_update(self):
        today = datetime.date.today()
        self.last_update = str(today.year) + "-" + str(today.month) + "-" + str(today.day)

    def validate_investigation_directory(self):
        if self.directory is None:
            path = os.path.join(os.getcwd(), "link")
            path = os.path.join(path, "investigations")
            self.directory = os.path.join(path, str(self.id) + "_" + self.name)

        # Create our investigation directory if necessary
        try:
            if not os.path.exists(self.directory):
                os.mkdir(self.directory)
                os.chmod(self.directory, 0o777)

        except OSError as e:
            print(e)
            return None

        # Create our documents directory if necessary
        docs = os.path.join(self.directory, "documents")
        try:
            if not os.path.exists(docs):
                os.mkdir(docs)
                os.chmod(docs, 0o777)

        except OSError as e:
            print(e)
            return None

        return docs


class Document(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), nullable=False)
    location = db.Column(db.String(500), nullable=False)
    description = db.Column(db.Text)

    investigation_id = db.Column(db.Integer, db.ForeignKey("Investigation.id"))
    submitter_id = db.Column(db.Integer, db.ForeignKey("Person.id"))
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))

    __tablename__ = "Document"

    def __repr__(self):
        return "<Investigation document %r>" % (self.docDescription)


class Submission(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(40), nullable=False)
    leader = db.Column(db.String(40), nullable=False)
    start_date = db.Column(db.Date, nullable=False)
    completion_date = db.Column(db.Date, nullable=False)
    data_location = db.Column(db.String(500), nullable=False)
    type = db.Column(db.String(50))

    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"))
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))

    sample_group = db.RelationshipProperty("SampleGroup", backref="submission", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="submission", lazy="dynamic")

    __tablename__ = "Submission"
    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Data Submission: %s from %s to %s by %s>" % (
        self.name, self.start_date, self.completion_date, self.leader)


class SequencingSubmission(Submission):
    id = db.Column(db.Integer, db.ForeignKey("Submission.id"), primary_key=True)

    flow_cell_id = db.Column(db.String(40), nullable=False)
    index_tag_cycles = db.Column(db.Integer, nullable=False)
    index_tag_cycles_2 = db.Column(db.Integer, nullable=False)
    read_cycles = db.Column(db.Integer, nullable=False)
    read_cycles_2 = db.Column(db.Integer, nullable=False)
    paired_end = db.Column(db.Enum("Yes", "No"), nullable=False)

    lane = db.RelationshipProperty("Lane", backref="submission", lazy="dynamic")

    __tablename__ = "SequencingSubmission"
    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == Submission.id)}

    def __repr__(self):
        return "<Sequencing Submission for flowcell %s: %s from %s to %s by %s>" % (
        self.flow_cell_id, self.name, self.start_date, self.completion_date, self.leader)


class SampleGroup(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(40), nullable=False)
    type = db.Column(db.String(50), nullable=False)

    submission_id = db.Column(db.Integer, db.ForeignKey("Submission.id"))
    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"))
    customer_id = db.Column(db.Integer, db.ForeignKey("Customer.id"))
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))
    investigation_id = db.Column(db.Integer, db.ForeignKey("Investigation.id"))
    source_pipeline_id = db.Column(db.Integer, db.ForeignKey("PipelineInstance.id"))
    current_pipeline_id = db.Column(db.Integer, db.ForeignKey("PipelineInstance.id"))

    source_pipeline_instance = db.RelationshipProperty("PipelineInstance", uselist=False, foreign_keys=[source_pipeline_id])
    # TODO: Assess if we really want to backref below - I'm not sure of it's use case tbh
    current_pipeline_instance = db.RelationshipProperty(
            "PipelineInstance",
            backref="sample_group",
            uselist=False, foreign_keys=[current_pipeline_id]
    )
    sample = db.RelationshipProperty("Sample", backref="sample_group", lazy="dynamic")

    __tablename__ = "SampleGroup"
    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<SampleGroup %s of type %s>" % (self.name, self.type)


class SequencingSampleGroup(SampleGroup):
    id = db.Column(db.Integer, db.ForeignKey("SampleGroup.id"), primary_key=True)

    # TODO Pipelines!

    # demultiplex_id = db.RelationshipProperty("Demultiplex", backref="project", lazy=False, uselist=False)
    #
    # sequencing_type = db.Column(db.Enum("exome", "RNAseq", "ChIPseq", "WGS", "other"))

    __tablename__ = "SequencingSampleGroup"
    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == SampleGroup.id)}

    def __repr__(self):
        return "<Sequencing SampleGroup: %s>" % self.name


class Sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    internal_sample_name = db.Column(db.String(50), nullable=False)
    customer_sample_name = db.Column(db.String(50), nullable=False)
    sample_type = db.Column(db.String(50))
    type = db.Column(db.String(50), nullable=False)

    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"))
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"))
    submission_id = db.Column(db.Integer, db.ForeignKey("Submission.id"))
    sample_group_id = db.Column(db.Integer, db.ForeignKey("SampleGroup.id"))
    customer_id = db.Column(db.Integer, db.ForeignKey("Customer.id"))
    tag_id = db.Column(db.Integer, db.ForeignKey("Tag.id"))

    __tablename__ = "Sample"
    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Sample: %s aka %s>" % (self.internal_sample_name, self.customer_sample_name)


class SequencingSample(Sample):
    id = db.Column(db.Integer, db.ForeignKey("Sample.id"), primary_key=True)

    adaptor_sequence = db.Column(db.String(200), nullable=False)

    lane_id = db.Column(db.Integer, db.ForeignKey("Lane.id"))

    index_tag = db.RelationshipProperty("Tag", backref="sample", lazy="dynamic")

    __tablename__ = "SequencingSample"
    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == Sample.id)}

    def __repr__(self):
        return "<Sequencing Sample %s aka %s>" % (self.internal_sample_name, self.customer_sample_name)


class Lane(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    number = db.Column(db.Integer, nullable=False)
    sequencing_concentration = db.Column(db.Float, nullable=False)
    phi_x_spiked = db.Column(db.Float, nullable=False)
    spike = db.Column(db.String(20), nullable=False)
    spike_ratio = db.Column(db.Float, nullable=False)

    submission_id = db.Column(db.Integer, db.ForeignKey("SequencingSubmission.id"))

    sample = db.RelationshipProperty("SequencingSample", backref="lane", lazy="dynamic")

    __tablename__ = "Lane"

    def __repr__(self):
        return "<Lane %s with conc %s>" % (self.number, self.sequencing_concentration)

    def set_sequencing_concentration(self, value):
        if isinstance(value, str):
            numeric = '0123456789-.'
            for i, c in enumerate(value):
                if c not in numeric:
                    break

            number = float(value[:i])
            unit = value[i:].lstrip().lower()

            if unit == "um" or unit == "umol":
                number *= 1000000

            elif unit == "nm" or unit == "nmol":
                number *= 1000

        else:
            number = value

        self.sequencing_concentration = number


class Tag(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    is_first_index = db.Column(db.Boolean, nullable=False)
    tag_id = db.Column(db.String(40), nullable=False)
    tag_library = db.Column(db.String(40), nullable=False)
    tag_sequence = db.Column(db.String(40), nullable=False)

    sample_id = db.Column(db.Integer, db.ForeignKey("SequencingSample.id"))

    __tablename__ = "Tag"

    def __repr__(self):
        return "<Tag %s from %s with sequence %s>" % (self.tag_id, self.tag_library, self.tag_sequence)


class Pipeline(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500), nullable=False)
    author = db.Column(db.String(50), nullable=False)
    version = db.Column(db.String(50), nullable=False)
    type = db.Column(db.Enum("I", "II", "III"), nullable=False)

    module = db.RelationshipProperty("PipelineModule", backref="pipeline", lazy="dynamic")
    instance = db.RelationshipProperty("PipelineInstance", backref="pipeline", lazy="dynamic")

    __tablename__ = "Pipeline"
    __table_args__ = (db.UniqueConstraint("name", "description", "author", "version", name="_unique"),)


class PipelineModule(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500), nullable=False)
    executor = db.Column(db.String(500), nullable=False)
    execution_index = db.Column(db.Integer, nullable=False)

    pipeline_id = db.Column(db.Integer, db.ForeignKey("Pipeline.id"), nullable=False)

    module_option = db.RelationshipProperty("PipelineModuleOption", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModule"
    __table_args__ = (db.UniqueConstraint("pipeline_id", "execution_index", name="_unique"),)


class PipelineModuleOption(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    display_name = db.Column(db.String(50), nullable=False)
    paramater_name = db.Column(db.String(50), nullable=False)
    user_interaction_type = db.Column(db.Enum("string", "boolean", "library"), nullable=False)
    default_value = db.Column(db.String(50), nullable=False)

    module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"), nullable=False)

    __tablename__ = "PipelineModuleOption"


class PipelineInstance(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    current_execution_index = db.Column(db.Integer, default=-1)

    pipeline_id = db.Column(db.Integer, db.ForeignKey("Pipeline.id"))
    current_module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"))

    module_instance = db.RelationshipProperty("PipelineModuleInstance", backref="pipeline_instance", lazy="dynamic")

    __tablename__ = "PipelineInstance"


class PipelineModuleInstance(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    current_status = db.Column(db.Enum("NOT_STARTED, RUNNING, WAITING, FINISHED, ERRORED"), default="NOT_STARTED")

    module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"), nullable=False)
    pipeline_instance_id = db.Column(db.Integer, db.ForeignKey("PipelineInstance.id"), nullable=False)

    module_option_value = db.RelationshipProperty("PipelineModuleOptionValue", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModuleInstance"


class PipelineModuleOptionValue(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    value = db.Column(db.String(50))

    pipeline_module_option_id = db.Column(db.Integer, db.ForeignKey("PipelineModuleOption.id"), nullable=False)
    pipeline_module_instance_id = db.Column(db.Integer, db.ForeignKey("PipelineModuleInstance.id"), nullable=False)

    __tablename__ = "PipelineModuleOptionValue"
