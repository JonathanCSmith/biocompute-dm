__author__ = 'jon'

import os
import datetime

from app import db
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


# class customerContact(db.Model):
#     __tablename__ = "customerContact"
#     contactID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     customerID = db.Column(db.Integer, db.ForeignKey("customer.id"))
#     name = db.Column(db.String(40))
#     email = db.Column(db.String(50))
#     tel = db.Column(db.String(30))
#
#     def __repr__(self):
#         return "<Customer Contact %r %r %r>" % (self.name, self.email, self.tel)


# class postAlignQC(db.Model):
#     __tablename__ = "postAlignQC"
#     postAlignQCID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     status = db.Column(ENUM("setup", "running", "complete", "finished"))
#     # seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
#     # sampleID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sampleData.sampleID"))
#     location = db.Column(db.String(500))
#     sourceLocation = db.Column(db.String(500))
#     JID = db.Column(INTEGER(10, unsigned=True))
#
#     def __repr__(self):
#         return "<Post Align QC %r %r %r>" % (self.seqProjectID, self.sampleID, self.status)

# class sftpAccount(db.Model):
#     __tablename__ = "sftpAccount"
#     sftpAccountID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     accountLocation = db.Column(db.String(500))
#     username = db.Column(db.String(20))
#     userContact = db.Column(db.String(100))
#     creationDate = db.Column(db.Date)
#     accountStatus = db.Column(ENUM("created", "deleted", "restore"))
#
#     def __repr__(self):
#         return "<SFTP Account %r %r %r>" % (self.username, self.creationDate, self.accountStatus)


# class state(db.Model):
#     __tablename__ = "state"
#     stateID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     itemID = db.Column(INTEGER(10, unsigned=True))
#     state = db.Column(ENUM("RED", "GREEN", "BLUE", "AMBER", "BROWN"))
#     type = db.Column(ENUM("masterProject", "sequencing", "flow cytometry", "analysis"))
#
#     def __repr__(self):
#         return "<State %r %r>" % (self.state, self.type)


# class transfer(db.Model):
#     __tablename__ = "transfer"
#     transferID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     # seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
#     transLocation = db.Column(db.String(500))
#     sftpAccountID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sftpAccount.sftpAccountID"))
#     dataStatus = db.Column(ENUM("origin", "sftp"))
#
#     def __repr__(self):
#         return "<Transfer %r %r %r>" % (self.seqProjectID, self.transLocation, self.dataStatus)


# class typeLinker(db.Model):
#     __tablename__ = "typeLinker"
#     linkID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     type = db.Column(ENUM("masterProject", "sequencing", "flow cytometry", "analysis"))
#     parentID = db.Column(INTEGER(10, unsigned=True))
#     childID = db.Column(INTEGER(10, unsigned=True))
#
#     def __repr__(self):
#         return "<Linker %r %r %r>" % (self.type, self.parentID, self.childID)


# class fastQC(db.Model):
#     __tablename__ = "fastQC"
#     fastQCID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     status = db.Column(ENUM("setup", "running", "complete", "finished"))
#     # seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
#     # sampleID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sampleData.sampleID"))
#     location = db.Column(db.String(500))
#     sourceLocation = db.Column(db.String(500))
#     JID = db.Column(INTEGER(10, unsigned=True))
#
#     def __repr__(self):
#         return "<FastQC %r %r %r>" % (self.seqProjectID, self.sampleID, self.status)


# Permissions wrapper & environment contextualiser
class Group(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.String(50), unique=True, nullable=False)

    member = db.RelationshipProperty("Person", backref="group", lazy="dynamic")
    investigation = db.RelationshipProperty("Investigation", backref="group", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="group", lazy="dynamic")
    submission = db.RelationshipProperty("Submission", backref="group", lazy="dynamic")
    project = db.RelationshipProperty("Project", backref="group", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="group", lazy="dynamic")

    def __repr__(self):
        return "<Group %r>" % self.name

    def add_administrator(self, person):
        person.set_user_role("Group Admin")
        self.member.append(person)

    def remove_administrator(self, person):
        person.set_user_role("Member")


# Person table, abstract parent for individuals interacting with the software
class Person(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)

    login_name = db.Column(db.String(50), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
    role = db.Column(db.String(20), default="Member")
    type = db.Column(db.String(50), nullable=False)

    group_id = db.Column(db.Integer, db.ForeignKey("Group.member_id"))

    investigation = db.RelationshipProperty("Investigation", backref="submitter", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="submitter", lazy="dynamic")

    _password = db.Column(db.String(160))

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
    id = db.Column(db.Integer, db.ForeignKey("person.id"), primary_key=True)

    submission = db.RelationshipProperty("Submission", backref="submitter", lazy="dynamic")
    project = db.RelationshipProperty("Project", backref="submitter", lazy="dynamic")

    __mapper_args__ = {"polymorphic_identity": "User", "inherit_condition": (id == Person.id)}

    def __repr__(self):
        return "<User %r %r>" % (self.username, self.email)


# Customer table
class Customer(Person):
    id = db.Column(db.Integer, db.ForeignKey("person.id"), primary_key=True)

    __mapper_args__ = {"polymorphic_identity": "Customer", "inherit_condition": (id == Person.id)}

    def __repr__(self):
        return "<Customer %r %r>" % (self.username, self.email)


class Investigation(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.String(40), nullable=False)
    leader = db.Column(db.String(40), nullable=False)
    directory = db.Column(db.String(500), nullable=False)
    description = db.Column(db.TEXT)
    open_date = db.Column(db.Date, nullable=False)
    last_update = db.Column(db.Date, nullable=False)
    status = db.Column(db.String(20))

    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)

    project = db.RelationshipProperty("Project", backref="investigation", lazy="dynamic")
    document = db.RelationshipProperty("Document", backref="investigation", lazy="dynamic")

    def __init__(self, investigation_name, investigation_lead):
        self.name = investigation_name
        self.lead = investigation_lead

        today = datetime.date.today()
        self.open_date = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
        self.last_update = self.open_date

    def __repr__(self):
        return "<Investigation %r %r>" % (self.investigation_name, self.investigation_lead)

    def set_last_update(self):
        today = datetime.date.today()
        self.last_update = str(today.year) + "-" + str(today.month) + "-" + str(today.day)

    def validate_investigation_directory(self):
        if self.investigation_directory is None:
            path = os.path.join(os.getcwd(), "link")
            path = os.path.join(path, "investigations")
            self.directory = os.path.join(path, str(self.id) + "_" + self.name)

        # Create our investigation directory if necessary
        try:
            if not os.path.exists(self.investigation_directory):
                os.mkdir(self.investigation_directory)
                os.chmod(self.investigation_directory, 0o777)

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

    name = db.Column(db.String(50), nullable=False)
    location = db.Column(db.String(500), nullable=False)
    description = db.Column(db.Text)

    investigation_id = db.Column(db.Integer, db.ForeignKey("Investigation.id"), nullable=False)
    submitter_id = db.Column(db.Integer, db.ForeignKey("Submitter.id"), nullable=False)
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)

    def __repr__(self):
        return "<Investigation document %r>" % (self.docDescription)


class Submission(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.String(40), nullable=False)
    leader = db.Column(db.String(40), nullable=False)
    start_date = db.Column(db.Date, nullable=False)
    completion_date = db.Column(db.Date, nullable=False)
    data_location = db.Column(db.String(500), nullable=False)
    type = db.Column(db.String(50))

    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)

    project = db.RelationshipProperty("Project", backref="run", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="run", lazy="dynamic")

    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Data Submission: %s from %s to %s by %s>" % (self.name, self.start_date, self.completion_date, self.leader)


class SequencingSubmission(Submission):
    id = db.Column(db.Integer, db.ForeignKey("Run.id"), primary_key=True)

    flow_cell_id = db.Column(db.String(40), nullable=False)
    index_tag_cycles = db.Column(db.Integer, nullable=False)
    index_tag_cycles_2 = db.Column(db.Integer, nullable=False)
    read_cycles = db.Column(db.Integer, nullable=False)
    read_cycles_2 = db.Column(db.Integer, nullable=False)
    paired_end = db.Column(db.Enum("Yes", "No"), nullable=False)

    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == Submission.id)}

    def __repr__(self):
        return "<Sequencing Submission for flowcell %s: %s from %s to %s by %s>" % (self.flow_cell_id, self.name, self.start_date, self.completion_date, self.leader)


class Project(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.String(40), nullable=False)
    type = db.Column(db.String(50), nullable=False)

    submission_id = db.Column(db.Integer, db.ForeignKey("Submission.id"), nullable=False)
    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    customer_id = db.Column(db.Integer, db.ForeignKey("Customer.id"), nullable=False)
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)
    investigation_id = db.Column(db.Integer, db.ForeignKey("Investigation.id"))

    sample = db.RelationshipProperty("Sample", backref="project", lazy="dynamic")

    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Project %s of type %s>" % (self.name, self.type)


class SequencingProject(Project):
    id = db.Column(db.Integer, db.ForeignKey("project.id"), primary_key=True)

    # TODO Pipelines!

    # demultiplex_id = db.RelationshipProperty("Demultiplex", backref="project", lazy=False, uselist=False)
    #
    # sequencing_type = db.Column(db.Enum("exome", "RNAseq", "ChIPseq", "WGS", "other"))

    lane = db.RelationshipProperty("Lane", backref="project", lazy="dynamic")

    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == Project.id)}

    def __repr__(self):
        return "<Sequencing Project: %s>" % self.name


class Sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    internal_sample_name = db.Column(db.String(50), nullable=False)
    customer_sample_name = db.Column(db.String(50), nullable=False)
    type = db.Column(db.String(50), nullable=False)

    submitter_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    group_id = db.Column(db.Integer, db.ForeignKey("Group.id"), nullable=False)
    run_id = db.Column(db.Integer, db.ForeignKey("Submission.id"), nullable=False)
    project_id = db.Column(db.Integer, db.ForeignKey("Project.id"), nullable=False)
    customer_id = db.Column(db.Integer, db.ForeignKey("Customer.id"), nullable=False)
    tag_id = db.Column(db.Integer, db.ForeignKey("Tag.id"), nullable=False)

    __mapper_args__ = {"polymorphic_on": type}

    def __repr__(self):
        return "<Sample: %s aka %>s" % (self.internal_sample_name, self.customer_sample_name)


class SequencingSample(Sample):
    id = db.Column(db.Integer, db.ForeignKey("sample.id"), primary_key=True)

    adaptor_sequence = db.Column(db.String(200), nullable=False)

    lane = db.RelationshipProperty("Lane", backref="sample", uselist=False, lazy=False)
    index_tag = db.RelationshipProperty("Tag", backref="sample", lazy="dynamic")

    __mapper_args__ = {"polymorphic_identity": "Sequencing", "inherit_condition": (id == Sample.id)}

    def __repr__(self):
        return "<Sequencing Sample %s aka %s>" % (self.internal_sample_name, self.customer_sample_name)


class Lane(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    number = db.Column(db.Integer, nullable=False)
    sequencing_concentration = db.Column(db.Float, nullable=False)
    phi_x_spiked = db.Column(db.Float, nullable=False)
    spike = db.Column(db.String(20), nullable=False)
    spike_ratio = db.Column(db.Float, nullable=False)

    project_id = db.Column(db.Integer, db.ForeignKey("SequencingProject.id"), nullable=False)
    sample_id = db.Column(db.Integer, db.ForeignKey("SequencingSample.id"), nullable=False)

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

    is_first_index = db.Column(db.Boolean, nullable=False)
    tag_id = db.Column(db.String(40), nullable=False)
    tag_library = db.Column(db.String(40), nullable=False)
    tag_sequence = db.Column(db.String(40), nullable=False)

    sample_id = db.Column(db.Integer, db.ForeignKey("SequencingSample.id"), nullable=False)

    def __repr__(self):
        return "<Tag %s from %s with sequence %s>" % (self.tag_id, self.tag_library, self.tag_sequence)


class FlowCytometryProject(Project):
    id = db.Column(db.Integer, db.ForeignKey("project.id"), primary_key=True)

    __tablename__ = "flow_cytometry_project"
    __mapper_args__ = {"polymorphic_identity": "flow_cytometry", "inherit_condition": (id == Project.id)}

    def __init__(self):
        super(Project, self).__init__()


class Demultiplex(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    project_id = db.Column(db.Integer, db.ForeignKey("sequencing_project.id"))
    argument = db.RelationshipProperty("DemultiplexArgument", backref="demultiplex", lazy="dynamic")

    type = db.Column(db.Enum("BCL2", "Cassava"))
    status = db.Column(db.Enum("Setting Up", "Running", "Run Complete", "Finished"))
    location = db.Column(db.String(500))
    job_id = db.Column(db.Integer)


class DemultiplexArgument(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    parent_id = db.Column(db.Integer, db.ForeignKey(Demultiplex.id))

    key = db.Column(db.String(60))
    value = db.Column(db.String(60))
