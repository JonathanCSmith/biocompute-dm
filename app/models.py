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


# class demultiplex(db.Model):
#     __tablename__ = "demultiplexer"
#     demuxID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
#     # seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
#     status = db.Column(ENUM("setup", "running", "complete", "finished"))
#     location = db.Column(db.String(500))
#     sourceLocation = db.Column(db.String(500))
#     JID = db.Column(INTEGER(10, unsigned=True))
#
#     def __repr__(self):
#         return "<Demultiplexing %r %r %r>" % (self.demuxID, self.seqProjectID, self.status)


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


# User table, currently unused
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(100), unique=True)
    email = db.Column(db.String(120), unique=True)
    _password = db.Column(db.String(160))
    role = db.Column(db.String(20), default="User")

    def __repr__(self):
        return "<User %r %r>" % (self.username, self.email)

    def __init__(self, username, password, email):
        self.username = username.title()
        self.email = email.lower()
        self.set_password(password)
        self.role = False

    def set_password(self, password):
        self._password = generate_password_hash(password)

    def check_password(self, pwd):
        return check_password_hash(self._password, pwd)

    def get_user_role(self):
        return self.role

    def set_user_role(self, role):
        self.role = role


class Investigation(db.Model):
    __tablename__ = "investigation"
    investigation_id = db.Column(db.Integer, primary_key=True)
    investigation_name = db.Column(db.String(40))
    investigation_lead = db.Column(db.String(40))
    investigation_directory = db.Column(db.String(500))
    description = db.Column(db.TEXT)
    open_date = db.Column(db.Date)
    last_update = db.Column(db.Date)
    status = db.Column(db.String(20))

    project = db.RelationshipProperty("Project", backref="investigation")
    document = db.RelationshipProperty("Document", backref="investigation")

    def __init__(self, investigation_name, investigation_lead):
        self.investigation_name = investigation_name
        self.investigation_lead = investigation_lead

        today = datetime.date.today()
        self.open_date = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
        self.last_update = self.open_date

    def set_last_update(self):
        today = datetime.date.today()
        self.last_update = str(today.year) + "-" + str(today.month) + "-" + str(today.day)

    def validate_investigation_directory(self):
        if self.investigation_directory is None:
            self.investigation_directory = os.path.join("./link/investigations",
                                                        str(self.investigation_id) + "_" + self.investigation_name)

        try:
            if not os.path.exists(self.investigation_directory):
                os.mkdir(self.investigation_directory)
                os.chmod(self.investigation_directory, 0o777)
        except OSError as e:
            print(e)
            return None

        docs = os.path.join(self.investigation_directory, "documents")
        try:
            if not os.path.exists(docs):
                os.mkdir(docs)
                os.chmod(docs, 0o777)
        except OSError as e:
            print(e)
            return None

        return docs

    def __repr__(self):
        return "<Investigation %r %r>" % (self.investigation_name, self.investigation_lead)


class Document(db.Model):
    __tablename__ = "document"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50))
    description = db.Column(db.Text)
    location = db.Column(db.String(500))
    investigation_id = db.Column(db.Integer, db.ForeignKey("investigation.investigation_id"))

    def __init__(self, name, description, location):
        self.name = name
        self.description = description
        self.location = location

    def __repr__(self):
        return "<Investigation document %r>" % (self.docDescription)


class Run(db.Model):
    __tablename__ = "run"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(40))
    start_date = db.Column(db.Date)
    completion_date = db.Column(db.Date)
    data_location = db.Column(db.String(500))
    project = db.RelationshipProperty("Project", backref="run", lazy="dynamic")
    type = db.Column(db.String(50))
    __mapper_args__ = {"polymorphic_on": type}


class SequencingRun(Run):
    __tablename__ = "sequencing_run"
    id = db.Column(db.Integer, db.ForeignKey("run.id"), primary_key=True)
    __mapper_args__ = {"polymorphic_identity": "sequencing", "inherit_condition": (id == Run.id)}
    flow_cell_id = db.Column(db.String(40))
    genomics_lead = db.Column(db.String(40))
    index_tag_cycles = db.Column(db.Integer)
    index_tag_cycles_2 = db.Column(db.Integer)
    read_cycles = db.Column(db.Integer)
    read_cycles_2 = db.Column(db.Integer)
    paired_end = db.Column(db.Enum("Yes", "No"))


class Project(db.Model):
    __tablename__ = "project"
    id = db.Column(db.Integer, primary_key=True)
    investigation_id = db.Column(db.Integer, db.ForeignKey("investigation.investigation_id"))
    run_id = db.Column(db.Integer, db.ForeignKey("run.id"))
    name = db.Column(db.String(40))
    type = db.Column(db.String(50))
    customer_id = db.Column(db.Integer, db.ForeignKey("customer.id"))
    sample = db.RelationshipProperty("SequencingSample", backref="project", lazy="dynamic")
    __mapper_args__ = {"polymorphic_on": type}

    def __init__(self):
        return


class SequencingProject(Project):
    __tablename__ = "sequencing_project"
    id = db.Column(db.Integer, db.ForeignKey("project.id"), primary_key=True)
    __mapper_args__ = {"polymorphic_identity": "sequencing", "inherit_condition": (id == Project.id)}
    sequencing_run = db.Column(db.Integer, db.ForeignKey("sequencing_run.id"))
    sequencing_type = db.Column(db.Enum("exome", "RNAseq", "ChIPseq", "WGS", "other"))
    lane = db.RelationshipProperty("Lane", backref="project", lazy="dynamic")

    def __init__(self):
        super(Project, self).__init__()


class Lane(db.Model):
    __tablename__ = "lane"
    id = db.Column(db.Integer, primary_key=True)
    project_id = db.Column(db.Integer, db.ForeignKey("project.id"))
    sample_id = db.Column(db.Integer, db.ForeignKey("sequencing_sample.id"))
    number = db.Column(db.Integer)
    sequencing_concentration = db.Column(db.Float)
    phi_x_spiked = db.Column(db.Float)
    spike = db.Column(db.String(20))
    spike_ratio = db.Column(db.Float)

    def __init__(self, number):
        self.number = number

    def set_sequencing_concentration(self, value):
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

        self.sequencing_concentration = number


class Sample(db.Model):
    __tablename__ = "sample"
    id = db.Column(db.Integer, primary_key=True)
    internal_sample_name = db.Column(db.String(50))
    customer_sample_name = db.Column(db.String(50))
    type = db.Column(db.String(50))
    project_id = db.Column(db.Integer, db.ForeignKey("project.id"))
    customer_id = db.Column(db.Integer, db.ForeignKey("customer.id"))
    __mapper_args__ = {"polymorphic_on": type}


class SequencingSample(Sample):
    __tablename__ = "sequencing_sample"
    id = db.Column(db.Integer, db.ForeignKey("sample.id"), primary_key=True)
    __mapper_args__ = {"polymorphic_identity": "sequencing", "inherit_condition": (id == Sample.id)}

    lane = db.RelationshipProperty("Lane", backref="sample", uselist=False, lazy="dynamic")

    index_tag = db.RelationshipProperty("Tag", backref="sequencing_sample", lazy="dynamic")
    adaptor_sequence = db.Column(db.String(200))


class Tag(db.Model):
    __tablename__ = "tag"
    id = db.Column(db.Integer, primary_key=True)
    sample_id = db.Column(db.Integer, db.ForeignKey("sequencing_sample.id"))
    is_first_index = db.Column(db.Boolean)
    tag_id = db.Column(db.String(40))
    tag_library = db.Column(db.String(40))
    tag_sequence = db.Column(db.String(40))

    def __init__(self, id, kit, seq, is1):
        self.is_first_index = is1
        self.tag_id = id
        self.tag_library = kit
        self.tag_sequence = seq


class Customer(db.Model):
    __tablename__ = "customer"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    project = db.RelationshipProperty("Project", backref="customer", lazy="dynamic")
    sample = db.RelationshipProperty("Sample", backref="customer", lazy="dynamic")

    def __init__(self, name):
        self.name = name


class FlowCytometryProject(Project):
    __tablename__ = "flow_cytometry_project"
    id = db.Column(db.Integer, db.ForeignKey("project.id"), primary_key=True)
    __mapper_args__ = {"polymorphic_identity": "flow_cytometry", "inherit_condition": (id == Project.id)}

    def __init__(self):
        super(Project, self).__init__()
