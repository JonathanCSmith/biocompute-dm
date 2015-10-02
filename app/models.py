__author__ = 'jon'

import datetime

from app import db
from sqlalchemy.dialects.mysql import INTEGER, ENUM, TINYINT, SMALLINT
from flask.ext.login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


# User table, currently unused
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
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


class IDtagLibs(db.Model):
    __tablename__ = "IDtagLibs"
    tagID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True, nullable=False)
    libraryName = db.Column(db.String(40))
    libraryTagID = db.Column(db.String(20))
    tagSequence = db.Column(db.String(40))

    def __repr__(self):
        return "<ID Tag %r %r %r>" % (self.libraryName, self.libraryTagID, self.tagSequence)


class customer(db.Model):
    __tablename__ = "customer"
    customerID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    name = db.Column(db.String(40))

    def __repr__(self):
        return "<Customer %r>" % (self.name)


class customerContact(db.Model):
    __tablename__ = "customerContact"
    contactID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    customerID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("customer.customerID"))
    name = db.Column(db.String(40))
    email = db.Column(db.String(50))
    tel = db.Column(db.String(30))

    def __repr__(self):
        return "<Customer Contact %r %r %r>" % (self.name, self.email, self.tel)


class flowCytProject(db.Model):
    __tablename__ = "flowCytProject"
    flowCytProjectID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    flowCytProjectName = db.Column(db.String(40))

    def __repr__(self):
        return "<Flow Cyt Project %r>" % (self.flowCytProjectName)


class laneData(db.Model):
    __tablename__ = "laneData"
    laneID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    laneNumber = db.Column(TINYINT(4))
    seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
    sequencingConc = db.Column(db.Float)
    readlClusterDensity = db.Column(INTEGER(11))
    PhiXspiked = db.Column(db.Float)
    spike = db.Column(db.String(20))
    spikeRation = db.Column(db.Float)

    def __repr__(self):
        return "<Lane data %r %r>" % (self.seqProjectID, self.laneID)


class MasterProject(db.Model):
    __tablename__ = "masterProject"
    masterProjectID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    status = db.Column(db.String(20))
    projectName = db.Column(db.String(40))
    projectLead = db.Column(db.String(40))
    description = db.Column(db.TEXT)
    openDate = db.Column(db.Date)
    lastUpdate = db.Column(db.Date)

    def __init__(self, projectName, projectLead):
        self.projectName = projectName
        self.projectLead = projectLead

        today = datetime.date.today()
        self.openDate = str(today.year) + "-" + str(today.month) + "-" + str(today.day)
        self.lastUpdate = self.openDate

    def set_last_update(self):
        today = datetime.date.today()
        self.lastUpdate = str(today.year) + "-" + str(today.month) + "-" + str(today.day)

    def __repr__(self):
        return "<Master project %r %r>" % (self.projectName, self.projectLead)


class masterProjectDocuments(db.Model):
    __tablename__ = "masterProjectDocuments"
    documentID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    docDescription = db.Column(db.Text)
    docLocation = db.Column(db.String(500))
    masterProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("masterProject.masterProjectID"))

    def __repr__(self):
        return "<Master project document %r>" % (self.docDescription)


class postAlignQC(db.Model):
    __tablename__ = "postAlignQC"
    postAlignQCID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    status = db.Column(ENUM("setup", "running", "complete", "finished"))
    seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
    sampleID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sampleData.sampleID"))
    location = db.Column(db.String(500))
    sourceLocation = db.Column(db.String(500))
    JID = db.Column(INTEGER(10, unsigned=True))

    def __repr__(self):
        return "<Post Align QC %r %r %r>" % (self.seqProjectID, self.sampleID, self.status)


class sampleData(db.Model):
    __tablename__ = "sampleData"
    sampleID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    sampleName = db.Column(db.String(50))
    tagID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("IDtagLibs.tagID"))
    laneID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("laneData.laneID"))
    tagSequence = db.Column(db.String(20))
    tagKit = db.Column(db.String(50))
    analysisID = db.Column(INTEGER(10, unsigned=True))
    adaptorSequence = db.Column(db.String(200))

    def __repr__(self):
        return "<Sample Data %r>" % (self.sampleName)


class SequencingProject(db.Model):
    __tablename__ = "seqProject"
    seqProjectID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    seqProjectName = db.Column(db.String(40))
    masterProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("masterProject.masterProjectID"))
    seqRunID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqRun.seqRunID"))
    customerID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("customer.customerID"))
    exptType = db.Column(ENUM("exome", "RNAseq", "ChIPseq", "other", "WGS"))

    def __repr__(self):
        return "<Sequencing Project %r %r>" % (self.seqProjectName, self.exptType)


class seqRun(db.Model):
    __tablename__ = "seqRun"
    seqRunID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    seqRunName = db.Column(db.String(40))
    flowcellID = db.Column(db.String(20))
    startDate = db.Column(db.Date)
    completionDate = db.Column(db.Date)
    genomicsLead = db.Column(db.String(20))
    dataLocation = db.Column(db.String(40))
    indexTagCycles = db.Column(TINYINT(4))
    readCycles = db.Column(SMALLINT(6))
    pairedEnd = db.Column(ENUM("Y", "N"))

    def __repr__(self):
        return "<Sequence Run %r %r %r>" % (self.seqRunName, self.startDate, self.completionDate)


class sftpAccount(db.Model):
    __tablename__ = "sftpAccount"
    sftpAccountID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    accountLocation = db.Column(db.String(500))
    username = db.Column(db.String(20))
    userContact = db.Column(db.String(100))
    creationDate = db.Column(db.Date)
    accountStatus = db.Column(ENUM("created", "deleted", "restore"))

    def __repr__(self):
        return "<SFTP Account %r %r %r>" % (self.username, self.creationDate, self.accountStatus)


class state(db.Model):
    __tablename__ = "state"
    stateID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    itemID = db.Column(INTEGER(10, unsigned=True))
    state = db.Column(ENUM("RED", "GREEN", "BLUE", "AMBER", "BROWN"))
    type = db.Column(ENUM("masterProject", "sequencing", "flow cytometry", "analysis"))

    def __repr__(self):
        return "<State %r %r>" % (self.state, self.type)


class transfer(db.Model):
    __tablename__ = "transfer"
    transferID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
    transLocation = db.Column(db.String(500))
    sftpAccountID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sftpAccount.sftpAccountID"))
    dataStatus = db.Column(ENUM("origin", "sftp"))

    def __repr__(self):
        return "<Transfer %r %r %r>" % (self.seqProjectID, self.transLocation, self.dataStatus)


class typeLinker(db.Model):
    __tablename__ = "typeLinker"
    linkID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    type = db.Column(ENUM("masterProject", "sequencing", "flow cytometry", "analysis"))
    parentID = db.Column(INTEGER(10, unsigned=True))
    childID = db.Column(INTEGER(10, unsigned=True))

    def __repr__(self):
        return "<Linker %r %r %r>" % (self.type, self.parentID, self.childID)


class demultiplex(db.Model):
    __tablename__ = "demultiplexer"
    demuxID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
    status = db.Column(ENUM("setup", "running", "complete", "finished"))
    location = db.Column(db.String(500))
    sourceLocation = db.Column(db.String(500))
    JID = db.Column(INTEGER(10, unsigned=True))

    def __repr__(self):
        return "<Demultiplexing %r %r %r>" % (self.demuxID, self.seqProjectID, self.status)


class fastQC(db.Model):
    __tablename__ = "fastQC"
    fastQCID = db.Column(INTEGER(10, unsigned=True), primary_key=True, autoincrement=True)
    status = db.Column(ENUM("setup", "running", "complete", "finished"))
    seqProjectID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("seqProject.seqProjectID"))
    sampleID = db.Column(INTEGER(10, unsigned=True), db.ForeignKey("sampleData.sampleID"))
    location = db.Column(db.String(500))
    sourceLocation = db.Column(db.String(500))
    JID = db.Column(INTEGER(10, unsigned=True))

    def __repr__(self):
        return "<FastQC %r %r %r>" % (self.seqProjectID, self.sampleID, self.status)
