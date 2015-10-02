from sqlalchemy import *
from sqlalchemy.dialects.mysql import *
from migrate import *


from migrate.changeset import schema
pre_meta = MetaData()
post_meta = MetaData()
IDtagLibs = Table('IDtagLibs', post_meta,
    Column('tagID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('libraryName', String(length=40)),
    Column('libraryTagID', String(length=20)),
    Column('tagSequence', String(length=40)),
)

customer = Table('customer', post_meta,
    Column('customerID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('name', String(length=40)),
)

customerContact = Table('customerContact', post_meta,
    Column('contactID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('customerID', INTEGER(display_width=10, unsigned=True)),
    Column('name', String(length=40)),
    Column('email', String(length=50)),
    Column('tel', String(length=30)),
)

demultiplexer = Table('demultiplexer', post_meta,
    Column('demuxID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('status', ENUM('setup', 'running', 'complete', 'finished')),
    Column('location', String(length=500)),
    Column('sourceLocation', String(length=500)),
    Column('JID', INTEGER(display_width=10, unsigned=True)),
)

fastQC = Table('fastQC', post_meta,
    Column('fastQCID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('status', ENUM('setup', 'running', 'complete', 'finished')),
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('sampleID', INTEGER(display_width=10, unsigned=True)),
    Column('location', String(length=500)),
    Column('sourceLocation', String(length=500)),
    Column('JID', INTEGER(display_width=10, unsigned=True)),
)

flowCytProject = Table('flowCytProject', post_meta,
    Column('flowCytProjectID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('flowCytProjectName', String(length=40)),
)

laneData = Table('laneData', post_meta,
    Column('laneID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('laneNumber', TINYINT(display_width=4)),
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('sequencingConc', Float),
    Column('readlClusterDensity', INTEGER(display_width=11)),
    Column('PhiXspiked', Float),
    Column('spike', String(length=20)),
    Column('spikeRation', Float),
)

masterProject = Table('masterProject', post_meta,
    Column('masterProjectID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('status', String(length=20)),
    Column('projectName', String(length=40)),
    Column('projectLead', String(length=40)),
    Column('description', TEXT),
    Column('openDate', Date),
    Column('lastUpdate', Date),
)

masterProjectDocuments = Table('masterProjectDocuments', post_meta,
    Column('documentID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('docDescription', Text),
    Column('docLocation', String(length=500)),
    Column('masterProjectID', INTEGER(display_width=10, unsigned=True)),
)

postAlignQC = Table('postAlignQC', post_meta,
    Column('postAlignQCID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('status', ENUM('setup', 'running', 'complete', 'finished')),
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('sampleID', INTEGER(display_width=10, unsigned=True)),
    Column('location', String(length=500)),
    Column('sourceLocation', String(length=500)),
    Column('JID', INTEGER(display_width=10, unsigned=True)),
)

sampleData = Table('sampleData', post_meta,
    Column('sampleID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('sampleName', String(length=50)),
    Column('tagID', INTEGER(display_width=10, unsigned=True)),
    Column('laneID', INTEGER(display_width=10, unsigned=True)),
    Column('tagSequence', String(length=20)),
    Column('tagKit', String(length=50)),
    Column('analysisID', INTEGER(display_width=10, unsigned=True)),
    Column('adaptorSequence', String(length=200)),
)

seqProject = Table('seqProject', post_meta,
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('seqProjectName', String(length=40)),
    Column('masterProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('seqRunID', INTEGER(display_width=10, unsigned=True)),
    Column('customerID', INTEGER(display_width=10, unsigned=True)),
    Column('exptType', ENUM('exome', 'RNAseq', 'ChIPseq', 'other', 'WGS')),
)

seqRun = Table('seqRun', post_meta,
    Column('seqRunID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('seqRunName', String(length=40)),
    Column('flowcellID', String(length=20)),
    Column('startDate', Date),
    Column('completionDate', Date),
    Column('genomicsLead', String(length=20)),
    Column('dataLocation', String(length=40)),
    Column('indexTagCycles', TINYINT(display_width=4)),
    Column('readCycles', SMALLINT(display_width=6)),
    Column('pairedEnd', ENUM('Y', 'N')),
)

sftpAccount = Table('sftpAccount', post_meta,
    Column('sftpAccountID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('accountLocation', String(length=500)),
    Column('username', String(length=20)),
    Column('userContact', String(length=100)),
    Column('creationDate', Date),
    Column('accountStatus', ENUM('created', 'deleted', 'restore')),
)

state = Table('state', post_meta,
    Column('stateID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('itemID', INTEGER(display_width=10, unsigned=True)),
    Column('state', ENUM('RED', 'GREEN', 'BLUE', 'AMBER', 'BROWN')),
    Column('type', ENUM('masterProject', 'sequencing', 'flow cytometry', 'analysis')),
)

transfer = Table('transfer', post_meta,
    Column('transferID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('seqProjectID', INTEGER(display_width=10, unsigned=True)),
    Column('transLocation', String(length=500)),
    Column('sftpAccountID', INTEGER(display_width=10, unsigned=True)),
    Column('dataStatus', ENUM('origin', 'sftp')),
)

typeLinker = Table('typeLinker', post_meta,
    Column('linkID', INTEGER(display_width=10, unsigned=True), primary_key=True, nullable=False),
    Column('type', ENUM('masterProject', 'sequencing', 'flow cytometry', 'analysis')),
    Column('parentID', INTEGER(display_width=10, unsigned=True)),
    Column('childID', INTEGER(display_width=10, unsigned=True)),
)

user = Table('user', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('username', String(length=100)),
    Column('email', String(length=120)),
    Column('_password', String(length=160)),
    Column('role', String(length=20), default=ColumnDefault('User')),
)


def upgrade(migrate_engine):
    # Upgrade operations go here. Don't create your own engine; bind
    # migrate_engine to your metadata
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['IDtagLibs'].create()
    post_meta.tables['customer'].create()
    post_meta.tables['customerContact'].create()
    post_meta.tables['demultiplexer'].create()
    post_meta.tables['fastQC'].create()
    post_meta.tables['flowCytProject'].create()
    post_meta.tables['laneData'].create()
    post_meta.tables['masterProject'].create()
    post_meta.tables['masterProjectDocuments'].create()
    post_meta.tables['postAlignQC'].create()
    post_meta.tables['sampleData'].create()
    post_meta.tables['seqProject'].create()
    post_meta.tables['seqRun'].create()
    post_meta.tables['sftpAccount'].create()
    post_meta.tables['state'].create()
    post_meta.tables['transfer'].create()
    post_meta.tables['typeLinker'].create()
    post_meta.tables['user'].create()


def downgrade(migrate_engine):
    # Operations to reverse the above upgrade go here.
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['IDtagLibs'].drop()
    post_meta.tables['customer'].drop()
    post_meta.tables['customerContact'].drop()
    post_meta.tables['demultiplexer'].drop()
    post_meta.tables['fastQC'].drop()
    post_meta.tables['flowCytProject'].drop()
    post_meta.tables['laneData'].drop()
    post_meta.tables['masterProject'].drop()
    post_meta.tables['masterProjectDocuments'].drop()
    post_meta.tables['postAlignQC'].drop()
    post_meta.tables['sampleData'].drop()
    post_meta.tables['seqProject'].drop()
    post_meta.tables['seqRun'].drop()
    post_meta.tables['sftpAccount'].drop()
    post_meta.tables['state'].drop()
    post_meta.tables['transfer'].drop()
    post_meta.tables['typeLinker'].drop()
    post_meta.tables['user'].drop()
