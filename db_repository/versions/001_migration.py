from sqlalchemy import *
from sqlalchemy.dialects.mysql import *
from migrate import *


from migrate.changeset import schema
pre_meta = MetaData()
post_meta = MetaData()
document = Table('document', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('description', Text),
    Column('location', String(length=500)),
    Column('investigation_id', Integer),
)

flow_cytometry_project = Table('flow_cytometry_project', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
)

investigation = Table('investigation', post_meta,
    Column('investigation_id', Integer, primary_key=True, nullable=False),
    Column('investigation_name', String(length=40)),
    Column('investigation_lead', String(length=40)),
    Column('investigation_directory', String(length=500)),
    Column('description', TEXT),
    Column('open_date', Date),
    Column('last_update', Date),
    Column('status', String(length=20)),
)

project = Table('project', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('investigation_id', Integer),
    Column('project_name', String(length=40)),
    Column('type', String(length=50)),
)

sequencing_project = Table('sequencing_project', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('sequencing_run', Integer),
    Column('customer_id', Integer),
    Column('sequencing_type', Enum('exome', 'RNAseq', 'ChIPseq', 'WGS', 'other')),
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
    post_meta.tables['document'].create()
    post_meta.tables['flow_cytometry_project'].create()
    post_meta.tables['investigation'].create()
    post_meta.tables['project'].create()
    post_meta.tables['sequencing_project'].create()
    post_meta.tables['user'].create()


def downgrade(migrate_engine):
    # Operations to reverse the above upgrade go here.
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['document'].drop()
    post_meta.tables['flow_cytometry_project'].drop()
    post_meta.tables['investigation'].drop()
    post_meta.tables['project'].drop()
    post_meta.tables['sequencing_project'].drop()
    post_meta.tables['user'].drop()
