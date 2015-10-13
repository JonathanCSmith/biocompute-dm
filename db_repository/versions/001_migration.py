from sqlalchemy import *
from sqlalchemy.dialects.mysql import *
from migrate import *


from migrate.changeset import schema
pre_meta = MetaData()
post_meta = MetaData()
customer = Table('customer', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('name', String(length=100)),
)

document = Table('document', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('name', String(length=50)),
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

lane = Table('lane', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('project_id', Integer),
    Column('number', Integer),
    Column('sequencing_concentration', Float),
    Column('phi_x_spiked', Float),
    Column('spike', String(length=20)),
    Column('spikeRatio', Float),
)

project = Table('project', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('investigation_id', Integer),
    Column('run_id', Integer),
    Column('project_name', String(length=40)),
    Column('type', String(length=50)),
    Column('customer_id', Integer),
)

run = Table('run', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('name', String(length=40)),
    Column('start_date', Date),
    Column('completion_date', Date),
    Column('data_location', String(length=500)),
    Column('type', String(length=50)),
)

sample = Table('sample', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('internal_sample_name', String(length=50)),
    Column('customer_sample_name', String(length=50)),
    Column('type', String(length=50)),
    Column('project_id', Integer),
    Column('customer_id', Integer),
)

sequencing_project = Table('sequencing_project', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('sequencing_run', Integer),
    Column('sequencing_type', Enum('exome', 'RNAseq', 'ChIPseq', 'WGS', 'other')),
)

sequencing_run = Table('sequencing_run', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('flow_cell_id', String(length=40)),
    Column('genomics_lead', String(length=40)),
    Column('index_tag_cycles', Integer),
    Column('index_tag_cycles_2', Integer),
    Column('read_cycles', Integer),
    Column('read_cycles_2', Integer),
    Column('paired_end', Enum('Yes', 'No')),
)

sequencing_sample = Table('sequencing_sample', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('adaptor_sequence', String(length=200)),
)

tag = Table('tag', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('sample_id', Integer),
    Column('is_first_index', Boolean),
    Column('tag_id', String(length=40)),
    Column('tag_library', String(length=40)),
    Column('tag_sequence', String(length=40)),
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
    post_meta.tables['customer'].create()
    post_meta.tables['document'].create()
    post_meta.tables['flow_cytometry_project'].create()
    post_meta.tables['investigation'].create()
    post_meta.tables['lane'].create()
    post_meta.tables['project'].create()
    post_meta.tables['run'].create()
    post_meta.tables['sample'].create()
    post_meta.tables['sequencing_project'].create()
    post_meta.tables['sequencing_run'].create()
    post_meta.tables['sequencing_sample'].create()
    post_meta.tables['tag'].create()
    post_meta.tables['user'].create()


def downgrade(migrate_engine):
    # Operations to reverse the above upgrade go here.
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['customer'].drop()
    post_meta.tables['document'].drop()
    post_meta.tables['flow_cytometry_project'].drop()
    post_meta.tables['investigation'].drop()
    post_meta.tables['lane'].drop()
    post_meta.tables['project'].drop()
    post_meta.tables['run'].drop()
    post_meta.tables['sample'].drop()
    post_meta.tables['sequencing_project'].drop()
    post_meta.tables['sequencing_run'].drop()
    post_meta.tables['sequencing_sample'].drop()
    post_meta.tables['tag'].drop()
    post_meta.tables['user'].drop()
