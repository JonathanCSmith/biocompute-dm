"""empty message

Revision ID: ea6764b46df4
Revises: 61e049b8e58f
Create Date: 2016-05-03 16:36:16.322752

"""

# revision identifiers, used by Alembic.
revision = 'ea6764b46df4'
down_revision = '61e049b8e58f'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.create_table('CustomerGroup',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('parent_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['Group.id'], ),
    sa.ForeignKeyConstraint(['parent_id'], ['Group.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('UserGroup',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['Group.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('CustomerDataGroupMapping',
    sa.Column('customer_id', sa.Integer(), nullable=True),
    sa.Column('data_group_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['CustomerGroup.id'], ),
    sa.ForeignKeyConstraint(['data_group_id'], ['DataGroup.id'], )
    )
    op.create_table('CustomerProjectMapping',
    sa.Column('customer_id', sa.Integer(), nullable=True),
    sa.Column('project_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['CustomerGroup.id'], ),
    sa.ForeignKeyConstraint(['project_id'], ['Project.id'], )
    )
    op.create_table('CustomerSampleMapping',
    sa.Column('customer_id', sa.Integer(), nullable=True),
    sa.Column('sample_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['CustomerGroup.id'], ),
    sa.ForeignKeyConstraint(['sample_id'], ['Sample.id'], )
    )
    op.create_table('CustomerDataItemMapping',
    sa.Column('customer_id', sa.Integer(), nullable=True),
    sa.Column('data_item_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['CustomerGroup.id'], ),
    sa.ForeignKeyConstraint(['data_item_id'], ['DataItem.id'], )
    )
    op.drop_table('SampleToCustomer')
    op.drop_table('DataGroupToCustomer')
    op.drop_table('ProjectToCustomer')
    op.add_column('Group', sa.Column('type', sa.String(length=50), nullable=False))
    op.drop_constraint('Group_ibfk_1', 'Group', type_='foreignkey')
    op.drop_column('Group', 'parent_id')
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('Group', sa.Column('parent_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True))
    op.create_foreign_key('Group_ibfk_1', 'Group', 'Group', ['parent_id'], ['id'])
    op.drop_column('Group', 'type')
    op.create_table('ProjectToCustomer',
    sa.Column('project_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.Column('customer_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['Customer.id'], name='ProjectToCustomer_ibfk_1'),
    sa.ForeignKeyConstraint(['project_id'], ['Project.id'], name='ProjectToCustomer_ibfk_2'),
    mysql_default_charset='latin1',
    mysql_engine='InnoDB'
    )
    op.create_table('DataGroupToCustomer',
    sa.Column('data_group_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.Column('customer_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['Customer.id'], name='DataGroupToCustomer_ibfk_1'),
    sa.ForeignKeyConstraint(['data_group_id'], ['DataGroup.id'], name='DataGroupToCustomer_ibfk_2'),
    mysql_default_charset='latin1',
    mysql_engine='InnoDB'
    )
    op.create_table('SampleToCustomer',
    sa.Column('sample_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.Column('customer_id', mysql.INTEGER(display_width=11), autoincrement=False, nullable=True),
    sa.ForeignKeyConstraint(['customer_id'], ['Customer.id'], name='SampleToCustomer_ibfk_1'),
    sa.ForeignKeyConstraint(['sample_id'], ['Sample.id'], name='SampleToCustomer_ibfk_2'),
    mysql_default_charset='latin1',
    mysql_engine='InnoDB'
    )
    op.drop_table('CustomerDataItemMapping')
    op.drop_table('CustomerSampleMapping')
    op.drop_table('CustomerProjectMapping')
    op.drop_table('CustomerDataGroupMapping')
    op.drop_table('UserGroup')
    op.drop_table('CustomerGroup')
    ### end Alembic commands ###
