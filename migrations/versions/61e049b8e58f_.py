"""empty message

Revision ID: 61e049b8e58f
Revises: e57d5b5489cc
Create Date: 2016-05-03 12:15:20.729577

"""

# revision identifiers, used by Alembic.
revision = '61e049b8e58f'
down_revision = 'e57d5b5489cc'

from alembic import op
import sqlalchemy as sa


def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.create_table('ProjectDataGroupGrouping',
    sa.Column('data_group_id', sa.Integer(), nullable=True),
    sa.Column('project_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['data_group_id'], ['DataGroup.id'], ),
    sa.ForeignKeyConstraint(['project_id'], ['Project.id'], )
    )
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('ProjectDataGroupGrouping')
    ### end Alembic commands ###
