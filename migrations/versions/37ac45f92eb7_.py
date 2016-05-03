"""empty message

Revision ID: 37ac45f92eb7
Revises: 86e41d66b10a
Create Date: 2016-04-28 16:38:58.368524

"""

# revision identifiers, used by Alembic.
revision = '37ac45f92eb7'
down_revision = '86e41d66b10a'

from alembic import op
import sqlalchemy as sa


def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('DataItem', sa.Column('group_id', sa.Integer(), nullable=True))
    op.create_foreign_key(None, 'DataItem', 'Group', ['group_id'], ['id'])
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_constraint(None, 'DataItem', type_='foreignkey')
    op.drop_column('DataItem', 'group_id')
    ### end Alembic commands ###
