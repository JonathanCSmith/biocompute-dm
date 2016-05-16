"""empty message

Revision ID: b9f1fa575c0c
Revises: 9a2a8df1bda0
Create Date: 2016-05-13 13:53:51.691627

"""

# revision identifiers, used by Alembic.
revision = 'b9f1fa575c0c'
down_revision = '9a2a8df1bda0'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('DataGroup', 'running')
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('DataGroup', sa.Column('running', mysql.SMALLINT(display_width=6), autoincrement=False, nullable=True))
    ### end Alembic commands ###