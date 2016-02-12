from biocomputedm.database import SurrogatePK, Model, reference_col, relationship, Table, Column, Integer, ForeignKey, \
    String
from biocomputedm.extensions import db
from flask.ext.login import current_user


class Submission(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
    validated = Column(db.Boolean, default=False, nullable=False)

    group_id = reference_col("Group")
    submitter_id = reference_col("User")

    __tablename__ = "Submission"

    def __init__(self, name, description):
        db.Model.__init__(self, name=name, description=description)

    def __repr__(self):
        return "<Submission %s>" % self.name


def get_submissions_query_by_user():
    if current_user.is_authenticated:
        return Submission.query.filter_by(submitter=current_user).filter_by(validated=True)


sample_grouping_association_table = Table("sample_grouping",
                                          Column("sample_group_id", Integer, ForeignKey("SampleGroup.id")),
                                          Column("sample_id", Integer, ForeignKey("Sample.id")))


class SampleGroup(SurrogatePK, Model):
    creator_id = reference_col("User")
    group_id = reference_col("Group")

    creator = relationship("User", uselist=False)
    group = relationship("Group", uselist=False)

    __tablename__ = "SampleGroup"


class Sample(SurrogatePK, Model):
    sample_groups = relationship("SampleGroup", secondary=sample_grouping_association_table, backref="samples")

    __tablename__ = "Sample"
