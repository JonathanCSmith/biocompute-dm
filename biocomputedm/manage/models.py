from biocomputedm.database import *
from biocomputedm.extensions import db
from biocomputedm.pipelines.models import PipelineInstance
from flask.ext.login import current_user


class DataSource(SurrogatePK, Model):
    current_pipeline = relationship(PipelineInstance, uselist=False)
    past_pipelines = relationship(PipelineInstance)
    type = Column(String(50), nullable=False)

    __tablename__ = "DataSource"
    __mapper_args__ = {"polymorphic_on": type}

    def __init__(self, **kwargs):
        db.Model.__init__(self, **kwargs)


class Submission(DataSource):
    id = reference_col("DataSource", primary_key=True)

    name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
    validated = Column(db.Boolean, default=False, nullable=False)

    group_id = reference_col("Group")
    submitter_id = reference_col("User")

    __tablename__ = "Submission"
    __mapper_args__ = {"polymorphic_identity": "Submission", "inherit_condition": (id == DataSource.id)}

    def __init__(self, name, description):
        DataSource.__init__(self, name=name, description=description)

    def __repr__(self):
        return "<Submission %s>" % self.name


def get_submissions_query_by_user():
    if current_user.is_authenticated:
        return Submission.query.filter_by(submitter=current_user).filter_by(validated=True)


sample_grouping_association_table = Table("sample_grouping",
                                          Column("sample_group_id", Integer, ForeignKey("SampleGroup.id")),
                                          Column("sample_id", Integer, ForeignKey("Sample.id")))


class SampleGroup(DataSource):
    id = reference_col("DataSource", primary_key=True)

    creator_id = reference_col("User")
    group_id = reference_col("Group")

    creator = relationship("User", uselist=False)
    group = relationship("Group", uselist=False)

    __tablename__ = "SampleGroup"
    __mapper_args__ = {"polymorphic_identity": "SampleGroup", "inherit_condition": (id == DataSource.id)}

    def __init__(self, creator, group):
        DataSource.__init__(self, creator=creator, group=group)

    def __repr__(self):
        return "<SampleGroup created by %s for %s>" % (self.creator.name, self.group.name)


class Sample(SurrogatePK, Model):
    sample_groups = relationship("SampleGroup", secondary=sample_grouping_association_table, backref="samples")

    __tablename__ = "Sample"
