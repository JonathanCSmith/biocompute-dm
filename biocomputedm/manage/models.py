import datetime

from biocomputedm import utils
from biocomputedm.database import *
from biocomputedm.extensions import db
from flask.ext.login import current_user

project_sample_association_table = Table("ProjectSampleGrouping",
                                         Column("sample_id", Integer, ForeignKey("Sample.id")),
                                         Column("project_id", Integer, ForeignKey("Project.id")))


class Submission(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())
    validated = Column(SmallInteger(), default=False)

    user_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)
    data_group_id = reference_col("DataGroup", nullable=True)

    data_group = relationship("DataGroup", uselist=False)

    __tablename__ = "Submission"

    def __init__(self, **kwargs):
        db.Model.__init__(self, **kwargs)


class Sample(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    user_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)
    pipeline_source_id = reference_col("PipelineInstance", nullable=True)

    pipeline_source = relationship("PipelineInstance", uselist=False, backref="samples")
    data = relationship("DataItem", backref="sample", lazy="dynamic")

    __tablename__ = "Sample"

    def __init__(self, name, pipeline):
        db.Model.__init__(self, name=name)
        self.pipeline_source = pipeline
        self.save()

        pipeline.user.samples.append(self)
        pipeline.user.group.samples.append(self)
        pipeline.save()

    def __repr__(self):
        return "<Sample name: %s>" % self.name


def get_samples_query_by_user():
    if current_user.is_authenticated:
        return current_user.group.samples

    return None


class DataGroup(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    user_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)

    data = relationship("DataItem", backref="data_group", cascade="all, delete-orphan", lazy='dynamic')

    __tablename__ = "DataGroup"

    def __init__(self, name, user, group, source_pipeline):
        db.Model.__init__(self, name=name)
        self.save()
        self.update(pipeline_source=source_pipeline)
        user.data_groups.append(self)
        user.save()
        group.data_groups.append(self)
        group.update()


class DataItem(SurrogatePK, Model):
    name = Column(String(100), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())
    unlocalised_path = Column(String(100), nullable=False)

    data_group_id = reference_col("DataGroup")
    sample_id = reference_col("Sample", nullable=True)

    __tablename__ = "DataItem"

    def __init__(self, name, unlocalised_path, data_group):
        db.Model.__init__(self, name=name, unlocalised_path=unlocalised_path)
        data_group.data.append(self)
        data_group.save()


# class DataSource(SurrogatePK, Model):
#     source_pipeline_id = reference_col("PipelineInstance", nullable=True)
#     running_pipeline_id = reference_col("PipelineInstance", nullable=True)
#
#     source_pipeline = relationship(PipelineInstance, foreign_keys=[source_pipeline_id], uselist=False)
#     running_pipeline = relationship(PipelineInstance, foreign_keys=[running_pipeline_id], uselist=False)
#     run_pipelines = relationship(PipelineInstance, foreign_keys=[PipelineInstance.data_consigner_id], backref=backref("data_consigner", uselist=False))
#
#     type = Column(String(50), nullable=False)
#
#     __tablename__ = "DataSource"
#     __mapper_args__ = {"polymorphic_on": type}
#
#     def __init__(self, **kwargs):
#         db.Model.__init__(self, **kwargs)
#
#
# class Submission(DataSource):
#     id = reference_col("DataSource", primary_key=True)
#
#     name = Column(String(50), nullable=False)
#     description = Column(String(500), nullable=False)
#     validated = Column(SmallInteger, default=False, nullable=False)
#
#     group_id = reference_col("Group")
#     submitter_id = reference_col("User")
#
#     __tablename__ = "Submission"
#     __mapper_args__ = {"polymorphic_identity": "Submission", "inherit_condition": (id == DataSource.id)}
#
#     def __init__(self, name, description):
#         DataSource.__init__(self, name=name, description=description)
#
#     def __repr__(self):
#         return "<Submission %s>" % self.name
#
#
# def get_submissions_query_by_user():
#     if current_user.is_authenticated:
#         return current_user.group.submissions.filter_by(validated=True)
#
#     return None
#
#
# class SampleGroup(DataSource):
#     id = reference_col("DataSource", primary_key=True)
#
#     name = Column(String(50), nullable=False)
#     modifiable = Column(SmallInteger, default=True)
#
#     creator_id = reference_col("User", nullable=True)
#     group_id = reference_col("Group", nullable=True)
#
#     __tablename__ = "SampleGroup"
#     __mapper_args__ = {"polymorphic_identity": "SampleGroup", "inherit_condition": (id == DataSource.id)}
#
#     def __init__(self, name, creator, group, pipeline):
#         DataSource.__init__(self, name=name)
#         self.update(source_pipeline=pipeline)
#         creator.sample_groups.append(self)
#         creator.update()
#         group.sample_groups.append(self)
#         group.update()
#
#     def __repr__(self):
#         return "<SampleGroup created by %s for %s>" % (self.creator.name, self.group.name)
#
#
# def get_sample_groups_query_by_user():
#     if current_user.is_authenticated:
#         return current_user.group.sample_groups
#
#     return None


class Project(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    group_id = reference_col("Group")
    creator_id = reference_col("User")

    documents = relationship("Document", backref="project")
    samples = relationship("Sample", secondary=project_sample_association_table, lazy="dynamic", backref=backref("projects", lazy="dynamic"))

    __tablename__ = "Project"

    def __init__(self, name, description, creator):
        today = datetime.date.today()
        db.Model.__init__(self, name=name, description=description, creator=creator, group=creator.group,
                          open_date=str(today.year) + "-" + str(today.month) + "-" + str(today.day))

    def __repr__(self):
        return "<Project: %s, %s>" % (self.name, self.description)


class Document(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    project_identifier = reference_col("Project", nullable=True)

    __tablename__ = "Document"

    def __init__(self, name, description):
        db.Model.__init__(self, name=name, description=description)

    def __repr__(self):
        return "<Project document %r>" % (self.description)
