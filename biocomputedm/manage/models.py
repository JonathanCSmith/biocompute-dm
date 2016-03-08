import os

from biocomputedm import utils
from biocomputedm.database import *
from biocomputedm.extensions import db
from biocomputedm.pipelines.models import PipelineInstance
from flask.ext.login import current_user
from sqlalchemy import update

pipeline_instance_association_table = Table("PipelineGrouping",
                                            Column("pipeline_instance_id", Integer, ForeignKey("PipelineInstance.id")),
                                            Column("sample_id", Integer, ForeignKey("Sample.id")))

sample_grouping_association_table = Table("SampleGrouping",
                                          Column("sample_group_id", Integer, ForeignKey("SampleGroup.id")),
                                          Column("sample_id", Integer, ForeignKey("Sample.id")))


class DataSource(SurrogatePK, Model):
    current_pipeline = relationship(PipelineInstance, uselist=False, backref="current_data_source")
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
        return current_user.group.submissions.filter_by(validated=True)

    return None


class SampleGroup(DataSource):
    id = reference_col("DataSource", primary_key=True)
    name = Column(String(50), nullable=False)

    creator_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)

    __tablename__ = "SampleGroup"
    __mapper_args__ = {"polymorphic_identity": "SampleGroup", "inherit_condition": (id == DataSource.id)}

    def __init__(self, name, creator, group):
        DataSource.__init__(self, name=name)
        creator.sample_groups.append(self)
        creator.update()
        group.sample_groups.append(self)
        group.update()

    def __repr__(self):
        return "<SampleGroup created by %s for %s>" % (self.creator.name, self.group.name)


def get_sample_groups_query_by_user():
    if current_user.is_authenticated:
        return current_user.group.sample_groups

    return None


class Sample(SurrogatePK, Model):
    name = Column(String(50), nullable=False)

    user_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)
    submission_source_id = reference_col("Submission", nullable=True)
    pipeline_source_id = reference_col("PipelineInstance", nullable=True)

    submission_source = relationship("Submission", uselist=False)
    pipeline_source = relationship("PipelineInstance", uselist=False)
    sample_groups = relationship("SampleGroup", secondary=sample_grouping_association_table, lazy="dynamic",
                                 backref=backref("samples", lazy="dynamic"))
    pipeline_runs = relationship("PipelineInstance", secondary=pipeline_instance_association_table, lazy="dynamic",
                                 backref=backref("samples", lazy="dynamic"))

    __tablename__ = "Sample"

    def __init__(self, name, submission, pipeline):
        db.Model.__init__(self, name=name)
        self.submission_source = submission
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


# Reference data set
class ReferenceData(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(db.String(500), nullable=False)
    version = Column(String(50), nullable=False)
    current = Column(Boolean(), default=False)

    __tablename__ = "ReferenceData"
    __table_args__ = (db.UniqueConstraint("name", "description", "version", name="_unique"),)

    def __init__(self, name, description, version):
        db.Model.__init__(self, name=name, description=description, version=version)

    def __repr__(self):
        return "<Reference Data with name: %s, description; %s and version: %s" % (
            self.name, self.description, self.version)


def refresh_reference_data_library():
    # Mark all as legacy on refresh then re-add
    db.session.execute(update(ReferenceData, values={ReferenceData.current: False}))
    db.session.commit()

    # HPC Side as we need the paths to be correct
    path = utils.get_path("reference_data", "webserver")
    directories = os.listdir(path)
    has_new = False
    for directory in directories:
        directory_path = os.path.join(path, directory)
        if not os.path.isdir(directory_path):
            continue

        file = os.path.join(directory_path, directory + ".json")
        if not os.path.isfile(file):
            continue

        from biocomputedm.manage.helpers import resource_helper as template_helper
        if not template_helper.validate(file):
            continue

        has_new |= template_helper.build(file)

    return has_new