from biocomputedm import utils
from biocomputedm.database import *
from biocomputedm.extensions import db

project_sample_association_table = Table("ProjectSampleGrouping",
                                         Column("sample_id", Integer, ForeignKey("Sample.id")),
                                         Column("project_id", Integer, ForeignKey("Project.id")))

project_data_group_association_table = Table("ProjectDataGroupGrouping",
                                         Column("data_group_id", Integer, ForeignKey("DataGroup.id")),
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
    name = Column(String(100), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    user_id = reference_col("User", nullable=True)
    group_id = reference_col("Group", nullable=True)
    pipeline_source_id = reference_col("PipelineInstance", nullable=True)

    pipeline_source = relationship("PipelineInstance", uselist=False, backref=backref("samples", uselist=True))
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
    group_id = reference_col("Group", nullable=True)
    sample_id = reference_col("Sample", nullable=True)

    __tablename__ = "DataItem"

    def __init__(self, name, unlocalised_path, data_group, group):
        db.Model.__init__(self, name=name, unlocalised_path=unlocalised_path)
        data_group.data.append(self)
        data_group.save()
        group.data_items.append(self)
        group.save()


class Project(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text)
    creation_date = Column(Date, default=utils.get_current_date())
    updated_date = Column(Date, onupdate=utils.get_current_date())

    group_id = reference_col("Group")
    creator_id = reference_col("User")

    documents = relationship("Document", backref="project")
    samples = relationship("Sample", secondary=project_sample_association_table, lazy="dynamic", backref=backref("projects", lazy="dynamic"))
    pipeline_outputs = relationship("DataGroup", secondary=project_data_group_association_table, lazy="dynamic")

    __tablename__ = "Project"

    def __init__(self, name, description, creator):
        db.Model.__init__(self, name=name, description=description, user=creator, group=creator.group)

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
