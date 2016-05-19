import os

from biocomputedm import utils
from biocomputedm.database import *
from biocomputedm.extensions import db
from sqlalchemy import update


class Pipeline(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
    author = Column(String(50), nullable=False)
    version = Column(String(50), nullable=False)
    type = Column(Enum("I", "II", "III"), nullable=False)
    regex_type = Column(Enum("AND", "OR"), nullable=False)
    regex = Column(String(100), nullable=False)
    documentation = Column(String(50), nullable=False)
    executable = Column(SmallInteger(), default=False)

    modules = relationship("PipelineModule", backref="pipeline", lazy="dynamic")
    instances = relationship("PipelineInstance", backref="pipeline", lazy="dynamic")

    __tablename__ = "Pipeline"
    __table_args__ = (db.UniqueConstraint("name", "author", "version", name="_unique"),)

    def __init__(self, name, description, author, version, type, regex_type, regex, documentation):
        db.Model.__init__(
            self,
            name=name,
            description=description,
            author=author,
            version=version,
            type=type,
            regex_type=regex_type,
            regex=regex,
            documentation=documentation
        )

    def __repr__(self):
        return "<Pipeline: %s, author: %s, version: %s>" % (self.name, self.author, self.version)

    def generate_unique(self):
        return self.name + " " + self.author + " " + self.version

    def generate_display(self):
        return self.name + "<br/>" + self.author + "<br/>" + self.version


class PipelineModule(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(Text, nullable=False)
    executor = Column(String(500), nullable=False)
    execution_index = Column(db.Integer, nullable=False)

    pipeline_id = reference_col("Pipeline")

    options = relationship("PipelineModuleOption", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModule"
    __table_args__ = (db.UniqueConstraint("pipeline_id", "execution_index", name="_unique"),)

    def __init__(self, name, description, executor, execution_index, pipeline):
        db.Model.__init__(self,
                          name=name,
                          description=description,
                          executor=os.path.join(
                              os.path.join(os.path.join(utils.get_path("scripts", "hpc"), "pipelines"), pipeline.name),
                              executor),
                          execution_index=execution_index)
        pipeline.modules.append(self)

    def __repr__(self):
        return "<Pipeline Module: %s, description %s>" % (self.name, self.description)


class PipelineModuleOption(SurrogatePK, Model):
    display_name = Column(String(50), nullable=False)
    description = Column(Text, nullable=False)
    parameter_name = Column(String(50), nullable=False)
    user_interaction_type = Column(Enum("file", "string", "boolean", "reference", "enum"), nullable=False)
    default_value = Column(String(100), nullable=False)
    necessary = Column(SmallInteger, default=False)

    module_id = reference_col("PipelineModule")

    __tablename__ = "PipelineModuleOption"

    def __init__(self,
                 display_name,
                 description,
                 parameter_name,
                 user_interaction_type,
                 default_value,
                 necessary,
                 module):
        db.Model.__init__(self,
                          display_name=display_name,
                          description=description,
                          parameter_name=parameter_name,
                          user_interaction_type=user_interaction_type,
                          default_value=default_value,
                          necessary=necessary)
        module.options.append(self)

    def __repr__(self):
        return "<Pipeline Option: %s>" % self.display_name


class PipelineInstance(SurrogatePK, Model):
    current_execution_index = Column(Integer, default=-1)
    current_execution_status = Column(Enum("NOT_STARTED", "WAITING", "RUNNING", "STOPPED", "FINISHED", "ERROR"), default="NOT_STARTED")
    execution_type = Column(Enum("Per Module", "Continuous"), default="Continuous")
    options_type = Column(Enum("Custom", "Default"), default="Default")

    pipeline_id = reference_col("Pipeline")
    user_id = reference_col("User")
    group_id = reference_col("Group", nullable=True)

    client_data_group_id = reference_col("DataGroup", nullable=True)
    sample_output_data_group_id = reference_col("DataGroup", nullable=True)
    pipeline_output_data_group_id = reference_col("DataGroup", nullable=True)

    client = relationship("DataGroup", backref="pipeline_instances", foreign_keys=[client_data_group_id], uselist=False)
    sample_output = relationship("DataGroup", backref=backref("pipeline_source", uselist=False), foreign_keys=[sample_output_data_group_id], uselist=False)
    pipeline_output = relationship("DataGroup", foreign_keys=[pipeline_output_data_group_id], uselist=False)

    module_instances = relationship("PipelineModuleInstance", backref="pipeline_instance", lazy="dynamic", cascade="all, delete-orphan")

    __tablename__ = "PipelineInstance"

    def __init__(self, pipeline, execution_type, options_type, user, client):
        db.Model.__init__(self, pipeline=pipeline, execution_type=execution_type, options_type=options_type)
        self.user = user
        self.client = client
        self.save()
        user.group.pipeline_instances.append(self)
        user.group.save()

    def __repr__(self):
        return "<Pipeline Instance for %s at module %s>" % (self.pipeline.name, self.current_execution_index)


class PipelineModuleInstance(SurrogatePK, Model):
    module_id = reference_col("PipelineModule")
    pipeline_instance_id = reference_col("PipelineInstance")
    module_output_data_group_id = reference_col("DataGroup", nullable=True)

    module = relationship("PipelineModule", uselist=False)
    instance = relationship("PipelineInstance", uselist=False)
    option_values = relationship("PipelineModuleOptionValue", backref="module", lazy="dynamic", cascade="all, delete-orphan")
    module_output = relationship("DataGroup", uselist=False)

    __tablename__ = "PipelineModuleInstance"

    def __init__(self, module, instance):
        db.Model.__init__(self, module=module, instance=instance)

    def __repr__(self):
        return "<Pipeline Module Instance for %s, current status: %s>" % (self.module.name, self.current_status)


class PipelineModuleOptionValue(SurrogatePK, Model):
    value = Column(String(100))

    option_id = reference_col("PipelineModuleOption")
    module_instance_id = reference_col("PipelineModuleInstance")

    option = relationship("PipelineModuleOption", uselist=False)
    module_instance = relationship("PipelineModuleInstance", uselist=False)

    __tablename__ = "PipelineModuleOptionValue"

    def __init__(self, option, module_instance):
        db.Model.__init__(self, option=option, module_instance=module_instance)

    def __repr__(self):
        return "<Option Value for: %s with value %s>" % (self.option.display_name, self.value)


def refresh_pipelines():
    # Mark all as legacy on refresh then re-add
    db.session.execute(update(Pipeline, values={Pipeline.executable: False}))
    db.session.commit()

    # Webserver side - as we need the execution paths to be correct for the executors
    path = os.path.join(utils.get_path("scripts", "webserver"), "pipelines")
    directories = os.listdir(path)
    has_new = False
    for directory in directories:
        directory_path = os.path.join(path, directory)
        if not os.path.isdir(directory_path):
            continue

        file = os.path.join(directory_path, directory + ".json")
        if not os.path.isfile(file):
            continue

        from biocomputedm.pipelines.helpers import pipeline_helper as template_helper
        if not template_helper.validate(file):
            continue

        has_new |= template_helper.build(file)

    return has_new


def get_current_module_instance(pipeline_instance):
    if pipeline_instance is None:
        return None

    current_execution_index = pipeline_instance.current_execution_index
    for mod in pipeline_instance.module_instances:
        if mod.module.execution_index == current_execution_index:
            return mod
