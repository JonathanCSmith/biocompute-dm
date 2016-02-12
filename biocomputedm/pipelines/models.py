import os

from biocomputedm.database import relationship, reference_col, Model, SurrogatePK, Column, Enum, String, Integer, \
    Boolean
from biocomputedm.extensions import db
from flask import current_app


class Pipeline(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
    author = Column(String(50), nullable=False)
    version = Column(String(50), nullable=False)
    type = Column(Enum("I", "II", "III"), nullable=False)
    executable = Column(Boolean(), default=False)

    modules = relationship("PipelineModule", backref="pipeline", lazy="dynamic")
    instances = relationship("PipelineInstance", backref="pipeline", lazy="dynamic")

    __tablename__ = "Pipeline"
    __table_args__ = (db.UniqueConstraint("name", "description", "author", "version", name="_unique"),)

    def __init__(self, name, description, author, version, type):
        db.Model.__init__(self, name=name, description=description, author=author, version=version, type=type)

    def __repr__(self):
        return "<Pipeline: %s, description: %s, author: %s, version: %s>" % (
            self.name, self.description, self.author, self.version)


class PipelineModule(SurrogatePK, Model):
    name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
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
                          executor=os.path.join(current_app.config["PIPELINES_PATH_ON_HPC"], executor),
                          execution_index=execution_index)
        pipeline.modules.append(self)

    def __repr__(self):
        return "<Pipeline Module: %s, description %s>" % (self.name, self.description)


class PipelineModuleOption(SurrogatePK, Model):
    display_name = Column(String(50), nullable=False)
    description = Column(String(500), nullable=False)
    parameter_name = Column(String(50), nullable=False)
    user_interaction_type = Column(Enum("file", "string", "boolean", "library"), nullable=False)
    default_value = Column(String(50), nullable=False)
    necessary = Column(Boolean, default_value=False)

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
    execution_type = Column(Enum("Per Module", "Continuous"), default="DEFAULT")
    options_type = Column(Enum("All", "Default"), default="Default")

    pipeline_id = reference_col("Pipeline")

    module_instances = relationship("PipelineModuleInstance", backref="pipeline_instance", lazy="dynamic")

    __tablename__ = "PipelineInstance"

    def __init__(self, pipeline):
        db.Model.__init__(self, pipeline=pipeline)

    def __repr__(self):
        return "<Pipeline Instance for %s at module %s>" % (self.pipeline.name, self.current_execution_index)


class PipelineModuleInstance(SurrogatePK, Model):
    current_status = Column(Enum("NOT_STARTED", "RUNNING", "WAITING", "FINISHED", "ERROR"), default="NOT_STARTED")

    module_id = reference_col("PipelineModule")
    pipeline_instance_id = reference_col("PipelineInstance")

    module = relationship("PipelineModule", uselist=False)
    instance = relationship("PipelineInstance", uselist=False)
    option_values = relationship("PipelineModuleOptionValue", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModuleInstance"

    def __init__(self, module, instance):
        db.Model.__init__(self, module=module, instance=instance)

    def __repr__(self):
        return "<Pipeline Module Instance for %s, current status: %s>" % (self.module.name, self.current_status)


class PipelineModuleOptionValue(SurrogatePK, Model):
    value = Column(String(50))

    option_id = reference_col("PipelineModuleOption")
    module_instance_id = reference_col("PipelineModuleInstance")

    option = relationship("PipelineModuleOption", uselist=False)
    module_instance = relationship("PipelineModuleInstance", uselist=False)

    __tablename__ = "PipelineModuleOptionValue"

    def __init__(self, option, module_instance):
        db.Model.__init__(self, option=option, module_instance=module_instance)

    def __repr__(self):
        return "<Option Value for: %s with value %s>" % (self.option.display_name, self.value)
