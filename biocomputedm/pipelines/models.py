import os
import uuid

from biocomputedm.extensions import db
from flask import current_app


class Pipeline(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500), nullable=False)
    author = db.Column(db.String(50), nullable=False)
    version = db.Column(db.String(50), nullable=False)
    type = db.Column(db.Enum("I", "II", "III"), nullable=False)

    module = db.RelationshipProperty("PipelineModule", backref="pipeline", lazy="dynamic")
    instance = db.RelationshipProperty("PipelineInstance", backref="pipeline", lazy="dynamic")

    __tablename__ = "Pipeline"
    __table_args__ = (db.UniqueConstraint("name", "description", "author", "version", name="_unique"),)


def create_pipeline(name, description, author, version, type):
    pipeline = Pipeline()
    pipeline.name = name
    pipeline.description = description
    pipeline.author = author
    pipeline.version = version
    pipeline.type = type

    db.session.add(pipeline)
    db.session.commit()
    return pipeline


class PipelineModule(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500), nullable=False)
    executor = db.Column(db.String(500), nullable=False)
    execution_index = db.Column(db.Integer, nullable=False)

    pipeline_id = db.Column(db.Integer, db.ForeignKey("Pipeline.id"), nullable=False)

    module_option = db.RelationshipProperty("PipelineModuleOption", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModule"
    __table_args__ = (db.UniqueConstraint("pipeline_id", "execution_index", name="_unique"),)


def create_module(name, description, executor, order_index, pipeline):
    module = PipelineModule()
    module.name = name
    module.description = description

    # Full path to executor - we have to use the config value here because it may be located at a different place
    # on the remote
    p = os.path.join(current_app.config["PIPELINES_PATH_ON_HPC"], pipeline.name)
    p = os.path.join(p, executor)
    module.executor = p

    module.execution_index = order_index
    pipeline.module.append(module)

    db.session.add(module)
    db.session.add(pipeline)
    db.session.commit()
    return module


class PipelineModuleOption(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    display_name = db.Column(db.String(50), nullable=False)
    paramater_name = db.Column(db.String(50), nullable=False)
    user_interaction_type = db.Column(db.Enum("string", "boolean", "library"), nullable=False)
    default_value = db.Column(db.String(50), nullable=False)

    module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"), nullable=False)

    __tablename__ = "PipelineModuleOption"


def create_option(display_name, parameter_name, default_value, user_interaction_type, module):
    option = PipelineModuleOption()
    option.display_name = display_name
    option.paramater_name = parameter_name
    option.user_interaction_type = user_interaction_type
    option.default_value = default_value
    module.module_option.append(option)

    db.session.add(option)
    db.session.add(module)
    db.session.commit()
    return option


class PipelineInstance(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    current_execution_index = db.Column(db.Integer, default=-1)

    pipeline_id = db.Column(db.Integer, db.ForeignKey("Pipeline.id"))
    current_module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"))

    module_instance = db.RelationshipProperty("PipelineModuleInstance", backref="pipeline_instance", lazy="dynamic")

    __tablename__ = "PipelineInstance"


class PipelineModuleInstance(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    current_status = db.Column(db.Enum("NOT_STARTED", "RUNNING", "WAITING", "FINISHED", "ERRORED"),
                               default="NOT_STARTED")

    module_id = db.Column(db.Integer, db.ForeignKey("PipelineModule.id"), nullable=False)
    pipeline_instance_id = db.Column(db.Integer, db.ForeignKey("PipelineInstance.id"), nullable=False)

    module_option_value = db.RelationshipProperty("PipelineModuleOptionValue", backref="module", lazy="dynamic")

    __tablename__ = "PipelineModuleInstance"


class PipelineModuleOptionValue(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    display_key = db.Column(db.String(32), default=lambda: uuid.uuid4().hex, unique=True)

    value = db.Column(db.String(50))

    pipeline_module_option_id = db.Column(db.Integer, db.ForeignKey("PipelineModuleOption.id"), nullable=False)
    pipeline_module_instance_id = db.Column(db.Integer, db.ForeignKey("PipelineModuleInstance.id"), nullable=False)

    __tablename__ = "PipelineModuleOptionValue"
