#!flask/bin/python3
import os.path
import config
import json

from migrate.versioning import api
from app import db, utils

__author__ = 'jon'

db.create_all()
if not os.path.exists(config.SQLALCHEMY_MIGRATE_REPO):
    api.create(config.SQLALCHEMY_MIGRATE_REPO, "database repository")
    api.version_control(config.SQLALCHEMY_DATABASE_URI, config.SQLALCHEMY_MIGRATE_REPO)
else:
    api.version_control(config.SQLALCHEMY_DATABASE_URI, config.SQLALCHEMY_MIGRATE_REPO, api.version(config.SQLALCHEMY_MIGRATE_REPO))

# Build our default user
from app.models import User, Pipeline

admin = User.query.filter_by(login_name=config.SITE_ADMIN_USERNAME).first()

if admin is None:
    admin = User()
    admin.login_name = config.SITE_ADMIN_USERNAME
    admin.set_password(config.SITE_ADMIN_PASSWORD)
    admin.email = config.SITE_ADMIN_EMAIL
    admin.set_role("Site Admin")

    from app.models import Group
    admin_group = Group()
    admin_group.name = "Site Admins"
    admin_group.member.append(admin)

    db.session.add(admin)
    db.session.add(admin_group)
    db.session.commit()

# Build our default pipelines
path = os.path.join(os.path.dirname(__file__), "link")
path = os.path.join(path, "pipelines")
directories = os.listdir(path)
for directory in directories:
    directory_path = os.path.join(path, directory)
    if not os.path.isdir(directory_path):
        continue

    file = os.path.join(directory_path, directory + ".json")
    if not os.path.isfile(file):
        continue

    from app.static.io import pipeline_mappings_template as template_helper
    if not template_helper.validate(file):
        continue

    json_instance = json.loads(open(file).read())
    name = json_instance.get("name")
    description = json_instance.get("description")
    author = json_instance.get("author")
    version = json_instance.get("version")
    type = json_instance.get("pipeline_type")

    pipeline_instance = Pipeline.query.filter_by(name=name, description=description, author=author, version=version, type=type).first()
    if pipeline_instance is not None:
        continue

    pipeline_instance = utils.create_pipeline(name, description, author, version, type)

    for module in json_instance.get("modules"):
        mod = utils.create_module(module.get("name"), module.get("description"), module.get("executor"), module.get("index_in_execution_order"), pipeline_instance)

        for option in module.get("options"):
            opt = utils.create_option(option.get("display_name"), option.get("parameter_name"), option.get("default_value"), option.get("user_interaction_type"), mod)

    db.session.add(pipeline_instance)
    db.session.commit()
