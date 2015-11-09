#!flask/bin/python3
__author__ = 'jon'

import os.path
import config

from migrate.versioning import api
from app import db

db.create_all()
if not os.path.exists(config.SQLALCHEMY_MIGRATE_REPO):
    api.create(config.SQLALCHEMY_MIGRATE_REPO, "database repository")
    api.version_control(config.SQLALCHEMY_DATABASE_URI, config.SQLALCHEMY_MIGRATE_REPO)
else:
    api.version_control(config.SQLALCHEMY_DATABASE_URI, config.SQLALCHEMY_MIGRATE_REPO, api.version(config.SQLALCHEMY_MIGRATE_REPO))

# Build our default user
from app.models import User
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