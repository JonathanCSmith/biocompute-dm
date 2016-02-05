from biocomputedm.database import SurrogatePK, Model, reference_col
from biocomputedm.extensions import db
from flask.ext.login import current_user


class Submission(SurrogatePK, Model):
    name = db.Column(db.String(50), nullable=False)
    description = db.Column(db.String(500), nullable=False)
    validated = db.Column(db.Boolean, default=False, nullable=False)

    group_id = reference_col("Group")
    submitter_id = reference_col("User")

    __tablename__ = "Submissions"

    def __init__(self, name, description):
        db.Model.__init__(self, name=name, description=description)

    def __repr__(self):
        return "<Submission %s>" % self.name


def get_submissions_query_by_user():
    if current_user.is_authenticated:
        return Submission.query.filter_by(submitter=current_user)
