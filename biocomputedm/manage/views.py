from biocomputedm.decorators import login_required
from flask import Blueprint, render_template, redirect, url_for

manage = Blueprint("manage", __name__, static_folder="static", template_folder="templates")


# @manage.route("/submissions", methods=["GET", "POST"])
# @manage.route("/submissions/<int:page>", methods=["GET", "POST"])
# @login_required("ANY")
# def submissions(page=1):
#     s = utils.get_allowed_submissions_query().paginate(page=page, per_page=20)
#     return render_template("submissions.html", title="My Runs", page=page, obs=s)


@manage.route("/user_profile")
@login_required("ANY")
def user_profile():
    return redirect(url_for("empty"))


# Move data from external into landing_zone
# TODO: Can this be user specific sftp zones?
@manage.route("/upload_data")
@login_required("ANY")
def upload_data():
    return render_template("upload.html", tile="Upload your data")


# Move data from external into landing_zone
# TODO: Can this be user specific sftp zones?
@manage.route("/download")
@login_required("ANY")
def download():
    return redirect(url_for("empty"))


# Move data from landing_zone into submission (with id)
# TODO: Should we delete the original
@manage.route("/new_submission", methods=["GET", "POST"])
@login_required("ANY")
def new_submission():
    return render_template("new_submission.html", title="New Data Submission")


@manage.route("/submissions")
@manage.route("/submissions/<int:page>")
@login_required("ANY")
def submissions(page=1):
    return redirect(url_for("empty"))


@manage.route("/process")
@manage.route("/process/<sg_id>")
@login_required("ANY")
def process(sg_id=""):
    return redirect(url_for("empty"))


@manage.route("/processed")
@manage.route("/processed/<int:page>")
@login_required("ANY")
def processed(page=1):
    return redirect(url_for("empty"))


@manage.route("/analyse")
@manage.route("/analyse/<sg_id>")
@login_required("ANY")
def analyse(sg_id=""):
    return redirect(url_for("empty"))


@manage.route("/analysed")
@manage.route("/analysed/<int:page>")
@login_required("ANY")
def analysed(page=1):
    return redirect(url_for("empty"))



