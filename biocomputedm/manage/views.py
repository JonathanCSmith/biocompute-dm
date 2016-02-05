import os
import tarfile

import time

import zipfile

from biocomputedm.decorators import login_required
from flask import Blueprint, render_template, redirect, url_for
from flask import current_app
from flask.ext.login import current_user

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
@manage.route("/upload_data/f:<file_uploaded>")
@manage.route("/upload_data/<int:page>")
@manage.route("/upload_data")
@login_required("ANY")
def upload_data(page=1, file_uploaded=""):
    # Current user information
    folder = current_user.display_key

    # Build path to the users sftp dir
    path = os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], folder)
    path = os.path.join(path, "landing_zone")

    # list of the available files
    filepaths = next(os.walk(path))[2]
    if len(filepaths) > 20:
        filepaths = filepaths[(page - 1) * 20:(page * 20) - 1]

        if page == 1:
            p = False
        else:
            p = True

        if (page * 20) - 1 < len(filepaths):
            n = True
        else:
            n = False

    else:
        p = False
        n = False

    files = []
    for file in filepaths:
        s = os.stat(os.path.join(path, file))
        files.append({
            "name": file,
            "size": s.st_size,
            "date": time.ctime(s.st_ctime)
        })

    return render_template("upload.html", tile="Upload your data", files=files, page=page, has_prev=p, has_next=n)


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
