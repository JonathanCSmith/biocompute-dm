import os
import subprocess
import time

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.manage import forms
from biocomputedm.manage.models import Submission, get_submissions_query_by_user
from biocomputedm.pipelines.models import Pipeline
from flask import Blueprint, render_template, redirect, url_for
from flask import abort
from flask import current_app
from flask import flash
from flask import request
from flask.ext.login import current_user

manage = Blueprint("manage", __name__, static_folder="static", template_folder="templates")


@manage.route("/manage_message", methods=["GET", "POST"])
def message():
    if request.method == "GET":
        return abort(404)

    # This post contains a message destined for the server - used for hooks etc
    else:
        msg = request.form
        display_key = msg["submission"]
        if display_key is not None and display_key is not "":
            submission = Submission.query.filter_by(display_key=display_key).first()
            if submission is not None:
                submission.validated = True
                submission.update()

        return ""


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


@manage.route("/uploads/<int:page>")
@manage.route("/uploads")
@login_required("ANY")
def uploads(page=1):
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

    return render_template("uploads.html", tile="My Uploads", files=files, page=page, has_prev=p, has_next=n)


# Move data from external into landing_zone
# TODO: Can this be user specific sftp zones?
@manage.route("/download")
@login_required("ANY")
def download():
    return redirect(url_for("empty"))


# Move data from landing_zone into submission (with id)
@manage.route("/new_submission", methods=["GET", "POST"])
@login_required("ANY")
def new_submission():
    # Build the submission form
    form = forms.NewSubmissionForm()

    # Current user information
    folder = current_user.display_key

    # Build path to the users sftp dir
    directory_path = os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], folder)
    directory_path = os.path.join(directory_path, "landing_zone")

    # list of the available files
    filepaths = next(os.walk(directory_path))[2]
    files = []
    for file in filepaths:
        s = os.stat(os.path.join(directory_path, file))
        files.append({
            "name": file,
            "size": s.st_size,
            "date": time.ctime(s.st_ctime)
        })

    # If we are just looking at the page don't perform any of the validation work
    if request.method == "GET":
        return render_template("new_submission.html", title="New Data Submission", form=form, files=files)

    else:
        if form.validate_on_submit():
            # Check to see that at least 1 upload was selected - empty submissions have no use
            ids = request.form.getlist("do_select")
            if ids is None or len(ids) == 0:
                flash("No data uploads were selected.", "warning")
                return render_template("new_submission.html", title="New Data Submission", form=form, files=files)

            else:
                # Create the submission entry
                submission = Submission(name=str(form.submission_name.data),
                                        description=str(form.submission_description.data))

                # Get the user and group
                group = current_user.group
                group.submissions.append(submission)
                current_user.submissions.append(submission)
                submission.save()
                current_user.update()
                group.update()

                # Create the directory to hold the submission
                output_directory_path = os.path.join(current_app.config["SUBMISSIONS_PATH_ON_WEBSERVER"],
                                                     submission.display_key)
                utils.make_directory(output_directory_path)

                # Submit the directory and uploaded file information to our m.u.d. script
                script_path = os.path.join(current_app.config["SCRIPTS_PATH"], "mud.sh")
                sources = ""
                for i in ids:
                    source_path = os.path.join(directory_path, i)
                    sources = sources + source_path + ","

                # Execute our move, unpack and delete script asynchronously so as to not interrupt webserving
                sources = sources[:-1]  # remove that pesky extra comma :D
                process = subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-d=" + output_directory_path,
                            "-s=" + sources,
                            "-i=" + submission.display_key
                        ]
                )  # We are allowing this to execute on it's own - no need to monitor

                # Meanwhile, here we will inform the user and display confirmation
                flash("Submission Successful.", "success")
                return render_template("submission_complete.html", title="Successful Job Submission")

        else:
            utils.flash_errors(form)
            return render_template("new_submission.html", title="New Data Submission", form=form, files=files)


@manage.route("/submissions/<int:page>")
@manage.route("/submissions")
@login_required("ANY")
def submissions(page=1):
    s = get_submissions_query_by_user().paginate(page=page, per_page=20)
    return render_template("submissions.html", title="My Submissions", page=page, obs=s)


@manage.route("/submission")
@manage.route("/submission/<sid>")
@login_required("ANY")
def submission(sid=""):
    if sid == "":
        flash("No submission id was provided", "warning")
        return redirect(url_for("index"))

    submission = Submission.query.filter_by(display_key=sid).first()
    if submission is None:
        flash("Invald submission id", "error")
        return redirect(url_for("index"))

    # Current submission information
    folder = submission.display_key

    # Build path to the users sftp dir
    directory_path = os.path.join(current_app.config["SUBMISSIONS_PATH_ON_WEBSERVER"], folder)

    # list of the available files
    filepaths = next(os.walk(directory_path))[2]
    files = []
    for file in filepaths:
        s = os.stat(os.path.join(directory_path, file))
        files.append({
            "name": file,
            "size": s.st_size,
            "date": time.ctime(s.st_ctime)
        })

    # list of the available type I pipelines
    pipelines = Pipeline.query.filter_by(type="I", executable=True)

    return render_template("submission.html", title="Submission", submission=submission, files=files, pipelines=pipelines)


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
