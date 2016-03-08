import os
import subprocess
import time

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.manage import forms
from biocomputedm.manage import models
from biocomputedm.manage.models import Submission, get_submissions_query_by_user, get_samples_query_by_user, \
    get_sample_groups_query_by_user, ReferenceData
from biocomputedm.pipelines.models import Pipeline
from flask import Blueprint, render_template, redirect, url_for
from flask import abort
from flask import current_app
from flask import flash
from flask import request
from flask.ext.login import current_user

manage = Blueprint("manage", __name__, static_folder="static", template_folder="templates")


@manage.route("/manage_message/<oid>", methods=["POST"])
def message(oid=""):
    if request.method == "GET":
        return abort(404)

    # This post contains a message destined for the server - used for hooks etc
    else:
        display_key = oid
        if display_key is not None and display_key is not "":
            submission = Submission.query.filter_by(display_key=display_key).first()
            if submission is not None:
                submission.validated = True
                submission.update()

        return "success"


@manage.route("/refresh_reference_data")
@login_required("Site Admin")
def refresh_reference_data():
    found = models.refresh_reference_data_library()
    if found:
        flash("Successfully loaded all reference data libraries", "success")
    else:
        flash("No reference data libraries were loaded as no new members were identified.", "warning")

    return redirect(url_for("admin.administrate"))


@manage.route("/display_reference_data")
@manage.route("/display_reference_data/<int:page>")
@login_required("ANY")
def display_reference_data(page=1):
    items = ReferenceData.query.paginate(page=page, per_page=20)
    return render_template("reference_libraries.html", title="Reference Data", page=page, obs=items)


@manage.route("/user_profile")
@login_required("ANY")
def user_profile():
    return redirect(url_for("empty"))


@manage.route("/display_data/<data_source_type>|<data_source>|<name>")
@manage.route("/display_data/<data_source_type>|<data_source>|<name>|<data_file>")
@login_required("ANY")
def display_data(data_source_type="", data_source="", name="", data_file=""):
    if data_source_type == "" or data_source == "" or name == "":
        flash("Could not identify the provided data set", "warning")
        return redirect(url_for("empty"))

    if data_source_type == "pipeline_output":
        pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=data_source).first()
        if pipeline_instance is None:
            flash("Could not identify the provided data set", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(os.path.join(os.path.join(utils.get_path("pipeline_data", "serve"), pipeline_instance.display_key), "pipeline_output"), name)

        return render_template("data_viewer.html", data_path=data_path, return_path="pipelines.display_pipeline_instance", oid=pipeline_instance.display_key, data_file="")

    elif data_source_type == "module_output":
        pipeline_instances = current_user.group.pipeline_instances.all()
        p_instance = None
        m_instance = None
        for pipeline_instance in pipeline_instances:
            module_instances = pipeline_instance.module_instances.all()
            for module_instance in module_instances:
                if module_instance.display_key == data_source:
                    p_instance = pipeline_instance
                    m_instance = module_instance
                    break

            if m_instance is not None:
                break

        if m_instance is None:
            flash("Could not locate the provided module instance", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(
            os.path.join(
                os.path.join(
                    os.path.join(
                        utils.get_path("pipeline_data", "serve"),
                        p_instance.display_key
                    ),
                    "modules_output"
                ),
                m_instance.module.name),
            name
        )

        return render_template("data_viewer.html", data_path=data_path, return_path="pipelines.module_instance", oid=m_instance.display_key, data_file="")

    elif data_source_type == "sample_data":
        sample = current_user.group.samples.filter_by(display_key=data_source).first()
        if sample is None:
            flash("Could not identify the provided data set", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(
            os.path.join(
                os.path.join(
                    utils.get_path("sample_data", "serve"),
                    sample.display_key
                ),
                data_file
            ),
            name
        )

        return render_template("data_viewer.html", data_path=data_path, return_path="manage.sample_data", oid=sample.display_key, data_file=data_file)

    else:
        flash("Could not identify the provided data set", "warning")
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


@manage.route("/submissions/<int:page>")
@manage.route("/submissions")
@login_required("ANY")
def submissions(page=1):
    items = get_submissions_query_by_user()
    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("submissions.html", title="Submissions", page=page, obs=items)


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
                output_directory_path = os.path.join(utils.get_path("submission_data", "webserver"),
                                                     submission.display_key)
                utils.make_directory(output_directory_path)

                # Submit the directory and uploaded file information to our m.u.d. script
                script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
                script_path = os.path.join(script_path, "mud.sh")
                sources = ""
                for i in ids:
                    source_path = os.path.join(directory_path, i)
                    sources = sources + source_path + ","
                sources = sources[:-1]  # remove that pesky extra comma :D

                # Execute our move, unpack and delete script asynchronously so as to not interrupt webserving
                subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-d=" + output_directory_path,
                            "-s=" + sources,
                            "-i=" + submission.display_key,
                            "-p=" + current_app.config["LOCAL_WEBSERVER_PORT"]
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                )  # We are allowing this to execute on it's own - no need to monitor

                # Meanwhile, here we will inform the user and display confirmation
                flash("Submission Successful.", "success")
                return render_template("submission_complete.html", title="Successful Job Submission")

        else:
            utils.flash_errors(form)
            return render_template("new_submission.html", title="New Data Submission", form=form, files=files)


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
    directory_path = os.path.join(utils.get_path("submission_data", "webserver"), folder)

    # list of the available files
    filepaths = next(os.walk(directory_path))
    files = []
    for file in filepaths[1]:
        try:
            s = os.stat(os.path.join(directory_path, file))
            files.append({
                "name": file,
                "size": s.st_size,
                "date": time.ctime(s.st_ctime)
            })
        except:
            pass

    for file in filepaths[2]:
        try:
            s = os.stat(os.path.join(directory_path, file))
            files.append({
                "name": file,
                "size": s.st_size,
                "date": time.ctime(s.st_ctime)
            })
        except:
            pass

    # list of the available type I pipelines
    pipelines = Pipeline.query.filter_by(type="I", executable=True)

    return render_template("submission.html", title="Submission", submission=submission, files=files,
                           pipelines=pipelines)


@manage.route("/samples/<int:page>")
@manage.route("/samples")
@login_required("ANY")
def samples(page=1):
    items = get_samples_query_by_user()
    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("samples.html", title="Samples", page=page, obs=items)


@manage.route("/sample/<oid>")
@manage.route("/sample/<oid>|<data_file>")
@login_required("ANY")
def sample(oid="", data_file=""):
    sample = current_user.group.samples.filter_by(display_key=oid).first()
    if sample is None:
        flash("Could not locate the provided sample", "warning")
        return redirect(url_for("empty"))

    data_path = os.path.join(utils.get_path("sample_data", "webserver"), sample.display_key)
    filepaths = next(os.walk(data_path))
    datasets = []
    for file in filepaths[1]:
        try:
            pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=file).first()
            if pipeline_instance is None:
                continue

            item = {
                "name": file,
                "pipeline": pipeline_instance
            }
            datasets.append(item)

        except:
            pass

    if len(datasets) == 0:
        datasets = None

    return render_template("sample.html",
                           title="Sample " + sample.name,
                           sample=sample,
                           data_source_type="sample_pipelines",
                           oid=sample.display_key,
                           datasets=datasets)


@manage.route("/sample_data/<oid>|<data_file>")
@login_required("ANY")
def sample_data(oid="", data_file=""):
    if oid == "" or data_file == "":
        flash("Could not locate the provided sample group", "warning")
        return redirect(url_for("empty"))

    sample = current_user.group.samples.filter_by(display_key=oid).first()
    if sample is None:
        flash("Could not locate the provided sample", "warning")
        return redirect(url_for("empty"))

    data_path = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), sample.display_key), data_file)
    filepaths = next(os.walk(data_path))
    datasets = []
    for file in filepaths[2]:
        try:
            item = {
                "name": file,
                "path": os.path.join(data_path, file)
            }
            datasets.append(item)

        except:
            pass

    if len(datasets) == 0:
        datasets = None

    return render_template("sample_data.html",
                           title="Sample Data " + sample.name,
                           sample=sample,
                           data_source_type="sample_data",
                           oid=sample.display_key,
                           datasets=datasets,
                           data_file=data_file)


@manage.route("/sample_groups/<int:page>")
@manage.route("/sample_groups")
@login_required("ANY")
def sample_groups(page=1):
    items = get_sample_groups_query_by_user()
    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("sample_groups.html", title="Sample Groups", page=page, obs=items)


@manage.route("/new_sample_group")
@login_required("ANY")
def new_sample_group():
    return abort(404)


@manage.route("/sample_group/<oid>")
@login_required("ANY")
def sample_group(oid=""):
    return abort(404)


@manage.route("/projects/<int:page>")
@manage.route("/projects")
@login_required("ANY")
def projects(page=1):
    return abort(404)


@manage.route("/new_project")
@login_required("ANY")
def new_project():
    return abort(404)


@manage.route("/project/<oid>")
@login_required("ANY")
def project(oid=""):
    return abort(404)
