import os
import subprocess
import time

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.manage import models
from biocomputedm.manage.models import Submission, get_submissions_query_by_user, get_samples_query_by_user, \
    get_sample_groups_query_by_user, ReferenceData, SampleGroup, Project, Document
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

        data_path = os.path.join(
                os.path.join(os.path.join(utils.get_path("pipeline_data", "serve"), pipeline_instance.display_key),
                             "pipeline_output"), name)

        return render_template("data_viewer.html", data_path=data_path,
                               return_path="pipelines.display_pipeline_instance", oid=pipeline_instance.display_key,
                               data_file="")

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

        return render_template("data_viewer.html", data_path=data_path, return_path="pipelines.module_instance",
                               oid=m_instance.display_key, data_file="")

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

        return render_template("data_viewer.html", data_path=data_path, return_path="manage.sample_data",
                               oid=sample.display_key, data_file=data_file)

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
    from biocomputedm.manage import forms
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


@manage.route("/new_sample_group", methods=["GET", "POST"])
@login_required("ANY")
def new_sample_group():
    from biocomputedm.manage import forms
    form = forms.NewSampleGroupForm()

    if request.method == "GET":
        return render_template("new_sample_group.html", title="New Sample Group", form=form)

    else:
        if form.validate_on_submit():
            pipeline = Pipeline.query.filter_by(display_key=str(form.pipeline.data)).first()
            if pipeline is None:
                flash("Could not locate the selected pipeline", "error")
                return redirect(url_for("index"))

            sample_group = SampleGroup.create(name=str(form.name.data),
                                              creator=current_user,
                                              group=current_user.group,
                                              pipeline=pipeline)

            return redirect(url_for("manage.sample_group", oid=sample_group.display_key))

        else:
            utils.flash_errors(form)
            return render_template("new_sample_group.html", title="New Sample Group", form=form)


@manage.route("/sample_group/<oid>", methods=["GET", "POST"])
@login_required("ANY")
def sample_group(oid=""):
    if oid == "":
        flash("Could not locate the provided sample group", "error")
        return redirect(url_for("index"))

    sample_group = current_user.group.sample_groups.filter_by(display_key=oid).first()
    if sample_group is None:
        flash("Could not locate the provided sample group", "error")
        return redirect(url_for("index"))

    from biocomputedm.manage import forms
    form = forms.UpdateSampleGroupForm()
    if request.method == "POST":
        # Check to see that at least 1 upload was selected - empty submissions have no use
        ids = request.form.getlist("do_select")
        if ids is not None and len(ids) != 0:
            for key in ids:
                new_sample = current_user.group.samples.filter_by(display_key=key).first()
                if new_sample is not None:
                    sample_group.samples.append(new_sample)
                    sample_group.save()

            flash("Sample Group was successfully updated", "success")

        else:
            flash("No samples were selected.", "warning")

    samples = current_user.group.samples.all()
    current_samples = sample_group.samples
    potential_samples = []
    if sample_group.modifiable:
        for sample in samples:
            if sample.pipeline_source.pipeline == sample_group.pipeline and sample not in current_samples:
                potential_samples.append(sample)

    # list of the available type I pipelines
    pipelines = Pipeline.query.filter((Pipeline.type == "II") | (Pipeline.type == "III")).filter_by(executable=True)

    return render_template("sample_group.html", sample_group=sample_group, samples=current_samples,
                           potential_samples=potential_samples, form=form, pipelines=pipelines)


@manage.route("/projects/<int:page>")
@manage.route("/projects")
@login_required("ANY")
def projects(page=1):
    items = current_user.group.projects
    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("projects.html", title="Projects", page=page, obs=items)


@manage.route("/new_project", methods=["GET", "POST"])
@login_required("ANY")
def new_project():
    from biocomputedm.manage import forms
    form = forms.NewProjectForm()
    if request.method == "GET":
        return render_template("new_project.html", title="New Project", form=form)

    else:
        if form.validate_on_submit():
            project = Project.create(name=str(form.investigation_name.data),
                                     description=str(form.investigation_description.data), creator=current_user)
            utils.make_directory(os.path.join(utils.get_path("project_data", "webserver"), project.display_key))
            flash("Investigation successfully registered!", "info")
            return redirect(url_for("investigations"))

        return render_template("new_project.html", title="New Project", form=form)


@manage.route("/project/<oid>", methods=["GET", "POST"])
@login_required("ANY")
def project(oid=""):
    if oid == "":
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    project = current_user.group.projects.filter_by(display_key=oid).first()
    if project is None:
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    from biocomputedm.manage import forms
    form = forms.UpdateProjectForm()
    if request.method == "POST":
        ids = request.form.getlst("do_select")
        if ids is not None and len(ids) != 0:
            for key in ids:
                new_sample_group = current_user.group.sample_groups.filter_by(display_key=key).first()
                if new_sample_group is not None:
                    project.sample_groups.append(new_sample_group)
                    project.save()

            flash("Project was successfully updated", "success")

        else:
            flash("No sample groups were selected", "warning")

    sample_groups = current_user.group.sample_groups.all()
    current_sample_groups = project.sample_groups
    potential_sample_groups = []
    for sample_group in sample_groups:
        if sample_group not in current_sample_groups:
            potential_sample_groups.append(sample_group)

    return render_template("project.html", title="Project", project=project, potential_sample_groups=potential_sample_groups, form=form)


@manage.route("/add_document/<oid>", methods=["GET", "POST"])
@login_required("ANY")
def add_document(oid=""):
    if oid == "":
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    project = current_user.group.projects.filter_by(display_key=oid).first()
    if project is None:
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    from biocomputedm.manage import forms
    form = forms.AddDocumentForm()
    if request.method == "POST":
        if form.validate_on_submit():
            # Handle maliciously named files (i.e. ../..)
            from werkzeug.utils import secure_filename
            filename = secure_filename(form.file_upload.data.filename)
            filepath = os.path.join(os.path.join(utils.get_path("project_data", "webserver"), project.display_key), filename)

            # Handle a document already existing
            if os.path.exists(filepath):
                flash("A document with this location (i.e. filename) already exists", "error")
                return redirect(url_for("investigations"))

            # Save the file to the given path
            form.file_upload.data.save(filepath)

            # Save the document to the db
            document = Document.create(name=str(form.file_upload.data.filename), description=str(form.description.data))
            project.documents.append(document)
            project.save()

            # Inform and redirect
            flash("Document uploaded successfully", "success")
            return redirect(url_for("manage.project", oid=oid))

    # Fail scenario
    return render_template("add_document.html", title="Add Document", form=form, oid=oid)


@manage.route("/remove_document/<oid>|<did>")
@login_required("ANY")
def remove_document(oid="", did=""):
    if oid == "" or did == "":
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    doc = None
    projects = current_user.groups.projects.all()
    for project in projects:
        documents = project.documents
        for document in documents:
            if document.display_key == did:
                doc = document;

    if doc is None:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    path = doc.location
    if os.path.exists(path):
        os.remove(path)

    doc.delete()
    return redirect(url_for("manage.project", oid=oid))
