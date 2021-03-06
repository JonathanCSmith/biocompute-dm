import os
import re
import subprocess
import time

from flask import Blueprint, render_template, redirect, url_for, abort, current_app, flash, request, send_from_directory
from flask.ext.login import current_user
from flask.ext.mail import Message

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.extensions import mail
from biocomputedm.manage.helpers.manage_helper import asynchronously_calculate_viable_pipelines_for_data_group, asynchronously_calculate_viable_pipelines_for_submission
from biocomputedm.manage.models import Submission, Project, Document, DataGroup, DataItem, Sample
from biocomputedm.pipelines.models import Pipeline

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

                # Build a data group for the submission
                data_group = DataGroup.create(
                    name="Data Group submitted by: " + submission.user.username,
                    user=submission.user,
                    group=submission.group,
                    source_pipeline=None
                )

                # Walk the submission directory to find data sets
                local_data_path = os.path.join(utils.get_path("submission_data", "webserver"), submission.display_key)
                paths = next(os.walk(local_data_path, followlinks=False))
                for folder in paths[1]:
                    data_item = DataItem.create(
                        name=folder,
                        unlocalised_path=os.path.join(submission.display_key),
                        data_group=data_group,
                        group=submission.group
                    )

                for file in paths[2]:
                    data_item = DataItem.create(
                        name=file,
                        unlocalised_path=os.path.join(submission.display_key),
                        data_group=data_group,
                        group=submission.group
                    )

                # Async submission to pipeline calculator
                asynchronously_calculate_viable_pipelines_for_submission(current_app._get_current_object(), submission.display_key)

                # Ensure the submission knows about the data group
                submission.update(validated=True, data_group=data_group)

                # Good place to let the user know their submission is done - technically the calc valid pipelines stuff isn't but here is where we have user context
                msg = Message(subject="Data Successfully Submitted",
                              body="Your data for " + submission.name + " has been successfully submitted.",
                              recipients=[submission.user.email])

                mail.send(msg)

        return "success"


@manage.route("/display_data/<item_id>|<data_type>")
@login_required("ANY")
def display_data(item_id="", data_type=""):
    if item_id == "" or data_type == "":
        flash("Could not identify the provided data set", "warning")
        return redirect(url_for("empty"))

    data_item = None
    file_path = None
    data_path = None
    if data_type == "sample":
        if current_user.get_role() == "Site Admin":
            data_item = DataItem.query.filter_by(display_key=item_id).first()

        else:
            data_item = current_user.group.data_items.filter_by(display_key=item_id).first()

        if data_item is None:
            flash("Could not identify the provided data set", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(os.path.join(utils.get_path("sample_data", "serve"), data_item.unlocalised_path), data_item.name)
        file_path = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), data_item.unlocalised_path), data_item.name)

    elif data_type == "pipeline":
        if current_user.get_role() == "Site Admin":
            data_item = DataItem.query.filter_by(display_key=item_id).first()

        else:
            data_item = current_user.group.data_items.filter_by(display_key=item_id).first()

        if data_item is None:
            flash("Could not identify the provided data set", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(os.path.join(utils.get_path("pipeline_data", "serve"), data_item.unlocalised_path), data_item.name)
        file_path = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), data_item.unlocalised_path), data_item.name)

    elif data_type == "module":
        if current_user.get_role() == "Site Admin":
            data_item = DataItem.query.filter_by(display_key=item_id).first()

        else:
            data_item = current_user.group.data_items.filter_by(display_key=item_id).first()

        if data_item is None:
            flash("Could not identify the provided data set", "warning")
            return redirect(url_for("empty"))

        data_path = os.path.join(os.path.join(utils.get_path("pipeline_data", "serve"), data_item.unlocalised_path), data_item.name)
        file_path = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), data_item.unlocalised_path), data_item.name)

    elif data_type == "module_alt":
        data_parts = item_id.split("___")
        if len(data_parts) is not 3:
            flash("Bad URL for module data provided", "warning")
            return redirect(url_for("empty"))

        sub_path = os.path.join(os.path.join(os.path.join(data_parts[0], "modules_output"), data_parts[1]), data_parts[2])
        data_path = os.path.join(utils.get_path("pipeline_data", "serve"), sub_path)
        file_path = os.path.join(utils.get_path("pipeline_data", "webserver"), sub_path)

    if data_path is None or file_path is None:
        flash("Could not identify the provided data set's path", "warning")
        return redirect(url_for("empty"))

    if not os.path.isfile(file_path):
        flash("Viewing nested data is currently not supported, please download the folder instead", "warning")
        return redirect(url_for("empty"))

    if data_type != "module_alt":
        stats = os.stat(file_path)
        if stats.st_size > 5242880: #5mb
            flash("The file is too big to download from the browser, please use the SFTP process instead.", "info")
            return redirect(url_for("empty"))

    return render_template("data_viewer.html", title="Data Display", data_item=data_item, data_path=data_path)


@manage.route("/staged_files")
@login_required("ANY")
def staged_files():
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    path = request.url_root
    if path.startswith("http://"):
        path = path.replace("http://", "")

    if path.startswith("https://"):
        path = path.replace("https://", "")

    if path.endswith("/"):
        path = path[:-1]

    if ":" in path:
        splitz = path.split(":")
        # if path.startswith("https://"):
        #     path = splitz[0] + ":" + splitz[1]
        # else:
        path = splitz[0]

    # TODO: Can make this a config option that allows us to exclude certain addresses for external port config
    # Allows us to have a different port for LAN and WAN - may be necessary with certain configs
    if path.startswith("10") or path.startswith("127"):
        port = "22"
    else:
        port = current_app.config["EXTERNAL_SFTP_PORT"]

    return render_template("staged_files.html", title="My Files", path=path, port=port)


@manage.route("/submissions/<int:page>")
@manage.route("/submissions")
@login_required("ANY")
def submissions(page=1):
    if current_user.get_role() == "Site Admin":
        items = Submission.query.filter_by(validated=True)
    else:
        items = current_user.group.submissions.filter_by(validated=True)

    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("submissions.html", title="Submissions", page=page, obs=items)


@manage.route("/new_submission", methods=["GET", "POST"])
@login_required("ANY")
def new_submission():
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    # Build the submission form
    from biocomputedm.manage import forms
    form = forms.NewSubmissionForm()

    # Current user information
    folder = current_user.display_key

    # Build path to the users sftp dir
    directory_path = os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], folder)
    directory_path = os.path.join(directory_path, "staged_files")

    # list of the available files in user directory
    filepaths = next(os.walk(directory_path))
    files = []

    for file in filepaths[1]:
        path = os.path.join(directory_path, file)
        s = os.stat(path)
        files.append({
            "name": file,
            "path": path,
            "size": os.path.getsize(path),
            "date": time.ctime(s.st_ctime)
        })

    for file in filepaths[2]:
        # rx = re.compile("^(.*)(?<!\.fastq)\.(7z|bz2|deb|gz|tar|tbz2|xz|tgz|rar|zip|Z)$")
        # if rx.match(file):
        path = os.path.join(directory_path, file)
        s = os.stat(path)
        files.append({
            "name": file,
            "path": path,
            "size": s.st_size,
            "date": time.ctime(s.st_ctime)
        })

    # List the available files in group directory
    directory_path = os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], current_user.group.name)
    directory_path = os.path.join(directory_path, "staged_files")
    filepaths = next(os.walk(directory_path))

    for file in filepaths[1]:
        path = os.path.join(directory_path, file)
        s = os.stat(path)
        files.append({
            "name": file,
            "path": path,
            "size": os.path.getsize(path),
            "date": time.ctime(s.st_ctime)
        })

    for file in filepaths[2]:
        # rx = re.compile("^(.*)(?<!\.fastq)\.(7z|bz2|deb|gz|tar|tbz2|xz|tgz|rar|zip|Z)$")
        # if rx.match(file):
        path = os.path.join(directory_path, file)
        s = os.stat(path)
        files.append({
            "name": file,
            "path": path,
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
                submission = Submission(name=str(form.submission_name.data), description=str(form.submission_description.data))

                # Unpack option
                unpack = bool(form.submission_unpack.data)
                if unpack:
                    unpack = "True"
                else:
                    unpack = "False"

                # Get the user and group
                group = current_user.group
                group.submissions.append(submission)
                current_user.submissions.append(submission)
                submission.save()
                current_user.update()
                group.update()

                # Create the directory to hold the submission
                output_directory_path = os.path.join(utils.get_path("submission_data", "webserver"), submission.display_key)
                utils.make_directory(output_directory_path)

                # Submit the directory and uploaded file information to our m.u.d. script
                script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
                script_path = os.path.join(script_path, "mud.sh")
                sources = ""
                for file in files:
                    for i in ids:
                        if file["name"] == i:
                            sources = sources + file["path"] + ","
                            break
                sources = sources[:-1]  # remove that pesky extra comma :D

                # Execute our move, unpack and delete script asynchronously so as to not interrupt webserving
                with open(os.devnull, "w") as fnull:
                    subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-d=" + output_directory_path,
                            "-s=" + sources,
                            "-i=" + submission.display_key,
                            "-p=" + current_app.config["LOCAL_WEBSERVER_PORT"],
                            "-u=" + unpack
                        ],
                        stdout=fnull,
                        stderr=fnull
                    )  # We are allowing this to execute on it's own - no need to monitor

                # In the meantime we will inform the user and display confirmation
                flash("Submission Successful.", "success")
                return render_template("submission_complete.html", title="Successful Job Submission")

        else:
            utils.flash_errors(form)
            return render_template("new_submission.html", title="New Data Submission", form=form, files=files)


@manage.route("/submission/<oid>|<int:option>")
@manage.route("/submission/<oid>")
@login_required("ANY")
def submission(oid="", option=0):
    if oid == "":
        flash("No submission id was provided", "warning")
        return redirect(url_for("index"))

    if current_user.get_role() == "Site Admin":
        submission = Submission.query.filter_by(display_key=oid).first()
    else:
        submission = current_user.group.submissions.filter_by(display_key=oid).first()

    if submission is None:
        flash("Invald submission id", "error")
        return redirect(url_for("index"))

    if option == 1:
        asynchronously_calculate_viable_pipelines_for_submission(current_app._get_current_object(), submission.display_key)
        flash("Your submission is currently recalculating it's available pipelines. This may take some time. Please wait 5 minutes before re-trying!", "success")
        return redirect(url_for("index"))

    # # list of the available type I pipelines
    # pipelines = Pipeline.query.filter_by(type="I", executable=True)
    #
    # # use the regex information to determine if a pipeline can run
    # valid_pipelines = []
    # for pipeline in pipelines:
    #     rxs = pipeline.regex.split("###")
    #
    #     if pipeline.regex_type == "OR":
    #         valid_pipeline = False
    #         for rx in rxs:
    #             found = False
    #             rx_matcher = re.compile(rx)
    #             for dirpath, dirnames, files in os.walk(os.path.join(utils.get_path("submission_data", "webserver"), oid)):
    #                 for file in files:
    #                     if rx_matcher.match(file):
    #                         found = True
    #                         break
    #
    #                 if found:
    #                     break
    #
    #             if found:
    #                 valid_pipeline = True
    #                 break
    #
    #         if valid_pipeline:
    #             valid_pipelines.append(pipeline)
    #             continue
    #
    #     elif pipeline.regex_type == "AND":
    #         valid_pipeline = True
    #         for rx in rxs:
    #             found = False
    #             rx_matcher = re.compile(rx)
    #             for dirpath, dirnames, files in os.walk(os.path.join(utils.get_path("submission_data", "webserver"), oid)):
    #                 for file in files:
    #                     if rx_matcher.match(file):
    #                         found = True
    #                         break
    #
    #                 if found:
    #                     break
    #
    #             if not found:
    #                 valid_pipeline = False
    #                 break
    #
    #         if valid_pipeline:
    #             valid_pipelines.append(pipeline)
    #             continue

    valid_pipelines = []
    pipeline_string = submission.data_group.valid_pipelines
    if pipeline_string is not None:
        pipeline_keys = pipeline_string.split(",")
        for pipeline_key in pipeline_keys:
            # We check for executable here as well as the list may be old
            pipeline = Pipeline.query.filter_by(display_key=pipeline_key, executable=True).first()
            if pipeline is not None:
                valid_pipelines.append(pipeline)

    running_pipelines = []
    if submission.data_group:
        for run_pipeline in submission.data_group.pipeline_instances:
            if run_pipeline.current_execution_status != "STOPPED" and run_pipeline.current_execution_status != "FINISHED":
                running_pipelines.append(run_pipeline)
                break

    return render_template("submission.html", title="Submission", submission=submission, pipelines=valid_pipelines, running_pipelines=running_pipelines)


@manage.route("/remove_submission/<oid>|<int:force>")
@login_required("ANY")
def remove_submission(oid="", force=0):
    if force != 1:
        return render_template(
            "confirm.html",
            message="Are you sure you wish to remove this submission?",
            oid=oid,
            url="manage.remove_submission"
        )

    if oid == "":
        flash("Could not identify the provided object", "warning")
        return redirect(url_for("manage.submissions"))

    if current_user.get_role() == "Site Admin":
        submission = Submission.query.filter_by(display_key=oid).first()
    else:
        submission = current_user.group.submissions.filter_by(display_key=oid).first()

    if submission is None:
        flash("Could not identify the provided object", "warning")
        return redirect(url_for("manage.submissions"))

    # Execute our delete script synchronously
    script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
    script_path = os.path.join(script_path, "delete.sh")
    source = os.path.join(utils.get_path("submission_data", "webserver"), submission.display_key)
    subprocess.Popen(
        [
            "sudo",
            script_path,
            "-s=" + source
        ]
    ).wait()

    submission.delete()

    flash("Submission deletion was successful", "success")
    return redirect(url_for("manage.submissions"))


@manage.route("/samples/<int:page>")
@manage.route("/samples")
@login_required("ANY")
def samples(page=1):
    if current_user.get_role() == "Site Admin":
        items = Sample.query
    else:
        items = current_user.group.samples

    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("samples.html", title="Samples", page=page, obs=items)


@manage.route("/sample/<oid>")
@login_required("ANY")
def sample(oid=""):
    if oid == "":
        flash("Could not locate the provided sample", "warning")
        return redirect(url_for("empty"))

    if current_user.get_role() == "Site Admin":
        sample = Sample.query.filter_by(display_key=oid).first()
    else:
        sample = current_user.group.samples.filter_by(display_key=oid).first()

    if sample is None:
        flash("Could not locate the provided sample", "warning")
        return redirect(url_for("empty"))

    return render_template("sample.html", title="Sample " + sample.name, sample=sample)


@manage.route("/delete_sample/<oid>")
@manage.route("/delete_sample/<oid>|<int:force>")
@login_required("Site Admin")
def delete_sample(oid="", force=0):
    if force != 1:
        return render_template(
            "confirm.html",
            message="Are you sure you wish to remove this sample?",
            oid=oid,
            url="manage.delete_sample"
        )

    if oid == "":
        flash("No instance identifiers were provided.", "warning")
        return redirect(url_for("empty"))

    if current_user.get_role() == "Site Admin":
        sample = Sample.query.filter_by(display_key=oid).first()
    else:
        flash("You do not have permission to do this.", "warning")
        return redirect(url_for("index"))

    if sample is None:
        flash("Could not identify the sample", "warning")
        return redirect(url_for("index"))

    for data in sample.data:
        data_group = data.data_group
        data.delete()
        if data_group.data is None or not data_group.data:
            data_group.delete()

    sample.delete()

    subprocess.Popen(
        [
            "sudo",
            os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "delete.sh"),
            "-s=" + os.path.join(utils.get_path("sample_data", "webserver"), oid)
        ]
    )

    flash("Sample was deleted successfully", "success")
    return redirect(url_for("manage.samples"))


@manage.route("/data_groups/<int:page>")
@manage.route("/data_groups")
@login_required("ANY")
def data_groups(page=1):
    if current_user.get_role() == "Site Admin":
        items = DataGroup.query.filter(DataGroup.pipeline_source != None)
    else:
        items = current_user.group.data_groups.filter(DataGroup.pipeline_source != None)

    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("data_groups.html", title="Data Groups", page=page, obs=items)


# @manage.route("/new_sample_group", methods=["GET", "POST"])
# @login_required("ANY")
# def new_sample_group():
#     from biocomputedm.manage import forms
#     form = forms.NewSampleGroupForm()
#
#     samples = current_user.group.samples.all()
#
#     if request.method == "GET":
#         return render_template("new_sample_group.html", title="New Sample Group", form=form, samples=samples)
#
#     else:
#         if form.validate_on_submit():
#             # Check that at least one sample was selected
#             ids = request.form.getlist("do_select")
#             if ids is None or len(ids) == 0:
#                 flash("No samples were selected to add to the group", "warning")
#                 return render_template("new_sample_group.html", title="New Sample Group", form=form, samples=samples)
#
#             else:
#                 sample_group = SampleGroup(
#                         name=str(form.name.data),
#                         creator=current_user,
#                         group=current_user.group,
#                         pipeline=None
#                 )
#
#                 for i in ids:
#                     sample = Sample.query.filter_by(display_key=i).first()
#                     sample_group.samples.append(sample)
#
#                 flash("Sample Group creation was successful", "success")
#                 return redirect(url_for("manage.sample_group", oid=sample_group.display_key))
#
#         else:
#             utils.flash_errors(form)
#             return render_template("new_sample_group.html", title="New Sample Group", form=form, samples=samples)


# TODO: Pagination of samples?
@manage.route("/data_group/<oid>|<data_type>|<int:option>")
@manage.route("/data_group/<oid>|<data_type>")
@login_required("ANY")
def data_group(oid="", data_type="", option=0):
    if oid == "" or data_type == "":
        flash("Could not locate the provided data group", "error")
        return redirect(url_for("index"))

    if current_user.get_role() == "Site Admin":
        data_group = DataGroup.query.filter_by(display_key=oid).first()
        if data_group is None:
            flash("Could not locate the provided data group", "error")
            return redirect(url_for("index"))

        if option == 1:
            asynchronously_calculate_viable_pipelines_for_data_group(current_app._get_current_object(), data_group.display_key)
            flash("Your data group is currently recalculating it's available pipelines. This may take some time. Please wait 5 minutes before re-trying!")
            return redirect(url_for("index"))

        valid_pipelines = None

    else:
        data_group = current_user.group.data_groups.filter_by(display_key=oid).first()
        if data_group is None:
            flash("Could not locate the provided data group", "error")
            return redirect(url_for("index"))

        if option == 1:
            asynchronously_calculate_viable_pipelines_for_data_group(current_app._get_current_object(), data_group.display_key)
            flash("Your data group is currently recalculating it's available pipelines. This may take some time. Please wait 5 minutes before re-trying!", "success")
            return redirect(url_for("index"))

        # pipelines = Pipeline.query.filter((Pipeline.type == "II") | (Pipeline.type == "III")).filter_by(executable=True)
        #
        # # use the regex information to determine if a pipeline can run
        # valid_pipelines = []
        # root_path = utils.get_path("sample_data", "webserver")
        # files = []
        # for data in data_group.data:
        #     data_path = os.path.join(os.path.join(root_path, data.unlocalised_path), data.name)
        #     if not os.path.isfile(data_path):
        #         for dirpath, dirnames, additional_files in os.walk(data_path):
        #             files.extend(additional_files)
        #
        #     else:
        #         files.append(data.name)
        #
        # for pipeline in pipelines:
        #     rxs = pipeline.regex.split("###")
        #
        #     if pipeline.regex_type == "OR":
        #         valid_pipeline = False
        #         for rx in rxs:
        #             found = False
        #             rx_matcher = re.compile(rx)
        #
        #             for file in files:
        #                 if rx_matcher.match(file):
        #                     found = True
        #                     break
        #
        #             if found:
        #                 valid_pipeline = True
        #                 break
        #
        #         if valid_pipeline:
        #             valid_pipelines.append(pipeline)
        #             continue
        #
        #     elif pipeline.regex_type == "AND":
        #         valid_pipeline = True
        #         for rx in rxs:
        #             found = False
        #             rx_matcher = re.compile(rx)
        #             for file in files:
        #                 if rx_matcher.match(file):
        #                     found = True
        #                     break
        #
        #             if not found:
        #                 valid_pipeline = False
        #                 break
        #
        #         if valid_pipeline:
        #             valid_pipelines.append(pipeline)
        #             continue

        valid_pipelines = []
        key_string = data_group.valid_pipelines
        if key_string is not None:
            pipeline_keys = key_string.split(",")
            for pipeline_key in pipeline_keys:
                # We check for executable here as well as the list may be old
                pipeline = Pipeline.query.filter_by(display_key=pipeline_key, executable=True).first()
                if pipeline is not None:
                    valid_pipelines.append(pipeline)

    running_pipelines = []
    for run_pipeline in data_group.pipeline_instances:
        if run_pipeline.current_execution_status != "STOPPED" and run_pipeline.current_execution_status != "FINISHED":
            running_pipelines.append(run_pipeline)
            break

    return render_template("data_group.html", data_group=data_group, pipelines=valid_pipelines, running_pipelines=running_pipelines, data_type=data_type)


@manage.route("/projects/<int:page>")
@manage.route("/projects")
@login_required("ANY")
def projects(page=1):
    if current_user.get_role() == "Site Admin":
        items = Project.query
    else:
        items = current_user.group.projects

    if items is not None:
        items = items.paginate(page=page, per_page=20)

    return render_template("projects.html", title="Projects", page=page, obs=items)


@manage.route("/new_project", methods=["GET", "POST"])
@login_required("ANY")
def new_project():
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    from biocomputedm.manage import forms
    form = forms.NewProjectForm()
    if request.method == "GET":
        return render_template("new_project.html", title="New Project", form=form)

    else:
        if form.validate_on_submit():
            project = Project.create(name=str(form.investigation_name.data), description=str(form.investigation_description.data), creator=current_user)
            utils.make_directory(os.path.join(utils.get_path("project_data", "webserver"), project.display_key))
            flash("Project successfully registered.", "info")
            return redirect(url_for("manage.project", oid=project.display_key))

        return render_template("new_project.html", title="New Project", form=form)


@manage.route("/project/<oid>|<did>")
@manage.route("/project/<oid>")
@login_required("ANY")
def project(oid="", did=""):
    if oid == "":
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    if current_user.get_role() == "Site Admin":
        project = Project.query.filter_by(display_key=oid).first()

    else:
        project = current_user.group.projects.filter_by(display_key=oid).first()

    if project is None:
        flash("Could not identify the provided project.", "error")
        return redirect(url_for("index"))

    if did != "":
        for document in project.documents:
            if document.display_key == did:
                return send_from_directory(os.path.join(utils.get_path("project_data", "webserver"), project.display_key), document.name, as_attachment=True, attachment_filename=document.name)

    return render_template("project.html", title="Project", project=project)


@manage.route("/link_to_project/<oid>|<data_type>|<int:page>", methods=["GET", "POST"])
@manage.route("/link_to_project/<oid>|<data_type>", methods=["GET", "POST"])
@login_required("ANY")
def link_to_project(page=1, oid="", data_type=""):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    if oid == "":
        flash("Could not link the provided information to a project", "error")
        return redirect(url_for("index"))

    if data_type != "sample" and data_type != "sample_data" and data_type != "pipeline":
        flash("Could not link the provided information to a project", "error")
        return redirect(url_for("index"))

    if request.method == "GET":
        from biocomputedm.manage import forms
        form = forms.SelectProjectForm()

        projects = current_user.group.projects
        if projects is not None:
            projects = projects.paginate(page=page, per_page=20)

        if projects is None:
            flash("There are no projects available to you, please create one first.", "error")
            return redirect(url_for("index"))

        return render_template("link_to_project.html", page=page, projects=projects, oid=oid, data_type=data_type, form=form)

    else:
        # Check to see that at least 1 upload was selected - empty submissions have no use
        ids = request.form.getlist("do_select")
        if ids is None or len(ids) == 0:
            flash("No project was selected.", "warning")

            projects = current_user.group.projects
            if projects is not None:
                projects = projects.paginate(page=page, per_page=20)

            if projects is None or len(projects) == 0:
                flash("There are no projects available to you, please create one first.", "error")
                return redirect(url_for("index"))

            return render_template("link_to_project.html", page=page, projects=projects, oid=oid, data_type=data_type)

        else:
            # Per project appending of the data
            has_real_project = False
            for i in ids:
                project = current_user.group.projects.filter_by(display_key=i).first()
                if project is None:
                    continue

                has_real_project = True

                # Switch behaviour based on the data item
                if data_type == "sample":
                    sample = current_user.group.samples.filter_by(display_key=oid).first()
                    if sample is None:
                        flash("Could not identify the sample to link to the project", "warning")
                        return redirect(url_for("index"))

                    project.samples.append(sample)

                elif data_type == "sample_data":
                    data_group = current_user.group.data_groups.filter_by(display_key=oid).first()
                    if data_group is None:
                        flash("Could not identify the sample set to link to the project", "warning")
                        return redirect(url_for("index"))

                    for customer in project.customers:
                        customer.data_groups.append(data_group)
                        customer.save()

                        for data in data_group.data:
                            customer.data_items.append(data)

                    has_real_data = False
                    for data in data_group.data:
                        if data.sample:
                            project.samples.append(data.sample)
                            has_real_data = True

                        for customer in project.customers:
                            customer.samples.append(data.sample)
                            customer.data_items.append(data)
                            customer.save()

                    if not has_real_data:
                        flash("None of the data items in the provided group had samples to link", "warning")
                        return redirect(url_for("index"))

                else:
                    data_group = current_user.group.data_groups.filter_by(display_key=oid).first()
                    if data_group is None:
                        flash("Could not identify the data group to link to the project.", "warning")
                        return redirect(url_for("index"))

                    project.pipeline_outputs.append(data_group)
                    for customer in project.customers:
                        customer.data_groups.append(data_group)
                        customer.save()

                project.save()

            if not has_real_project:
                flash("No valid projects were provided", "warning")
                return redirect(url_for("index"))

            flash("Projects were updated successfully", "success")
            return redirect(url_for("index"))


@manage.route("/add_document/<oid>", methods=["GET", "POST"])
@login_required("ANY")
def add_document(oid=""):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

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
            filename = secure_filename(form.file_upload.data.filename).replace(" ", "_")
            filepath = os.path.join(os.path.join(utils.get_path("project_data", "webserver"), project.display_key),
                                    filename)

            # Handle a document already existing
            if os.path.exists(filepath):
                flash("A document with this location (i.e. filename) already exists", "error")
                return redirect(url_for("projects"))

            # Save the file to the given path
            form.file_upload.data.save(filepath)

            # Save the document to the db
            document = Document.create(name=filename, description=str(form.description.data))
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
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    if oid == "" or did == "":
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    doc = None
    pro = None
    projects = current_user.group.projects.all()
    for project in projects:
        documents = project.documents
        for document in documents:
            if document.display_key == did:
                doc = document
                pro = project

    if doc is None:
        flash("Incorrect arguments for query provided", "error")
        return redirect(url_for("index"))

    filepath = os.path.join(os.path.join(utils.get_path("project_data", "webserver"), pro.display_key), doc.name)
    if os.path.exists(filepath):
        os.remove(filepath)

    doc.delete()
    return redirect(url_for("manage.project", oid=oid))


@manage.route("/copy_to_staging_drive/<oid>|<data_type>|<move_type>")
@manage.route("/copy_to_staging_drive/<oid>|<data_type>")
@login_required("ANY")
def copy_to_staging_drive(oid="", data_type="", move_type=""):
    if move_type == "":
        return render_template("transfer_target.html", oid=oid, data_type=data_type)

    if current_user.get_role() == "Site Admin":
        return redirect(url_for("activity"))

    if oid == "":
        flash("Could not identify the provided object.", "warning")
        return redirect(url_for("index"))

    if data_type != "pipeline_output" and data_type != "pipeline_sample_group" and data_type != "project_sample_group" and data_type != "project_pipeline_output" and data_type != "sample":
        flash("Could not identify the object type.", "warning")
        return redirect(url_for("index"))

    from biocomputedm.manage.helpers.manage_helper import copy_data_to_staging
    if move_type == "self":
        copy_data_to_staging(current_app._get_current_object(), oid, data_type, current_user.display_key)
    else:
        copy_data_to_staging(current_app._get_current_object(), oid, data_type, current_user.display_key, "yes")

    flash("Your files are now being moved to a secure location. This may take some time to complete. You will receive an email when this process is finished.", "success")
    return redirect(url_for("index"))
