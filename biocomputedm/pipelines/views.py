import os
import subprocess
from datetime import datetime

from biocomputedm import utils
from biocomputedm.admin.models import ReferenceData
from biocomputedm.decorators import login_required
from biocomputedm.manage.models import Submission, SampleGroup
from biocomputedm.pipelines import models
from biocomputedm.pipelines.forms import build_options_form, PipelinePropertiesForm
from biocomputedm.pipelines.models import Pipeline, PipelineInstance, PipelineModuleInstance, PipelineModuleOptionValue, \
    create_pipeline_instance
from flask import Blueprint, redirect, url_for, current_app, render_template
from flask import abort
from flask import flash
from flask import request
from flask import send_from_directory
from flask.ext.login import current_user
from werkzeug.utils import secure_filename

pipelines = Blueprint("pipelines", __name__, static_folder="static", template_folder="templates")


@pipelines.route("/pipelines_message/<oid>", methods=["POST"])
def message(oid=""):
    msg = request.form
    event = msg["event"]

    try:
        # Module has finished
        if event == "module_end":
            module_instance = PipelineModuleInstance.query.filter_by(display_key=oid).first()
            if module_instance is None:
                return "<html>at p w/ " + oid + "</html>"
                return abort(404)

            # Check that the pipeline isn't in blocking mode
            pipeline_instance = module_instance.pipeline_instance
            if pipeline_instance.current_execution_status != "RUNNING":
                return "<html>pipeline is in a false state</html>"
                return abort(404)

            # Display a list of module options to the user depending on the continue building flag
            # Note this is fall through from post so that we can iterate where necessary
            module_instances = pipeline_instance.module_instances.all()

            # Find the correct module template
            module = None
            for mod in module_instances:
                if mod.module.execution_index == pipeline_instance.current_execution_index:
                    module = mod
                    break

            # Check that this is our currently running module
            if module.display_key != oid:
                return "<html>old message!</html>"
                return abort(404)

            # TODO: Check padding for day
            if msg["sub"] == "0":
                sub = datetime.strptime(msg["sub"], '%a %b %d %H:%M:%S %Y')
                start = datetime.strptime(msg["start"], '%a %b %d %H:%M:%S %Y')
                end = datetime.strptime(msg["end"], '%a %b %d %H:%M:%S %Y')
                wait = start - sub
                duration = end - start
                module_instance.update(wait_time=wait, execution_time=duration)

            # If its the last module
            if pipeline_instance.current_execution_index == len(pipeline_instance.pipeline.modules.all()) - 1:
                pipeline_instance.update(current_execution_status="FINISHED")

                from biocomputedm.pipelines.helpers.pipeline_helper import finish_pipeline_instance
                finish_pipeline_instance(current_app._get_current_object(),
                                         module_instance.pipeline_instance.display_key,
                                         module_instance.pipeline_instance.data_consigner.display_key)

            # If we need to wait for options
            elif pipeline_instance.execution_type == "Per Module":
                pipeline_instance.update(current_execution_index=pipeline_instance.current_execution_index + 1,
                                         current_execution_status="WAITING")

                # TODO: Email notification

            # Otherwise continue pipeline execution
            else:
                pipeline_instance.update(current_execution_index=(pipeline_instance.current_execution_index + 1))

                from biocomputedm.pipelines.helpers.pipeline_helper import execute_module_instance
                execute_module_instance(current_app._get_current_object(),
                                        module_instance.pipeline_instance.display_key,
                                        module_instance.pipeline_instance.data_consigner.display_key)

        # Module has error
        elif event == "module_error":
            module_instance = PipelineModuleInstance.query.filter_by(display_key=oid).first()
            if module_instance is None:
                return "<html>at p with v = " + oid + " no m</html>"
                return abort(404)

            # Check that the pipeline isn't in blocking mode
            pipeline_instance = module_instance.pipeline_instance
            if pipeline_instance.current_execution_status != "RUNNING":
                return "<html>pipeline is in a false state</html>"
                return abort(404)

            # Display a list of module options to the user depending on the continue building flag
            # Note this is fall through from post so that we can iterate where necessary
            module_instances = pipeline_instance.module_instances.all()

            # Find the correct module template
            module = None
            for mod in module_instances:
                if mod.module.execution_index == pipeline_instance.current_execution_index:
                    module = mod
                    break

            # Check that this is our currently running module
            if module.display_key != oid:
                return "<html>old message!</html>"
                return abort(404)

            module_instance.pipeline_instance.update(current_execution_status="ERROR")

        # WHA
        else:
            return "<html>NOOP</html>"
            return abort(404)

        return "<html>success</html>"
        return abort(404)

    except Exception as e:
        return "<html>" + str(e) + "</html>"
        return abort(404)


@pipelines.route("/refresh_pipelines")
@login_required("Site Admin")
def refresh_pipelines():
    found = models.refresh_pipelines()
    if found:
        flash("Successfully loaded all pipelines from the pipeline directory.", "success")
    else:
        flash("No pipelines were loaded as none were found or they were already present.", "warning")

    return redirect(url_for("admin.administrate"))


@pipelines.route("/display_pipelines/")
@pipelines.route("/display_pipelines/<int:page>")
@pipelines.route("/display_pipelines/<int:display>")
@pipelines.route("/display_pipelines/<int:page>|<int:display>")
@login_required("ANY")
def display_pipelines(display=0, page=1):
    p = Pipeline.query.paginate(page=page, per_page=20)
    return render_template("pipelines.html", title="Pipelines", page=page, obs=p, display=display)


@pipelines.route("/display_pipeline/<pid>")
@pipelines.route("/display_pipeline/<pid>|<int:download>")
@login_required("ANY")
def display_pipeline(pid="", download=0):
    p = Pipeline.query.filter_by(display_key=pid).first()
    if p is None:
        flash("There was an error finding your pipeline", "error")
        return redirect(url_for("index"))

    if download == 1:
        directory = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "pipelines"), p.name)
        return send_from_directory(directory, p.documentation, as_attachment=True)

    return render_template("pipeline.html", title="Pipeline: " + p.name, pipeline=p)


@pipelines.route("/display_pipeline_instances/<int:page>")
@pipelines.route("/display_pipeline_instances")
@login_required("ANY")
def display_pipeline_instances(page=1):
    obs = current_user.group.pipeline_instances.paginate(page=page, per_page=20)
    return render_template("pipeline_instances.html", title="Pipeline Instances", page=page, obs=obs)


@pipelines.route("/display_pipeline_instance/<oid>")
@pipelines.route("/display_pipeline_instance/<oid>|<data_file>")
@login_required("ANY")
def display_pipeline_instance(oid="", data_file=""):
    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not locate the provided pipeline instance", "warning")
        return redirect(url_for("empty"))

    data_path = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key),
                             "pipeline_output")
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

    return render_template("pipeline_instance.html",
                           title="Pipeline Instance",
                           pipeline_instance=pipeline_instance,
                           data_source_type="pipeline_output",
                           oid=pipeline_instance.display_key,
                           datasets=datasets,
                           data_file="")


# @pipelines.route("/run_pipeline", methods=["GET", "POST"])
# @login_required("ANY")
# def run_pipeline():
#     form = RunPipelineForm()
#     if request.method == "GET":
#         return render_template("run_pipeline.html", form=form)
#
#     else:
#         form2 = SelectDataForPipelineForm()
#         # Build a list of uploaded files and display them for selection
#         if form.validate_on_submit() and form.submit.data:
#
#             # Build a list of uploaded files and display them for selection
#             if str(form.runtime_type.data) == "Upload":
#                 # Current user information
#                 folder = current_user.display_key
#
#                 # Build path to the users sftp dir
#                 directory_path = os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], folder)
#                 directory_path = os.path.join(directory_path, "landing_zone")
#
#                 # list of the available files
#                 filepaths = next(os.walk(directory_path))[2]
#                 files = []
#                 for file in filepaths:
#                     s = os.stat(os.path.join(directory_path, file))
#                     files.append({
#                         "name": file,
#                         "size": s.st_size,
#                         "date": time.ctime(s.st_ctime)
#                     })
#
#             # Build a list of data sources and display them for selection
#             else:
#                 items = current_user.group.data_source.all()
#
#             return render_template("run_pipeline.html", form=form, form2=form2, items=items)
#
#         elif form2.validate_on_submit() and form2.submit.data:
#



@pipelines.route("/build_pipeline_instance/<oid>|<pid>|<runtime_type>", methods=["GET", "POST"])
@login_required("ANY")
def build_pipeline_instance(oid="", pid="", runtime_type=""):
    if oid == "" or pid == "" or type == "":
        flash("The pipeline information provided was invalid", "error")
        return redirect(url_for("index"))

    # Retrieve our db records
    pipeline = Pipeline.query.filter_by(display_key=pid).first()
    if pipeline is None:
        flash("There pipeline information provided was invalid", "error")
        return redirect(url_for("index"))

    # Generate the form responsible for gathering the required information
    form = PipelinePropertiesForm(request.form)

    # Build the pipeline information and get information from the user about how we want to execute
    if request.method == "GET":
        return render_template("build_pipeline_instance.html", title="Pipeline Options", pid=pipeline.display_key,
                               oid=oid,
                               runtime_type=runtime_type, pipeline=pipeline, form=form)

    # The user has selected the pipeline and now wants to execute it using their selection
    else:
        # Get our properties
        execution_type = str(form.execution_field.data)
        options_type = str(form.options_field.data)

        # Check our other input args
        if execution_type != "Per Module" and execution_type != "Continuous":
            flash("The pipeline information provided regarding execution type was invalid", "error")
            return redirect(url_for("index"))

        if options_type != "Custom" and options_type != "Default":
            flash("The pipeline information provided regarding options type was invalid", "error")
            return redirect(url_for("index"))

        if runtime_type != "Submission" and runtime_type != "SampleGroup":
            flash("The pipeline information provided regarding runtime type was invalid", "error")
            return redirect(url_for("index"))

        # Type of datasource
        if runtime_type == "Submission":
            data_source = Submission.query.filter_by(display_key=oid).first()
        else:
            data_source = SampleGroup.query.filter_by(display_key=oid).first()

        if data_source is None:
            flash("The pipeline information provided was invalid", "error")
            return redirect(url_for("index"))

        # Check data source
        if data_source.running_pipeline is not None:
            flash(
                    "The provided data source is already attached to a pipeline. Please ensure it finishes or you quit the original pipeline before proceeding",
                    "warning")
            return redirect(url_for("index"))

        # Instance creation and assignment
        pipeline_instance = create_pipeline_instance(current_user, pipeline, data_source, execution_type,
                                                     options_type)

        # Create the directory to hold the submission
        pipeline_directory = utils.get_path("pipeline_data", "webserver")
        pipeline_directory = os.path.join(pipeline_directory, pipeline_instance.display_key)
        utils.make_directory(pipeline_directory)
        utils.make_directory(os.path.join(pipeline_directory, "samples_output"))
        modules_directory = os.path.join(pipeline_directory, "modules_output")
        utils.make_directory(modules_directory)
        utils.make_directory(os.path.join(pipeline_directory, "pipeline_output"))

        # Make the module instances and their directories
        for module in pipeline.modules:
            module_instance = PipelineModuleInstance.create(module=module, instance=pipeline_instance)
            pipeline_instance.module_instances.append(module_instance)
            pipeline_instance.save()
            utils.make_directory(os.path.join(modules_directory, module.name))

        pipeline_instance.update(current_execution_status="WAITING", current_execution_index=0)

        return redirect(url_for("pipelines.build_module_instance", pid=pipeline_instance.display_key, oid=oid, index=0))


@pipelines.route("/build_module_instance/<pid>|<oid>|<int:index>", methods=["GET", "POST"])
@login_required("ANY")
def build_module_instance(pid="", oid="", index=-1):
    if pid == "" or oid == "":
        flash("There was an error with your provided information.", "error")
        return redirect(url_for("index"))

    pipeline_instance = PipelineInstance.query.filter_by(display_key=pid).first()
    if pipeline_instance is None:
        flash("There was an error with your provided information.", "error")
        return redirect(url_for("index"))

    # We are out of modules
    module_instances = pipeline_instance.module_instances.all()
    if index >= len(module_instances):
        flash("There are no more modules to build!", "warning")
        pipeline_instance.update(current_execution_status="FINISHED")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=pipeline_instance.display_key))

    # Find the correct module template
    module_instance = None
    for mod in module_instances:
        if mod.module.execution_index == index:
            module_instance = mod
            break

    # Check that we haven't been provided with a faulty index
    if module_instance is None:
        flash("There was an error building the next pipeline module", "error")
        return redirect(url_for("index"))

    # Delete any existing options for this module
    values = module_instance.option_values
    for value in values:
        module_instance.option_values.remove(value)
        value.delete()

    # We only want to assign values for the module we are about to generate that are absolutely necessary
    options = []
    default_options = []
    if pipeline_instance.options_type == "Default":
        # Retrieve the module options
        possible_options = module_instance.module.options.all()
        for option in possible_options:
            if option.necessary:
                options.append(option)
            else:
                default_options.append(option)

    # We want to assign all values for the module we are about to generate
    else:
        # Retrieve the module options
        options = module_instance.module.options.all()

    # Handle the case where there are no options to assign to this module
    if len(options) == 0:
        for option in default_options:
            # Append the default option set
            option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
            option_value.update(value=option.default_value)
            module_instance.option_values.append(option_value)
            module_instance.save()

        # If we are out of modules to assign or we do not want to assign more information just yet
        if len(module_instances) == index + 1 or pipeline_instance.execution_type == "Per Module":
            from biocomputedm.pipelines.helpers.pipeline_helper import execute_module_instance
            execute_module_instance(current_app._get_current_object(), pid, oid)
            flash(
                    "The next pipeline step is queued for submission. It may take time before it is registered as running.",
                    "success")
            return redirect(url_for("pipelines.display_pipeline_instance", oid=pid))

        # Continue assigning information
        else:
            return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index + 1))

    # Build a form containing all of the options
    form = build_options_form(options, request.form)
    if form is None:
        flash("There was an error in the selected pipeline", "error")
        return redirect(url_for("index"))

    # Parse the options and inject them into our pipeline information
    if request.method == "POST":
        # We need to handle the case where a template has been requested
        # (i.e. a button other than the form submit button)
        file_options = []
        for option in options:
            if option.user_interaction_type == "file":
                file_options.append(option)

        # We have a request for a file template
        if len(file_options) != 0:
            # Note, this will only catch 1 - we will need to test if their data status
            # persists between the responses, if so we will get in a loop here
            for option in file_options:
                field = getattr(form, option.display_key + "_template")
                if field is not None and bool(field.data):
                    directory = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "pipelines"),
                                             pipeline_instance.pipeline.name)
                    return send_from_directory(directory, option.default_value, as_attachment=True)

        # Validate the form properties
        if form.validate_on_submit():
            # Add in any default fields?
            if pipeline_instance.options_type == "Default":
                for option in default_options:
                    # Append the default option set
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    option_value.update(value=option.default_value)
                    module_instance.option_values.append(option_value)
                    module_instance.save()

            for option in options:
                if option.user_interaction_type == "file":
                    field = getattr(form, option.display_key + "_upload")

                    # Update the db
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)

                    # Save the file - we cannot validate this so we have to hope that it is provisioned in the pipeline
                    directory = utils.get_path("pipeline_data", "webserver")
                    directory = os.path.join(directory, pipeline_instance.display_key)
                    filename = secure_filename(field.data.filename)
                    if filename == "":
                        filename = option.name + "_" + option_value.display_key + "_file"
                    filepath = os.path.join(directory, filename)
                    field.data.save(filepath)

                    # Build the filepath into our value
                    option_value.update(value=filename)
                    module_instance.option_values.append(option_value)
                    module_instance.save()

                elif option.user_interaction_type == "string" \
                        or option.user_interaction_type == "enum":
                    field = getattr(form, option.display_key)

                    # Update the db
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    option_value.update(value=str(field.data))
                    module_instance.option_values.append(option_value)
                    module_instance.save()

                elif option.user_interaction_type == "boolean":
                    field = getattr(form, option.display_key)

                    # Update the db
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    option_value.update(value=bool(field.data))
                    module_instance.option_values.append(option_value)
                    module_instance.save()

                elif option.user_interaction_type == "reference":
                    field = getattr(form, option.display_key)

                    # Inefficient but effective :P
                    value = str(field.data)
                    name = value.split(" (")[0]

                    # Lookup
                    reference_data = ReferenceData.query.filter_by(display_key=name).first()
                    path = os.path.join(os.path.join(utils.get_path("reference_data", "hpc"), reference_data.path), reference_data.name)

                    # Update the db
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    option_value.update(value=path)
                    module_instance.option_values.append(option_value)
                    module_instance.save()

                else:
                    continue

            # If we are out of modules to assign or we do not want to assign more information just yet
            if len(module_instances) == index + 1 or pipeline_instance.execution_type == "Per Module":
                from biocomputedm.pipelines.helpers.pipeline_helper import execute_module_instance
                execute_module_instance(current_app._get_current_object(), pid, oid)
                flash(
                        "The next pipeline step is queued for submission. It may take time before it is registered as running.",
                        "success")
                return redirect(url_for("pipelines.display_pipeline_instance", oid=pid))

            # Continue assigning information
            else:
                return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index + 1))

        else:
            utils.flash_errors(form)

    return render_template("build_module_instance.html", title="Pipeline Module Options", pid=pid, oid=oid, index=index,
                           form=form)


@pipelines.route("/module_instance/<oid>")
@pipelines.route("/module_instance/<oid>|<data_file>")
@login_required("ANY")
def module_instance(oid="", data_file=""):
    if oid == "":
        flash("No instance identifiers were provided!", "warning")
        return redirect(url_for("empty"))

    pipeline_instances = current_user.group.pipeline_instances.all()
    p_instance = None
    m_instance = None
    for pipeline_instance in pipeline_instances:
        module_instances = pipeline_instance.module_instances.all()
        for module_instance in module_instances:
            if module_instance.display_key == oid:
                m_instance = module_instance
                p_instance = pipeline_instance
                break

        if m_instance is not None:
            break

    if m_instance is None:
        flash("Could not locate the provided module instance", "warning")
        return redirect(url_for("empty"))

    data_path = os.path.join(
            os.path.join(
                    os.path.join(
                            utils.get_path("pipeline_data", "webserver"),
                            p_instance.display_key),
                    "modules_output"),
            m_instance.module.name)

    datasets = []
    if os.path.isdir(data_path):
        filepaths = next(os.walk(data_path))
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

    return render_template("module_instance.html",
                           title="Module Instance",
                           module_instance=m_instance,
                           data_source_type="module_output",
                           oid=m_instance.display_key,
                           datasets=datasets,
                           data_file="")


@pipelines.route("/continue_pipeline_instance/<oid>")
@login_required("ANY")
def continue_pipeline(oid=""):
    if oid == "":
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    # We were not waiting for input - something else has happened
    if pipeline_instance.current_execution_status != "WAITING":
        flash("Could not resume the current pipeline as it is in a false state", "error")
        return redirect(url_for("empty"))

    # We have the options already - from a state perspective I am not sure how this will arise but its best to handle
    if pipeline_instance.execution_type == "Continuous":
        from biocomputedm.pipelines.helpers.pipeline_helper import execute_module_instance
        execute_module_instance(current_app._get_current_object(),
                                pipeline_instance.display_key,
                                pipeline_instance.data_consigner.display_key)

        flash("Submitting the next module for execution!", "success")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

    # Need to obtain the options
    else:
        return redirect(url_for("pipelines.build_module_instance",
                                pid=oid,
                                oid=pipeline_instance.data_consigner.display_key,
                                index=pipeline_instance.current_execution_index))


@pipelines.route("/change_module/<oid>|<change_type>")
@pipelines.route("/change_module/<oid>|<change_type>|<int:force>")
@login_required("ANY")
def change_module(oid="", change_type="", force=0):
    if force != 1:
        return render_template("confirm.html",
                               message="Are you sure you wish to restart the current module?",
                               oid=oid,
                               url="pipelines.change_module",
                               type=change_type)

    if oid == "":
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

    # TODO - If module is running parse for job id and kill all

    if change_type == "back":
        if pipeline_instance.current_execution_index == 0:
            flash("Cannot go back a step when there are no previous steps!", "error")
            return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

        pipeline_instance.update(current_execution_status="STOPPED")

        # Find the correct module template
        module_instances = pipeline_instance.module_instances.all()
        module = None
        for mod in module_instances:
            if mod.module.execution_index == pipeline_instance.current_execution_index:
                module = mod
                break

        # Clean the current module directory
        subprocess.Popen(
                [
                    "sudo",
                    os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "cleanup"), "wipe_directory.sh"),
                    "-p=" + os.path.join(
                            os.path.join(
                                    os.path.join(
                                            utils.get_path("pipeline_data", "webserver"),
                                            pipeline_instance.display_key),
                                    "modules_output"),
                            module.module.name)
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
        )

        target_index = pipeline_instance.current_execution_index - 1

    elif change_type == "current":
        pipeline_instance.update(current_execution_status="STOPPED")

        target_index = pipeline_instance.current_execution_index

    elif change_type == "next":
        if pipeline_instance.current_execution_index >= len(pipeline_instance.module_instances.all()) - 1:
            flash("Cannot proceed to the next step when there are no more steps!", "error")
            return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

        pipeline_instance.update(current_execution_status="STOPPED")

        target_index = pipeline_instance.current_execution_index + 1

    else:
        flash("Could not modify the provided pipeline instance", "error")
        return redirect(url_for("index"))

    # Find the correct module template
    module_instances = pipeline_instance.module_instances.all()
    module = None
    for mod in module_instances:
        if mod.module.execution_index == target_index:
            module = mod
            break

    # Clean the module directory
    subprocess.Popen(
            [
                "sudo",
                os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "cleanup"), "wipe_directory.sh"),
                "-p=" + os.path.join(
                        os.path.join(
                                os.path.join(
                                        utils.get_path("pipeline_data", "webserver"),
                                        pipeline_instance.display_key),
                                "modules_output"),
                        module.module.name)
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
    )

    # Change the current execution index
    pipeline_instance.update(current_execution_index=target_index)

    # We have the options already
    if pipeline_instance.execution_type == "Continuous":
        from biocomputedm.pipelines.helpers.pipeline_helper import execute_module_instance
        execute_module_instance(current_app._get_current_object(),
                                pipeline_instance.display_key,
                                pipeline_instance.data_consigner.display_key)

        flash("Submitting the next module for execution!", "success")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

    # Need to obtain the options
    else:
        flash("Please enter the options for this module.", "success")
        return redirect(
                url_for("pipelines.build_module_instance",
                        pid=oid,
                        oid=pipeline_instance.data_consigner.display_key,
                        index=pipeline_instance.current_execution_index))


@pipelines.route("/finish_pipeline/<oid>")
@pipelines.route("/finish_pipeline/<oid>|<int:force>")
@login_required("ANY")
def finish_pipeline(oid="", force=0):
    if force != 1:
        return render_template("confirm.html",
                               message="Are you sure you wish to quit the current pipeline?",
                               oid=oid,
                               url="pipelines.finish_pipeline")

    if oid == "":
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    pipeline_instance.update(current_execution_status="STOPPED")

    # TODO - If module is running parse for job id and kill all

    current_data_source = pipeline_instance.data_consigner
    current_data_source.run_pipelines.append(pipeline_instance)
    current_data_source.update(currently_running_pipeline=None)

    flash("The pipeline was stopped and disassociated with your parent data set", "success")
    return redirect(url_for("index"))


@pipelines.route("/restart_pipeline/<oid>")
@pipelines.route("/restart_pipeline/<oid>|<int:force>")
@login_required("ANY")
def restart_pipeline(oid="", force=0):
    if force != 1:
        return render_template("confirm.html",
                               message="Are you sure you wish to restart the current pipeline?",
                               oid=oid,
                               url="pipelines.restart_pipeline")

    if oid == "":
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    pipeline_instance.update(current_execution_status="STOPPED")

    # TODO - If module is running parse for job id and kill all

    data_source = pipeline_instance.data_consigner
    data_source.run_pipelines.append(pipeline_instance)
    data_source.update(running_pipeline=None)

    flash("The previous pipeline has been removed, follow the instructions below to restart!", "success")
    return redirect(url_for("pipelines.build_pipeline_instance",
                            pid=pipeline_instance.pipeline.display_key,
                            oid=data_source.display_key,
                            runtime_type=data_source.type))
