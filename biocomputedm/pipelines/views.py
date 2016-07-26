import os
import subprocess
from datetime import datetime

from biocomputedm import utils
from biocomputedm.admin.models import ReferenceData
from biocomputedm.decorators import login_required
from biocomputedm.manage.models import DataGroup
from biocomputedm.pipelines import models
from biocomputedm.pipelines.forms import build_options_form, PipelinePropertiesForm
from biocomputedm.pipelines.models import Pipeline, PipelineInstance, PipelineModuleInstance, PipelineModuleOptionValue
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

            # Find the correct module template
            module_instances = pipeline_instance.module_instances.all()
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
                pipeline_instance.update(current_execution_status="STOPPED")

                from biocomputedm.pipelines.helpers.pipeline_helper import finish_pipeline_instance
                finish_pipeline_instance(current_app._get_current_object(), pipeline_instance.display_key)

            else:
                pipeline_instance.update(current_execution_index=pipeline_instance.current_execution_index + 1, current_execution_status="WAITING")

                if pipeline_instance.execution_type == "Continuous":
                    from biocomputedm.pipelines.helpers.pipeline_helper import execute_pipeline_module
                    execute_pipeline_module(current_app._get_current_object(), pipeline_instance.display_key)

        # Module has error
        elif event == "module_error":
            module_instance = PipelineModuleInstance.query.filter_by(display_key=oid).first()
            if module_instance is None:
                return "<html>at p with v = " + oid + " no m</html>"
                return abort(404)

            # Check that the pipeline isn't in blocking mode
            pipeline_instance = module_instance.pipeline_instance
            # if pipeline_instance.current_execution_status != "RUNNING":
            #     return "<html>pipeline is in a false state</html>"
            #     return abort(404)

            # Find the correct module template
            module_instances = pipeline_instance.module_instances.all()
            module = None
            for mod in module_instances:
                if mod.module.execution_index == pipeline_instance.current_execution_index:
                    module = mod
                    break

            # Check that this is our currently running module
            if module.display_key != oid:
                return "<html>old message!</html>"
                return abort(404)

            from biocomputedm.pipelines.helpers.pipeline_helper import parse_outputs
            parse_outputs(current_app._get_current_object(), pipeline_instance.display_key)

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
    if current_user.get_role() == "Site Admin":
        obs = PipelineInstance.query.paginate(page=page, per_page=20)
    else:
        obs = current_user.group.pipeline_instances.paginate(page=page, per_page=20)

    return render_template("pipeline_instances.html", title="Pipeline Instances", page=page, obs=obs)


@pipelines.route("/display_pipeline_instance/<oid>")
@login_required("ANY")
def display_pipeline_instance(oid=""):
    if oid == "":
        flash("Could not identify the provided pipeline run", "warning")
        return redirect(url_for("index"))

    if current_user.get_role() == "Site Admin":
        pipeline_instance = PipelineInstance.query.filter_by(display_key=oid).first()
    else:
        pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()

    if pipeline_instance is None:
        flash("Could not locate the provided pipeline instance", "warning")
        return redirect(url_for("empty"))

    return render_template("pipeline_instance.html", title="Pipeline Instance", pipeline_instance=pipeline_instance)


@pipelines.route("/build_pipeline_instance/<oid>|<pid>", methods=["GET", "POST"])
@login_required("ANY")
def build_pipeline_instance(oid="", pid=""):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    if oid == "" or pid == "":
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
        # Pipeline defaults display
        options = []
        for module in pipeline.modules:
            for option in module.options:
                options.append(option)


        return render_template(
            "build_pipeline_instance.html",
            title="Pipeline Options",
            pid=pipeline.display_key,
            oid=oid,
            pipeline=pipeline,
            options=options,
            form=form
        )

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

        # Get the data group to use as a source
        source_data_group = DataGroup.query.filter_by(display_key=oid).first()
        if source_data_group is None:
            flash("The pipeline information provided was invalid", "error")
            return redirect(url_for("index"))

        # Instance creation and assignment
        pipeline_instance = PipelineInstance.create(pipeline=pipeline, execution_type=execution_type, options_type=options_type, user=current_user, data_consignor=source_data_group)

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

        from biocomputedm.pipelines.helpers.pipeline_helper import initialise_running_pipeline
        initialise_running_pipeline(pipeline_instance.display_key, source_data_group.display_key)

        return redirect(url_for("pipelines.build_module_instance", pid=pipeline_instance.display_key, oid=oid, index=0))


@pipelines.route("/build_module_instance/<pid>|<oid>|<int:index>", methods=["GET", "POST"])
@login_required("ANY")
def build_module_instance(pid="", oid="", index=-1):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    if pid == "" or oid == "":
        flash("There was an error with your provided information.", "error")
        return redirect(url_for("index"))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=pid).first()
    if pipeline_instance is None:
        flash("There was an error with your provided information.", "error")
        return redirect(url_for("index"))

    # We are out of modules
    module_instances = pipeline_instance.module_instances.all()
    if index >= len(module_instances):
        flash("There are no more modules to build.", "warning")
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

            # Special case enums as we should take the first value here!
            if option.user_interaction_type == "enum":
                default_value = option.default_value.split(",")[0]

            else:
                default_value = option.default_value

            option_value.update(value=default_value)
            module_instance.option_values.append(option_value)
            module_instance.save()

        # If we are out of modules to assign or we do not want to assign more information just yet
        if len(module_instances) == index + 1 or pipeline_instance.execution_type == "Per Module":
            from biocomputedm.pipelines.helpers.pipeline_helper import execute_pipeline_module
            execute_pipeline_module(current_app._get_current_object(), pid)
            flash(
                "The next pipeline step is queued for submission. It may take time before it is registered as running.",
                "success"
            )

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
                    directory = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "pipelines"), pipeline_instance.pipeline.name)
                    return send_from_directory(directory, option.default_value, as_attachment=True)

        # Validate the form properties
        if form.validate_on_submit():
            # Add in any default fields?
            if pipeline_instance.options_type == "Default":
                for option in default_options:
                    # Append the default option set
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)

                    # Special case enums as we should take the first value here!
                    if option.user_interaction_type == "enum":
                        default_value = option.default_value.split(",")[0]

                    else:
                        default_value = option.default_value

                    option_value.update(value=default_value)
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
                    value = bool(field.data)
                    option_value.update(value=str(value))
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
                from biocomputedm.pipelines.helpers.pipeline_helper import execute_pipeline_module
                execute_pipeline_module(current_app._get_current_object(), pid)
                flash(
                    "The next pipeline step is queued for submission. It may take time before it is registered as running.",
                    "success"
                )
                return redirect(url_for("pipelines.display_pipeline_instance", oid=pid))

            # Continue assigning information
            else:
                return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index + 1))

        else:
            utils.flash_errors(form)

    return render_template("build_module_instance.html", title="Pipeline Module Options", pid=pid, oid=oid, index=index,
                           form=form)


@pipelines.route("/module_instance/<pid>|<oid>")
@login_required("ANY")
def module_instance(pid="", oid=""):
    if pid == "" or oid == "":
        flash("No instance identifiers were provided.", "warning")
        return redirect(url_for("empty"))

    if current_user.get_role() == "Site Admin":
        pipeline_instance = PipelineInstance.query.filter_by(display_key=pid).first()
    else:
        pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=pid).first()

    if pipeline_instance is None:
        flash("Could not identify the module parent", "warning")
        return redirect(url_for("index"))

    m_instance = None
    module_instances = pipeline_instance.module_instances.all()
    for module_instance in module_instances:
        if module_instance.display_key == oid:
            m_instance = module_instance
            break

    if m_instance is None:
        flash("Could not locate the provided module instance", "warning")
        return redirect(url_for("empty"))

    # Conditionally index the module files so that we can view them on the web page whilst executing - if the module is done (for whatever reason) instead, make use of the indexed files
    files = []
    if pipeline_instance.current_execution_status != "FINISHED" and pipeline_instance.current_execution_status != "STOPPED" and pipeline_instance.current_execution_status != "ERROR":
        local_pipeline_directory = os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key)
        local_module_directory = os.path.join(os.path.join(local_pipeline_directory, "modules_output"), m_instance.module.name)
        filepaths = next(os.walk(local_module_directory))
        for file in filepaths[1]:
            path = os.path.join(os.path.join(os.path.join(pipeline_instance.display_key, "modules_output"), m_instance.display_key), file)
            files.append({
                "name": file,
                "path": path
            })

        for file in filepaths[2]:
            path = os.path.join(local_module_directory, file)
            files.append({
                "name": file,
                "path": path
            })

    return render_template("module_instance.html", title="Module Instance", module_instance=m_instance, files=files)


@pipelines.route("/delete_pipeline_instance/<oid>")
@pipelines.route("/delete_pipeline_instance/<oid>|<int:force>")
@login_required("ANY")
def delete_pipeline_instance(oid="", force=0):
    if force != 1:
        return render_template(
            "confirm.html",
            message="Are you sure you wish to remove this pipeline run?",
            oid=oid,
            url="pipelines.delete_pipeline_instance"
        )

    if oid == "":
        flash("No instance identifiers were provided.", "warning")
        return redirect(url_for("empty"))

    if current_user.get_role() == "Site Admin":
        pipeline_instance = PipelineInstance.query.filter_by(display_key=oid).first()
    else:
        pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()

    if pipeline_instance is None:
        flash("Could not identify the module parent", "warning")
        return redirect(url_for("index"))

    data_group = pipeline_instance.pipeline_output
    pipeline_instance.delete()
    if data_group is not None:
        data_group.delete()

    subprocess.Popen(
        [
            "sudo",
            os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "io"), "delete.sh"),
            "-s=" + os.path.join(utils.get_path("pipeline_data", "webserver"), oid)
        ]
    )

    flash("Pipeline was deleted successfully", "success")
    return redirect(url_for("pipelines.display_pipeline_instances"))


@pipelines.route("/continue_pipeline_instance/<oid>")
@login_required("ANY")
def continue_pipeline(oid=""):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

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
        from biocomputedm.pipelines.helpers.pipeline_helper import execute_pipeline_module
        execute_pipeline_module(current_app._get_current_object(), pipeline_instance.display_key)

        flash("Submitting the next module for execution.", "success")
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
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

    if force != 1:
        return render_template("confirm.html",
                               message="Are you sure you wish to restart the current module?",
                               oid=oid,
                               url="pipelines.change_module",
                               type=change_type)

    if oid == "":
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("pipelines.display_pipelines"))

    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=oid).first()
    if pipeline_instance is None:
        flash("Could not load the provided pipeline instance", "error")
        return redirect(url_for("pipelines.display_pipeline_instance", oid=oid))

    # TODO - If module is running parse for job id and kill all

    if change_type == "back":
        if pipeline_instance.current_execution_index == 0:
            flash("Cannot go back a step when there are no previous steps.", "error")
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
            flash("Cannot proceed to the next step when there are no more steps.", "error")
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

        flash("Submitting the next module for execution.", "success")
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
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

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

    flash("The pipeline was stopped and disassociated with your parent data set", "success")
    return redirect(url_for("index"))


@pipelines.route("/restart_pipeline/<oid>")
@pipelines.route("/restart_pipeline/<oid>|<int:force>")
@login_required("ANY")
def restart_pipeline(oid="", force=0):
    if current_user.get_role() == "Site Admin":
        return redirect(url_for("content.activity"))

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

    data_source = pipeline_instance.data_consignor
    data_source.run_pipelines.append(pipeline_instance)
    data_source.update(running_pipeline=None)

    flash("The previous pipeline has been removed, follow the instructions below to restart.", "success")
    return redirect(url_for("pipelines.build_pipeline_instance",
                            pid=pipeline_instance.pipeline.display_key,
                            oid=data_source.display_key,
                            runtime_type=data_source.type))
