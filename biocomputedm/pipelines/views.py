import os
import subprocess

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.manage.models import Submission, SampleGroup
from biocomputedm.pipelines import models
from biocomputedm.pipelines.forms import build_options_form, PipelinePropertiesForm
from biocomputedm.pipelines.models import Pipeline, PipelineInstance, PipelineModuleInstance, PipelineModuleOptionValue, \
    create_pipeline_instance
from flask import Blueprint, redirect, url_for, current_app, render_template
from flask import flash
from flask import request
from flask import send_from_directory
from flask.ext.login import current_user
from werkzeug.utils import secure_filename

pipelines = Blueprint("pipelines", __name__, static_folder="static", template_folder="templates")


@pipelines.route("/refresh_pipelines")
@login_required("Site Admin")
def refresh_pipelines():
    found = models.refresh_pipelines()
    if found:
        flash("Successfully loaded all pipelines from the pipeline directory.", "success")
    else:
        flash("No pipelines were loaded as none were found or they were already present.", "warning")

    return redirect(url_for("admin.administrate"))


@pipelines.route("/display_pipelines")
@pipelines.route("/display_pipelines/<int:page>")
@login_required("Site Admin", "Group Admin")
def display_pipelines(page=1):
    p = Pipeline.query.paginate(page=page, per_page=20)
    return render_template("pipelines.html", title="Pipelines", page=page, obs=p)


@pipelines.route("/display_pipeline/<pid>")
@login_required("ANY")
def display_pipeline(pid=""):
    p = Pipeline.query.filter_by(display_key=pid).first()
    if p is None:
        flash("There was an error finding your pipeline", "error")
        return redirect(url_for("index"))

    return render_template("pipeline.html", pipeline=p)


@pipelines.route("/display_pipeline_instances/<int:page>")
@pipelines.route("/display_pipeline_instances")
@login_required("ANY")
def display_pipeline_instances(page=1):
    obs = current_user.group.pipeline_instances.query.paginate(page=page, per_page=20)
    return render_template("pipeline_instances.html", title="Pipeline Instances", page=page, obs=obs)


@pipelines.route("/display_pipeline_instance/<pid>")
@login_required("ANY")
def display_pipeline_instance(pid=""):
    pipeline_instance = current_user.group.pipeline_instances.filter_by(display_key=pid).first()
    if pipeline_instance is None:
        flash("Could not locate the provided pipeline instance", "error")
        return redirect(url_for("empty"))

    return render_template("pipeline_instance.html", title="Pipeline Instance", pipeline_instance=pipeline_instance)


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

        if runtime_type != "submission" and runtime_type != "sample_group":
            flash("The pipeline information provided regarding runtime type was invalid", "error")
            return redirect(url_for("index"))

        # Type of datasource
        if runtime_type == "submission":
            data_source = Submission.query.filter_by(display_key=oid).first()
        else:
            data_source = SampleGroup.query.filter_by(display_key=oid).first()

        if data_source is None:
            flash("The pipeline information provided was invalid", "error")
            return redirect(url_for("index"))

        # Instance creation and assignment
        pipeline_instance = create_pipeline_instance(current_user.group, pipeline, data_source, execution_type,
                                                     options_type)

        # Create the directory to hold the submission
        pipeline_directory = utils.get_path("pipeline_data", "webserver")
        pipeline_directory = os.path.join(pipeline_directory, pipeline_instance.display_key)
        utils.make_directory(pipeline_directory)
        utils.make_directory(os.path.join(pipeline_directory, "samples_output"))

        return redirect(url_for("pipelines.build_module_instance", pid=pipeline_instance.display_key, oid=oid, index=0))


# The behaviour of continue building is complex an d
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

    # Display a list of module options to the user depending on the continue building flag
    # Note this is fall through from post so that we can iterate where necessary
    module_instances = pipeline_instance.module_instances.all()
    modules = pipeline_instance.pipeline.modules.all()

    # Find the correct module template
    module = None
    for mod in modules:
        if mod.execution_index == index:
            module = mod
            break

    # Check that we haven't been provided with a faulty index
    if module is None:
        flash("There was an error building the next pipeline module", "error")
        return redirect(url_for("index"))

    # Check whether we need to create a module instance for this
    if len(module_instances) == index:
        module_instance = PipelineModuleInstance.create(module=module, instance=pipeline_instance)
        pipeline_instance.module_instances.append(module_instance)
        pipeline_instance.save()

    # This is likely to occur during a page refresh
    elif len(module_instances) == index + 1:
        module_instance = PipelineModuleInstance.query.filter_by(module=module).first()
        if module_instance is None:
            flash("There was an error building the next pipeline module", "error")
            return redirect(url_for("index"))

    # Bad request
    else:
        flash("There was an error building the next pipeline module", "error")
        return redirect(url_for("index"))

    # We only want to assign values for the module we are about to generate that are absolutely necessary
    if pipeline_instance.options_type == "Default":
        # Retrieve the module options
        possible_options = module.options.all()
        options = []
        for option in possible_options:
            if option.necessary:
                options.append(option)
                continue

    # We want to assign all values for the module we are about to generate
    else:
        # Retrieve the module options
        options = module.options.all()

    # Handle the case where there are no options to assign to this module
    if len(options) == 0:
        index += 1

        # If we are out of modules to assign or we do not want to assign more information just yet
        if len(modules) == index or pipeline_instance.execution_type == "Per Module":
            pipeline_instance.update(current_execution_status="NOT_STARTED", current_execution_index=0)
            return redirect(url_for("pipelines.execute_pipeline_instance", pid=pid, oid=oid))

        # Continue assigning information
        else:
            return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index))

    # Build a form containing all of the options
    form = build_options_form(options, request.form)
    if form is None:
        flash("There was an error in the selected pipeline", "error")
        return redirect(url_for("index"))

    # Parse the options and inject them into our pipeline information
    if request.method == "POST":
        # We need to handle page refreshes
        values = module_instance.option_values.all()
        for value in values:
            value.delete()

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
                    directory = os.path.join(utils.get_path("pipeline_scripts", "webserver"),
                                             pipeline_instance.pipeline.name)
                    return send_from_directory(directory, option.default_value, as_attachment=True)

        # Validate the form properties
        if form.validate_on_submit():
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

                elif option.user_interaction_type == "string":
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

                elif option.user_interaction_type == "library":
                    field = getattr(form, option.display_key)

                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    # TODO - this should be a db ref instead but we havent implemented libraries yet
                    option_value.update(value=str(field.data))
                    module_instance.option_values.append(option_value)
                    module_instance.save()

                else:
                    continue

            # Add in any default fields?
            if pipeline_instance.options_type == "Default":
                # Retrieve the module options
                possible_options = module.options.all()
                for option in possible_options:
                    # Check that we wont have created this option already!
                    if option.necessary:
                        continue

                    # Append the default option set
                    option_value = PipelineModuleOptionValue.create(option=option, module_instance=module_instance)
                    option_value.update(value=option.default_value)
                    module_instance.option_values.append(option_value)
                    module_instance.save()

            index += 1

            # If we are out of modules to assign or we do not want to assign more information just yet
            if len(modules) == index or pipeline_instance.execution_type == "Per Module":
                pipeline_instance.update(current_execution_status="NOT_STARTED", current_execution_index=0)
                return redirect(url_for("pipelines.execute_pipeline_instance", pid=pid, oid=oid))

            # Continue assigning information
            else:
                return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index))

        else:
            utils.flash_errors(form)

    return render_template("build_module_instance.html", title="Pipeline Module Options", pid=pid, oid=oid, index=index,
                           form=form)


@pipelines.route("/resume_pipeline/<pid>")
@login_required("ANY")
def resume_pipeline_instance(pid="", oid=""):
    return redirect(url_for("empty"))


@pipelines.route("/execute_pipeline_instance/<pid>|<oid>")
@login_required("ANY")
def execute_pipeline_instance(pid="", oid=""):
    # Check that our objects exist
    pipeline_instance = PipelineInstance.query.filter_by(display_key=pid).first()
    if pipeline_instance is None:
        flash("Could not identify the provided pipeline", "error")
        return redirect(url_for("index"))

    # Directories
    local_working_directory = os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key)
    working_directory = os.path.join(utils.get_path("pipeline_data", "hpc"), pipeline_instance.display_key)
    local_csv_path = os.path.join(local_working_directory, "data_map.csv")
    csv_path = os.path.join(working_directory, "data_map.csv")
    output_directory = os.path.join(working_directory, "samples_output")

    # Get the object and build it's data path file
    import csv
    if pipeline_instance.pipeline.type == "I":
        o = Submission.query.filter_by(display_key=oid).first()
        if o is None:
            flash("Could not identify the provided object", "error")
            return redirect(url_for("index"))

        # Build the csv
        with open(local_csv_path, "a", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([oid, os.path.join(utils.get_path("submission_data", "hpc"), oid), output_directory])

    else:
        o = SampleGroup.query.filter_by(display_key=oid).first()
        if o is None:
            flash("Could not identify the provided object", "error")
            return redirect(url_for("index"))

        # Build the csv
        with open(local_csv_path, "a", newline="") as csvfile:
            writer = csv.writer(csvfile)
            for sample in o.samples.query.all():
                writer.writerow(
                        [sample.display_key, os.path.join(utils.get_path("sample_data", "hpc"), sample.display_key),
                         os.path.join(output_directory, sample.display_key)])

    current_module_instance = models.get_current_module_instance(pipeline_instance)
    if current_module_instance is None:
        flash("Could not identify the current pipeline status", "error")
        return redirect(url_for("index"))

    # Behaviour on current module status
    if pipeline_instance.current_execution_status != "NOT_STARTED":
        flash("Attempted to start a pipeline that is not in the correct state.", "error")
        return redirect(url_for("empty"))

    # Build variables string
    vstring = "'"
    for value in current_module_instance.option_values:
        marker = value.option.parameter_name
        result = value.value
        vstring += marker + "=" + result + ","
    vstring = vstring[:-1] + "'"

    # Submit module w/ options into HPC - note the cwd
    shell_path = current_app.config["HPC_JOB_SUBMISSION_FILE"]
    executor_path = current_module_instance.module.executor
    with open(os.path.join(local_working_directory, "submission_out.log"), "wb") as out, \
            open(os.path.join(local_working_directory, "submission_error.log"), "wb") as err:
        subprocess.Popen(
                [
                    shell_path,
                    "-t=" + pipeline_instance.display_key,
                    "-s=" + executor_path,
                    "-l=" + local_working_directory,
                    "-w=" + working_directory,
                    "-i=" + csv_path,
                    "-v=" + vstring
                ],
                # cwd=os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key),
                stdout=out,
                stderr=err
        )

    pipeline_instance.update(current_execution_status="RUNNING")
    flash("Pipeline executing", "Success")
    return redirect(url_for("pipelines.display_pipeline_instance", pid=pid))
