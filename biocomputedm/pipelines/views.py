import os

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.manage.models import Submission, SampleGroup
from biocomputedm.pipelines.forms import build_options_form, PipelinePropertiesForm
from biocomputedm.pipelines.models import Pipeline, PipelineInstance, PipelineModuleInstance, PipelineModuleOptionValue
from flask import Blueprint, redirect, url_for, current_app, render_template
from flask import flash
from flask import request
from flask import send_from_directory
from werkzeug.utils import secure_filename

pipelines = Blueprint("pipelines", __name__, static_folder="static", template_folder="templates")


@pipelines.route("/show_pipelines")
@pipelines.route("/show_pipelines/<int:page>")
@login_required("Site Admin", "Group Admin")
def show_pipelines(page=1):
    p = Pipeline.query.paginate(page=page, per_page=20)
    return render_template("pipelines.html", title="Pipelines", page=page, obs=p)


@pipelines.route("/refresh_pipelines")
@login_required("Site Admin")
def refresh_pipelines():
    found = utils.refresh_pipelines()
    if found:
        flash("Successfully loaded all pipelines from the pipeline directory.", "success")
    else:
        flash("No pipelines were loaded as none were found or they were already present.", "warning")

    return redirect(url_for("admin.administrate"))


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
        return render_template("build_pipeline_instance.html", pid=pipeline.display_key, oid=oid,
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
        pipeline_instance = PipelineInstance.create(pipeline=pipeline, data_source=data_source)
        pipeline_instance.update(execution_type=execution_type, options_type=options_type)

        # Create the directory to hold the submission
        pipeline_directory = utils.get_path("pipeline_data", "webserver")
        pipeline_directory = os.path.join(pipeline_directory, pipeline_instance.display_key)
        utils.make_directory(pipeline_directory)

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
            return redirect(url_for("pipelines.run_pipeline_instance", pid=pid, oid=oid))

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
                    option_value.update(value=filepath)
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
                return redirect(url_for("pipelines.run_pipeline_instance", pid=pid, oid=oid))

            # Continue assigning information
            else:
                return redirect(url_for("pipelines.build_module_instance", pid=pid, oid=oid, index=index))

        else:
            utils.flash_errors(form)
            return render_template("build_module_instance.html", pid=pid, oid=oid, index=index)

    else:
        return render_template("build_module_instance.html", pid=pid, oid=oid, index=index, form=form)


@pipelines.route("/run_pipeline_instance")
@login_required("ANY")
def run_pipeline_instance():
    # Check that the current user owns the object provided

    # Pipeline initiation db

    # Pipeline initiation variables

    # Pipeline asynchronous executor
    # TODO: set execution index
    return redirect(url_for("empty"))


@pipelines.route("/test_pipeline")
@login_required("Site Admin")
def test_pipeline(pid="", unique_group_id=""):
    app = current_app

    # TODO REMOVE BELOW THIS: ITS TEST CODE
    p = Pipeline.query.filter_by(name="test_pipeline").first()
    pid = p.display_key
    # TODO REMOVE ABOVE THIS: ITS TEST CODE

    # Get the relevant pipeline and build an instance
    p = Pipeline.query.filter_by(display_key=pid).first()
    pi = utils.create_pipeline_instance(p)

    # Get the first module and build an instance
    m = p.module.filter_by(execution_index=0).first()
    mi = utils.create_module_instance(m, pi)

    # TODO: Generate module options, based on behaviours and query the server about it

    # Setup samples csv
    import csv
    if p.type == "I":
        # Create the working directory
        working_directory = os.path.join(app.config["SUBMISSIONS_PATH_ON_HPC"], pi.display_key)
        if not os.path.exists(working_directory):
            os.makedirs(working_directory)

        # Build a single line csv that just gives information about the submission
        csv_path = os.path.join(working_directory, "io.txt")
        with open(csv_path, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([unique_group_id, app.config["RAW_DATA_PATH_ON_HPC"] + "/" + unique_group_id,
                             app.config["SUBMISSIONS_PATH_ON_HPC"] + "/" + pi.display_key])

    elif type == "II":
        # Create the working directory
        working_directory = os.path.join(app.config["SAMPLES_PATH_ON_HPC"], pi.display_key)
        if not os.path.exists(working_directory):
            os.makedirs(working_directory)

        # Build a csv of relevant samples and their output directory
        # TODO: assess for per module
        from biocomputedm.models import SampleGroup
        sg = SampleGroup.query.filter_by(display_key=unique_group_id).first()
        csv_path = os.path.join(working_directory, "io.txt")
        with open(csv_path, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            for sample in sg.sample.all():
                sample_directory = os.path.join(working_directory, sample.display_key)
                if not os.path.exists(sample_directory):
                    os.makedirs(sample_directory)

                writer.writerow([sample.display_key,
                                 app.config[
                                     "SUBMISSIONS_PATH_ON_HPC"] + "/" + sg.source_pipeline.display_key + "/" + sample.display_key,
                                 app.config["SAMPLES_PATH_ON_HPC"] + "/" + pi.display_key + "/" + sample.display_key])

    else:
        # Create the working directory
        working_directory = os.path.join(app.config["INVESTIGATIONS_PATH_ON_HPC"], pi.display_key)
        if not os.path.exists(working_directory):
            os.makedirs(working_directory)

        # Build a csv of relevant samples and their output directory
        # TODO: As above
        from biocomputedm.models import Investigation
        i = Investigation.query.filter_by(display_key=unique_group_id).first()
        csv_path = os.path.join(working_directory, "io.txt")
        with open(csv_path, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            for sg in i.sample_group.all():
                for sample in sg.sample.all():
                    sample_directory = os.path.join(working_directory, sample.display_key)
                    if not os.path.exists(sample_directory):
                        os.makedirs(sample_directory)

                    writer.writerow([sample.display_key,
                                     app.config[
                                         "SAMPLES_PATH_ON_HPC"] + "/" + sg.source_pipeline.display_key + "/" + sample.display_key,
                                     app.config[
                                         "INVESTIGATIONS_PATH_ON_HPC"] + "/" + pi.display_key + "/" + sample.display_key])

    # TODO: Assess whether or not we need to modify r/w/x privileges here
    raise NameError("Blind stop")

    # Step 4) Submit job
    shell_path = app.config["HPC_JOB_SUBMISSION_FILE"]
    pipeline_path = m.executor
    process = subprocess.Popen(
            [
                shell_path,
                "-t=A_Ticket",
                "-j=A_JOB",
                "-s=" + pipeline_path,
                "-w=" + working_directory,
                "-i=" + csv_path,
                "-v='a=1,b=eleven'"
            ],
            cwd=working_directory,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE
    ).stdout

    # Lets read out our process information - this is relatively safe as it should just entail a job submission
    # Anything more complex may lead to the website hanging...
    while True:
        lines = process.readline()
        if lines == "" and process.poll() is not None:
            break

        if lines:
            print(lines.strip())

    return redirect(url_for("empty"))
