import os

from biocomputedm import utils
from biocomputedm.decorators import login_required
from biocomputedm.pipelines.models import Pipeline
from flask import Blueprint, redirect, url_for, current_app, render_template

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
    utils.refresh_pipelines()
    return redirect(url_for("admin.administrate"))


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
                                 app.config["SUBMISSIONS_PATH_ON_HPC"] + "/" + sg.source_pipeline.display_key + "/" + sample.display_key,
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
                                     app.config["SAMPLES_PATH_ON_HPC"] + "/" + sg.source_pipeline.display_key + "/" + sample.display_key,
                                     app.config["INVESTIGATIONS_PATH_ON_HPC"] + "/" + pi.display_key + "/" + sample.display_key])

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
