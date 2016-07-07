import os
import subprocess

import re

from biocomputedm import utils
from biocomputedm.admin.models import Person
from biocomputedm.decorators import async
from biocomputedm.manage.models import Submission, DataGroup, Project, Sample
from biocomputedm.pipelines.models import Pipeline, PipelineInstance


@async
def copy_data_to_staging(app, oid, type, user_key):
    try:
        with app.app_context():
            person = Person.query.filter_by(display_key=user_key).first()
            if person is None:
                app.logger.error("Could not identify the provided person: " + user_key)
                return

            # Submit the directory and uploaded file information to our script
            script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
            script_path = os.path.join(script_path, "copy_folder.sh")

            if type == "pipeline_output":
                pipeline_instance = person.group.pipeline_instances.filter_by(display_key=oid).first()
                if pipeline_instance is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                source_directory = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key), "pipeline_output")
                output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Pipeline_Data_From_" + pipeline_instance.pipeline.name)
                utils.make_directory(output_directory_path)

                # Execute our copy script
                subprocess.Popen(
                    [
                        "sudo",
                        script_path,
                        "-s=" + source_directory,
                        "-d=" + output_directory_path
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )  # We are allowing this to execute on it's own - no need to monitor

            elif type == "project_pipeline_output":
                project = Project.query.filter_by(display_key=oid).first()
                if project is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                parent_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Project_Pipeline_Data_From_" + project.name)
                utils.make_directory(parent_directory_path)

                for pipeline_output in project.pipeline_outputs:
                    pipeline_key = pipeline_output.data[0].unlocalised_path.split("/")[0]
                    pipeline = Pipeline.query.filter_by(display_key=pipeline_key).first()
                    source_directory = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_key, "pipeline_output"))
                    output_directory_path = os.path.join(parent_directory_path, pipeline.name)
                    utils.make_directory(output_directory_path)

                    # Execute our copy script
                    subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-s=" + source_directory,
                            "-d=" + output_directory_path
                        ],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )  # We are allowing this to execute on it's own - no need to monitor

            elif type == "pipeline_sample_group":
                pipeline_instance = PipelineInstance.filter_by(display_key=oid).first()
                if pipeline_instance is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                done = []
                for data in pipeline_instance.sample_output.data:
                    if data.unlocalised_path in done:
                        continue

                    done.append(data.unlocalised_path)
                    source_directory = os.path.join(utils.get_path("sample_data", "webserver"), data.unlocalised_path)
                    output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Sample_Data_For_Sample_" + data.sample.name)
                    utils.make_directory(output_directory_path)

                    # Execute our copy script
                    subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-s=" + source_directory,
                            "-d=" + output_directory_path
                        ],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )  # We are allowing this to execute on it's own - no need to monitor

            elif type == "project_sample_group":
                project = Project.query.filter_by(display_key=oid).first()
                if project is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                project_output_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Project_Sample_Data_From_" + project.name)
                utils.make_directory(project_output_path)

                for sample in project.samples:
                    source_directory = os.path.join(utils.get_path("sample_data", "webserver"), sample.display_key)
                    sample_output_path = os.path.join(project_output_path, "Sample_Data_For_" + sample.name)
                    utils.make_directory(sample_output_path)

                    # Execute our copy script
                    subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-s=" + source_directory,
                            "-d=" + sample_output_path
                        ],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )  # We are allowing this to execute on it's own - no need to monitor

            elif type == "sample":
                sample = Sample.filter_by(display_key=oid).first()
                if sample is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Sample_Data_For_Sample_" + sample.name)
                utils.make_directory(output_directory_path)
                source_directory = os.path.join(utils.get_path("sample_data", "webserver"), sample.display_key)

                # Execute our copy script
                subprocess.Popen(
                    [
                        "sudo",
                        script_path,
                        "-s=" + source_directory,
                        "-d=" + output_directory_path
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )  # We are allowing this to execute on it's own - no need to monitor

    except Exception as ex:
        app.logger.error("There was an exception when copying data from oid: " + oid + " type: " + type + " with error: " + str(ex))
        return


def calculate_viable_pipelines_for_submission(oid):
    # Get the submission
    submission = Submission.query.filter_by(display_key=oid).first()
    if submission is None:
        return
    data_group = submission.data_group

    # list of the available type I pipelines
    pipelines = Pipeline.query.filter_by(type="I", executable=True)

    # use the regex information to determine if a pipeline can run
    valid_pipelines = ""
    for pipeline in pipelines:
        rxs = pipeline.regex.split("###")

        if pipeline.regex_type == "OR":
            valid_pipeline = False
            for rx in rxs:
                found = False
                rx_matcher = re.compile(rx)
                for dirpath, dirnames, files in os.walk(os.path.join(utils.get_path("submission_data", "webserver"), oid)):
                    for file in files:
                        if rx_matcher.match(file):
                            found = True
                            break

                    if found:
                        break

                if found:
                    valid_pipeline = True
                    break

            if valid_pipeline:
                if valid_pipelines != "":
                    valid_pipelines += ","

                valid_pipelines += pipeline.display_key
                continue

        elif pipeline.regex_type == "AND":
            valid_pipeline = True
            for rx in rxs:
                found = False
                rx_matcher = re.compile(rx)
                for dirpath, dirnames, files in os.walk(os.path.join(utils.get_path("submission_data", "webserver"), oid)):
                    for file in files:
                        if rx_matcher.match(file):
                            found = True
                            break

                    if found:
                        break

                if not found:
                    valid_pipeline = False
                    break

            if valid_pipeline:
                if valid_pipelines != "":
                    valid_pipelines += ","

                valid_pipelines += pipeline.display_key
                continue

    data_group.update(valid_pipelines=valid_pipelines)


def calculate_viable_pipelines_for_data_group(oid):
    data_group = DataGroup.query.filter_by(display_key=oid).first()

    pipelines = Pipeline.query.filter((Pipeline.type == "II") | (Pipeline.type == "III")).filter_by(executable=True)

    # use the regex information to determine if a pipeline can run
    valid_pipelines = ""
    root_path = utils.get_path("sample_data", "webserver")
    files = []
    for data in data_group.data:
        data_path = os.path.join(os.path.join(root_path, data.unlocalised_path), data.name)
        if not os.path.isfile(data_path):
            for dirpath, dirnames, additional_files in os.walk(data_path):
                files.extend(additional_files)

        else:
            files.append(data.name)

    for pipeline in pipelines:
        rxs = pipeline.regex.split("###")

        if pipeline.regex_type == "OR":
            valid_pipeline = False
            for rx in rxs:
                found = False
                rx_matcher = re.compile(rx)

                for file in files:
                    if rx_matcher.match(file):
                        found = True
                        break

                if found:
                    valid_pipeline = True
                    break

            if valid_pipeline:
                if valid_pipelines != "":
                    valid_pipelines +=","

                valid_pipelines += pipeline.display_key
                continue

        elif pipeline.regex_type == "AND":
            valid_pipeline = True
            for rx in rxs:
                found = False
                rx_matcher = re.compile(rx)
                for file in files:
                    if rx_matcher.match(file):
                        found = True
                        break

                if not found:
                    valid_pipeline = False
                    break

            if valid_pipeline:
                if valid_pipelines != "":
                    valid_pipelines += ","

                valid_pipelines += pipeline.display_key
                continue

    data_group.update(valid_pipelines=valid_pipelines)


@async
def asynchronously_calculate_viable_pipelines_for_submission(app, oid):
    try:
        with app.app_context():
            calculate_viable_pipelines_for_submission(oid)

    except Exception as e:
        app.logger.error("There was an exception when calculating viable pipelines: " + str(e))
        return


@async
def asynchronously_calculate_viable_pipelines_for_data_group(app, oid):
    try:
        with app.app_context():
            calculate_viable_pipelines_for_data_group(oid)

    except Exception as e:
        app.logger.error("There was an exception when calculating viable pipelines: " + str(e))
        return

