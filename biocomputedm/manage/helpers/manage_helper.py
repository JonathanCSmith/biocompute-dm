import os
import subprocess

from biocomputedm import utils
from biocomputedm.admin.models import Person
from biocomputedm.decorators import async


@async
def copy_data_to_staging(app, oid, type, user_key):
    try:
        with app.app_context():
            # Submit the directory and uploaded file information to our script
            script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
            script_path = os.path.join(script_path, "copy.sh")

            person = Person.query.filter_by(display_key=user_key).first()
            if person is None:
                app.logger.error("Could not identify the provided person: " + user_key)
                return

            if type == "pipeline_output":
                pipeline_instance = person.group.pipeline_instances.filter_by(display_key=oid).first()
                if pipeline_instance is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Pipeline_Data_From_" + pipeline_instance.pipeline.name)
                utils.make_directory(output_directory_path)

                source_directory = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key), "pipeline_output")

                # Execute our copy script
                subprocess.Popen(
                    [
                        "sudo",
                        script_path,
                        "-s=" + source_directory,
                        "-d=" + output_directory_path,
                        "-f=folder"
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )  # We are allowing this to execute on it's own - no need to monitor

            elif type == "pipeline_sample_group":
                pipeline_instance = person.group.pipeline_instances.filter_by(display_key=oid).first()
                if pipeline_instance is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                for data in pipeline_instance.sample_output.data:
                    source_directory = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), data.unlocalised_path), data.name)

                    # Create the directory to hold the data
                    output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), data.sample.name)
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
                project = person.group.projects.query.filter_by(display_key=oid).first()
                if project is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                output_directory_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Project_Pipeline_Data_From_" + project.name)
                utils.make_directory(output_directory_path)

                for outputs in project.pipeline_outputs:
                    for data in outputs.data:
                        source_directory = os.path.join(os.path.join(utils.get_path("pipeline_data", "webserver"), data.unlocalised_path), data.name)

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
                project = person.group.projects.query.filter_by(display_key=oid).first()
                if project is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                project_output_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Project_Sample_Data_From_" + project.name)
                utils.make_directory(project_output_path)

                for sample in project.samples:
                    # Create the directory to hold the data
                    sample_output_path = os.path.join(project_output_path, "Sample_Data_For_" + sample.name)
                    utils.make_directory(sample_output_path)

                    for data in sample.data:
                        source_directory = os.path.join(os.path.join(utils.get_path("samples_data", "webserver"), data.unlocalised_path), data.name)

                        # Create the directory to hold the data
                        output_directory_path = os.path.join(sample_output_path, data.sample.name)
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

            elif type == "sample":
                sample = person.group.samples.filter_by(display_key=oid).first()
                if sample is None:
                    app.logger.error("Could not identify the provided object: " + oid + " with type: " + type)
                    return

                # Create the directory to hold the data
                sample_output_path = os.path.join(os.path.join(os.path.join(app.config["SFTP_USER_ROOT_PATH"], user_key), "staged_files"), "Sample_Data_For_" + sample.name)
                utils.make_directory(sample_output_path)

                for data in sample.data:
                    source_directory = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), data.unlocalised_path), data.name)
                    output_directory_path = os.path.join(sample_output_path, data.name)

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
