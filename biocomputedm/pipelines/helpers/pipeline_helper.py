import codecs
import json
import os
import re
import subprocess

import jsonschema
from biocomputedm import utils
from biocomputedm.decorators import async
from biocomputedm.manage.models import Sample, DataGroup, DataItem
from biocomputedm.pipelines.models import Pipeline, PipelineModule, PipelineModuleOption, PipelineInstance
from flask import current_app
from flask import flash

pipeline = \
    '''
    {
      "type": "object",
      "properties": {
        "name": {
          "type": "string"
        },
        "description": {
          "type": "string"
        },
        "pipeline_type": {
          "enum": ["I", "II", "III"]
        },
        "author": {
          "type": "string"
        },
        "version": {
          "type": "string"
        },
        "regex_type": {
          "enum": ["AND", "OR"]
        },
        "file_regex": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "minItems": 1,
          "uniqueItems": true
        },
        "documentation_file_name": {
          "type": "string"
        },
        "modules": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "name": {
                "type": "string"
              },
              "description": {
                "type": "string"
              },
              "executor": {
                "type": "string"
              },
              "options": {
                "type": "array",
                "items": {
                  "type": "object",
                  "properties": {
                    "display_name": {
                      "type": "string"
                    },
                    "parameter_name": {
                      "type": "string"
                    },
                    "default_value": {
                      "type": "string"
                    },
                    "user_interaction_type": {
                      "enum": ["boolean", "string", "reference", "file", "enum"]
                    },
                    "necessary": {
                      "type": "boolean"
                    },
                    "description": {
                      "type": "string"
                    }
                  },
                  "additionalProperties": false,
                  "required": [
                    "display_name",
                    "parameter_name",
                    "default_value",
                    "user_interaction_type",
                    "necessary",
                    "description"
                  ]
                }
              }
            },
            "additionalProperties": false,
            "required": [
              "name",
              "description",
              "executor",
              "options"
            ]
          }
        }
      },
      "additionalProperties": false,
      "required": [
        "name",
        "description",
        "pipeline_type",
        "author",
        "version",
        "regex_type",
        "file_regex",
        "documentation_file_name",
        "modules"
      ]
    }
    '''


# Validate whether a file provided can be used as a valid pipeline
def validate(file):
    json_instance = json.loads(codecs.open(file, mode="r", encoding="utf-8").read())

    try:
        jsonschema.validate(json_instance, json.loads(pipeline))

        # Validate documentation presence
        documentation = json_instance.get("documentation_file_name")
        if documentation is None or documentation == "":
            flash("File: " + file + " did not provide a valid documentation link", "error")
            return False

        pipeline_directory = os.path.dirname(file)
        if not os.path.isfile(os.path.join(pipeline_directory, documentation)):
            flash("File: " + file + " did not provide a valid documentation link", "error")
            return False

        # Validate regex presence
        rxs = json_instance.get("file_regex")
        for rx in rxs:
            if rx is None or rx == "":
                flash("File: " + file + " did not provide a valid file regex", "error")
                return False

            try:
                re.compile(rx)
            except re.error:
                flash("File: " + file + " did not provide a valid file regex", "error")
                return False

    except jsonschema.ValidationError as e:
        flash(e.message + " for file: " + file, "error")
        return False

    except jsonschema.SchemaError as e:
        flash(e.message + " for file: " + file, "error")
        return False

    except Exception as e:
        flash("File: " + file + " received: " + str(e), "error")
        return False

    return True


def build(file):
    json_instance = json.loads(codecs.open(file, mode="r", encoding="utf-8").read())
    name = json_instance.get("name")
    description = json_instance.get("description")
    author = json_instance.get("author")
    version = json_instance.get("version")
    type = json_instance.get("pipeline_type")
    rx_type = json_instance.get("regex_type")
    rxs = json_instance.get("file_regex")
    documentation = json_instance.get("documentation_file_name")

    regex = ""
    for rx in rxs:
        regex += rx
        regex += "###"
    regex = regex[:-3]  # Remove that pesky additional marker!

    pipeline = Pipeline.query.filter_by(
        name=name,
        description=description,
        author=author,
        version=version
    ).first()

    if pipeline is not None:
        pipeline.update(executable=True)
        return False

    pipeline = Pipeline.create(
        name=name,
        description=description,
        author=author,
        version=version,
        type=type,
        regex=regex,
        regex_type=rx_type,
        documentation=documentation
    )

    pipeline.update(executable=True)
    count = -1
    for module in json_instance.get("modules"):
        count += 1
        mod = PipelineModule.create(
            name=module.get("name"),
            description=module.get("description"),
            executor=module.get("executor"),
            execution_index=count,
            pipeline=pipeline
        )

        for option in module.get("options"):
            opt = PipelineModuleOption.create(
                display_name=option.get("display_name"),
                description=option.get("description"),
                parameter_name=option.get("parameter_name"),
                default_value=option.get("default_value"),
                user_interaction_type=option.get("user_interaction_type"),
                necessary=option.get("necessary"),
                module=mod
            )

    return True


def initialise_running_pipeline(running_pipeline_id, source_data_group_id):
    target_pipeline = None
    try:

        # Pipeline properties
        if running_pipeline_id == "":
            flash("Could not identify the running pipeline in order to initialise.", "error")
            return

        target_pipeline = PipelineInstance.query.filter_by(display_key=running_pipeline_id).first()
        if target_pipeline is None:
            flash("Could not identify the running pipeline in order to initialise.", "error")
            return

        if target_pipeline.current_execution_status != "NOT_STARTED":
            flash("Could not initialise the pipeline as it is not in the 'NOT_STARTED' phase.", "error")
            return

        # Source data group properties
        if source_data_group_id == "":
            flash("Could not identify the source data for initialising the pipeline.", "error")
            target_pipeline.update(current_execution_status="ERROR")
            return

        source_data_group = DataGroup.query.filter_by(display_key=source_data_group_id).first()
        if source_data_group is None:
            flash("Could not identify the source data for initialising the pipeline.", "error")
            target_pipeline.update(current_execution_status="ERROR")
            return

        # Build the pipeline's data
        local_pipeline_path = os.path.join(utils.get_path("pipeline_data", "webserver"), target_pipeline.display_key)
        remote_pipeline_path = os.path.join(utils.get_path("pipeline_data", "hpc"), target_pipeline.display_key)
        remote_output_path = os.path.join(remote_pipeline_path, "samples_output")
        local_csv_path = os.path.join(local_pipeline_path, "data_map.csv")
        import csv
        with open(local_csv_path, "a", newline="") as csv_file:
            writer = csv.writer(csv_file)
            for source_data in source_data_group.data:
                try:
                    # Special case I pipelines as they don't have samples yet
                    if target_pipeline.pipeline.type == "I":
                        key = source_data.display_key
                        input = os.path.join(os.path.join(utils.get_path("submission_data", "hpc"), source_data.unlocalised_path), source_data.name)
                        output = remote_output_path

                    else:
                        key = source_data.sample.display_key
                        input = os.path.join(utils.get_path("sample_data", "hpc"), source_data.unlocalised_path)
                        output = os.path.join(remote_output_path, source_data.unlocalised_path)

                    # Generate a csv row where columns 1 = sample, 2 = data input path, 3 = data output path
                    writer.writerow(
                        [
                            key,
                            input,
                            output,
                            "EMPTY INFORMATION"
                        ]
                    )

                except:
                    flash("Could not write row for: " + source_data.display_key, "error")
                    pass

        target_pipeline.update(current_execution_status="WAITING", current_execution_index=0)
        source_data_group.update(running=True)
        source_data_group.pipeline_instances.append(target_pipeline)
        source_data_group.save()

    except Exception as e:
        flash("There was an exception when executing the current pipeline: " + str(e), "error")
        if target_pipeline is not None:
            target_pipeline.update(current_execution_status="ERROR")

        return


@async
def execute_pipeline_module(app, running_pipeline_id=""):
    running_pipeline = None
    try:
        with app.app_context():

            # Pipeline properties
            if running_pipeline_id == "":
                app.logger.error("Could not identify the running pipeline in order to initialise.")
                return

            running_pipeline = PipelineInstance.query.filter_by(display_key=running_pipeline_id).first()
            if running_pipeline is None:
                app.logger.error("Could not identify the running pipeline in order to initialise.")
                return

            if running_pipeline.current_execution_status != "WAITING":
                app.logger.error("Could not execute the next pipeline module as it is not in the 'WAITING' phase.")
                return

            # Module properties
            from biocomputedm.pipelines.models import get_current_module_instance
            current_module_instance = get_current_module_instance(running_pipeline)
            if current_module_instance is None:
                app.logger.error("Could not identify the current module instance.")
                running_pipeline.update(current_execution_status="ERROR")
                return

            # Build the pipeline's data
            local_pipeline_path = os.path.join(utils.get_path("pipeline_data", "webserver"), running_pipeline.display_key)
            remote_pipeline_path = os.path.join(utils.get_path("pipeline_data", "hpc"), running_pipeline.display_key)

            # Build variables string
            vstring = ""
            for value in current_module_instance.option_values:
                marker = value.option.parameter_name
                result = value.value
                vstring += marker + "=\"" + result.replace(",", "%%___%%") + "\","

            if len(vstring) != 0:
                vstring = vstring[:-1]

            # Submit module w/ options into HPC - note the cwd
            shell_path = os.path.join(utils.get_path("scripts", "webserver"), "jobs")
            cleanup_script_path = os.path.join(os.path.join(utils.get_path("scripts", "hpc"), "jobs"), "cleanup.sh")
            if current_app.config["HPC_DEBUG"] == "False":
                shell_path = os.path.join(shell_path, "submit_job.sh")
                server = current_app.config["NETWORK_PATH_TO_WEBSERVER_FROM_HPC"]
            else:
                shell_path = os.path.join(shell_path, "fake_submit_job.sh")
                server = current_app.config["LOCAL_WEBSERVER_PORT"]

            local_modules_output_directory = os.path.join(local_pipeline_path, "modules_output")
            local_module_directory = os.path.join(local_modules_output_directory, current_module_instance.module.name)

            remote_modules_output_directory = os.path.join(remote_pipeline_path, "modules_output")
            remote_module_directory = os.path.join(remote_modules_output_directory, current_module_instance.module.name)

            pipeline_scripts = os.path.join(os.path.join(utils.get_path("scripts", "hpc"), "pipelines"), running_pipeline.pipeline.name)
            pipeline_output_directory = os.path.join(local_pipeline_path, "pipeline_output")
            with open(os.path.join(local_module_directory, current_module_instance.module.name + "_hpc_submission_out.txt"), "wb") as out, \
                    open(os.path.join(local_module_directory, current_module_instance.module.name + "_hpc_submission_error.txt"), "wb") as err:
                subprocess.Popen(
                    [
                        shell_path,
                        "-u=" + current_app.config["HPC_USERNAME"],
                        "-h=" + current_app.config["NETWORK_PATH_TO_HPC_FROM_WEBSERVER"],
                        "-m=" + current_module_instance.module.name,
                        "-t=" + current_module_instance.display_key,
                        "-e=" + current_module_instance.module.executor,
                        "-l=" + local_module_directory,
                        "-p=" + pipeline_scripts,
                        "-w=" + remote_pipeline_path,
                        "-mo=" + remote_module_directory,
                        "-po=" + pipeline_output_directory,
                        "-i=" + os.path.join(remote_pipeline_path, "data_map.csv"),
                        "-v=" + vstring,
                        "-c=" + cleanup_script_path,
                        "-s=" + server
                    ],
                    stdout=out,
                    stderr=err
                )

            running_pipeline.update(current_execution_status="RUNNING")

    except Exception as e:
        app.logger.error("There was an exception when executing the current pipeline: " + str(e))
        if running_pipeline is not None:
            running_pipeline.update(current_execution_status="ERROR")

        return


@async
def finish_pipeline_instance(app, running_pipeline_id):
    try:
        with app.app_context():

            # Pipeline properties
            if running_pipeline_id == "":
                app.logger.error("Could not identify the running pipeline in order to initialise.")
                return

            running_pipeline = PipelineInstance.query.filter_by(display_key=running_pipeline_id).first()
            if running_pipeline is None:
                app.logger.error("Could not identify the running pipeline in order to initialise.")
                return

            if running_pipeline.current_execution_status != "STOPPED":
                app.logger.error("Could not execute the next pipeline module as it is not in the 'WAITING' phase.")
                return

            if running_pipeline.current_execution_index != len(running_pipeline.pipeline.modules.all()) - 1:
                app.logger.error("Cannot finalise pipeline when it has not finished all of its modules.")
                return

            # Dummy container
            data_group = None

            # Walk the output directory to find samples
            local_pipeline_directory = os.path.join(utils.get_path("pipeline_data", "webserver"), running_pipeline.display_key);
            local_pipeline_sample_output_directory = os.path.join(local_pipeline_directory, "samples_output")
            filepaths = next(os.walk(local_pipeline_sample_output_directory))
            for source in filepaths[1]:

                # Verify that there is data within  the directory, if so index it
                sample_filepaths = next(os.walk(os.path.join(local_pipeline_sample_output_directory, source)))
                if sample_filepaths[1] or sample_filepaths[2]:
                    # Create a new data group for the outputs
                    if data_group is None:
                        data_group = DataGroup.create(
                            name="Data Group from the " + running_pipeline.pipeline.name + " pipeline",
                            user=running_pipeline.user,
                            group=running_pipeline.group,
                            source_pipeline=running_pipeline
                        )

                    # Create new samples
                    if running_pipeline.pipeline.type == "I":
                        s = Sample.create(name=source, pipeline=running_pipeline)
                        utils.make_directory(os.path.join(utils.get_path("sample_data", "webserver"), s.display_key))

                    # Link new data to old samples
                    else:
                        s = Sample.query.filter_by(display_key=source).first()
                        if s is None:
                            s = Sample.create(name=source, pipeline=running_pipeline)
                            utils.make_directory(os.path.join(utils.get_path("sample_data", "webserver"), s.display_key))

                    # Transfer the data using an sh
                    destination_path = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), s.display_key), running_pipeline.display_key)
                    utils.make_directory(destination_path)
                    script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
                    script_path = os.path.join(script_path, "move.sh")
                    subprocess.Popen(
                        [
                            "sudo",
                            script_path,
                            "-s=" + os.path.join(local_pipeline_sample_output_directory, source),
                            "-t=" + destination_path
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    ).wait()

                    output_files = next(os.walk(destination_path))
                    for file in output_files[1]:
                        # Create a new data item and append to sample and data group
                        item = DataItem.create(
                            name=file,
                            unlocalised_path=os.path.join(s.display_key, running_pipeline.display_key),
                            data_group=data_group,
                            group=running_pipeline.group
                        )

                        s.data.append(item)

                    for file in output_files[2]:
                        # Create a new data item and append to sample and data group
                        item = DataItem.create(
                            name=file,
                            unlocalised_path=os.path.join(s.display_key, running_pipeline.display_key),
                            data_group=data_group,
                            group=running_pipeline.group
                        )

                        s.data.append(item)

            # Build pipeline outputs (no sample association)
            local_pipeline_pipeline_output_directory = os.path.join(local_pipeline_directory, "pipeline_output")
            filepaths = next(os.walk(local_pipeline_pipeline_output_directory))
            if filepaths[1] or filepaths[2]:
                # Create a new data group for the outputs
                pipeline_data_group = DataGroup.create(
                    name="Runtime Data from the " + running_pipeline.pipeline.name + " pipeline",
                    user=running_pipeline.user,
                    group=running_pipeline.group,
                    source_pipeline=None
                )

                for file in filepaths[1]:
                    item = DataItem.create(
                        name=file,
                        unlocalised_path=os.path.join(running_pipeline.display_key, "pipeline_output"),
                        data_group=pipeline_data_group,
                        group=running_pipeline.group
                    )

                for file in filepaths[2]:
                    item = DataItem.create(
                        name=file,
                        unlocalised_path=os.path.join(running_pipeline.display_key, "pipeline_output"),
                        data_group=pipeline_data_group,
                        group=running_pipeline.group
                    )

                running_pipeline.update(pipeline_output=pipeline_data_group)

            # Build module outputs (no sample association)
            local_pipeline_module_output_directory = os.path.join(local_pipeline_directory, "modules_output")
            for module_instance in running_pipeline.module_instances:

                local_module_output_directory = os.path.join(local_pipeline_module_output_directory, module_instance.module.name)
                filepaths = next(os.walk(local_module_output_directory))
                if filepaths[1] or filepaths[2]:
                    # Create a new data group for the outputs
                    module_data_group = DataGroup.create(
                        name="Module Data from the " + running_pipeline.pipeline.name + " pipeline, module: " + module_instance.module.name,
                        user=running_pipeline.user,
                        group=running_pipeline.group,
                        source_pipeline=None
                    )

                    for file in filepaths[1]:
                        item = DataItem.create(
                            name=file,
                            unlocalised_path=os.path.join(os.path.join(running_pipeline.display_key, "modules_output"), module_instance.module.name),
                            data_group=module_data_group,
                            group=running_pipeline.group
                        )

                    for file in filepaths[2]:
                        item = DataItem.create(
                            name=file,
                            unlocalised_path=os.path.join(os.path.join(running_pipeline.display_key, "modules_output"), module_instance.module.name),
                            data_group=module_data_group,
                            group=running_pipeline.group
                        )

                    module_instance.update(module_output=module_data_group)

            running_pipeline.update(current_execution_status="FINISHED")
            source_data_group = DataGroup.query.filter_by(display_key=running_pipeline.data_consignor.display_key).first()
            source_data_group.update(running=True)

    except Exception as e:
        app.logger.error("There was an exception when executing the current pipeline: " + str(e))
        return
