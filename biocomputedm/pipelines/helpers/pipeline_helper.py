import codecs
import json
import os
import re
import subprocess

import jsonschema
from biocomputedm import utils
from biocomputedm.decorators import async
from biocomputedm.manage.models import Submission, SampleGroup, Sample
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
        "file_regex": {
          "type": "string"
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
        rx = json_instance.get("file_regex")
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
    rx = json_instance.get("file_regex")
    documentation = json_instance.get("documentation_file_name")
    pipeline = Pipeline.query.filter_by(name=name, description=description, author=author, version=version,
                                        type=type, regex=rx, documentation=documentation).first()
    if pipeline is not None:
        pipeline.update(executable=True)
        return False

    pipeline = Pipeline.create(name=name, description=description, author=author, version=version, type=type, regex=rx, documentation=documentation)
    pipeline.update(executable=True)
    count = -1
    for module in json_instance.get("modules"):
        count += 1
        mod = PipelineModule.create(name=module.get("name"),
                                    description=module.get("description"),
                                    executor=module.get("executor"),
                                    execution_index=count,
                                    pipeline=pipeline)

        for option in module.get("options"):
            opt = PipelineModuleOption.create(display_name=option.get("display_name"),
                                              description=option.get("description"),
                                              parameter_name=option.get("parameter_name"),
                                              default_value=option.get("default_value"),
                                              user_interaction_type=option.get("user_interaction_type"),
                                              necessary=option.get("necessary"),
                                              module=mod)

    return True


@async
def execute_module_instance(app, pid="", oid=""):
    try:
        with app.app_context():
            # Check that our objects exist
            pipeline_instance = PipelineInstance.query.filter_by(display_key=pid).first()
            if pipeline_instance is None:
                return

            # Check that we have a module to execute
            from biocomputedm.pipelines.models import get_current_module_instance
            current_module_instance = get_current_module_instance(pipeline_instance)
            if current_module_instance is None:
                return

            # Directories
            local_pipeline_directory = os.path.join(utils.get_path("pipeline_data", "webserver"),
                                                    pipeline_instance.display_key)
            remote_pipeline_directory = os.path.join(utils.get_path("pipeline_data", "hpc"),
                                                     pipeline_instance.display_key)
            local_csv_path = os.path.join(local_pipeline_directory, "data_map.csv")
            csv_path = os.path.join(remote_pipeline_directory, "data_map.csv")

            # Get the object and build it's data path file - only on first module, this stays constant otherwise
            if pipeline_instance.current_execution_index == 0:
                samples_output_directory = os.path.join(remote_pipeline_directory, "samples_output")

                import csv
                if pipeline_instance.pipeline.type == "I":
                    o = Submission.query.filter_by(display_key=oid).first()
                    if o is None:
                        return

                    # Build the csv
                    with open(local_csv_path, "a", newline="") as csvfile:
                        remote_input_path = os.path.join(utils.get_path("submission_data", "hpc"), oid)
                        local_input_path = os.path.join(utils.get_path("submission_data", "webserver"), oid)

                        writer = csv.writer(csvfile)
                        filepaths = next(os.walk(local_input_path))
                        for file in filepaths[1]:
                            try:
                                writer.writerow(
                                        [
                                            oid,
                                            os.path.join(remote_input_path, file),
                                            samples_output_directory,
                                            "EMPTY INFORMATION"
                                        ])
                            except:
                                pass

                        for file in filepaths[2]:
                            try:
                                writer.writerow(
                                        [
                                            oid,
                                            os.path.join(remote_input_path, file),
                                            samples_output_directory,
                                            "EMPTY INFORMATION"
                                        ])
                            except:
                                pass

                else:
                    o = SampleGroup.query.filter_by(display_key=oid).first()
                    if o is None:
                        return

                    # Build the csv
                    with open(local_csv_path, "a", newline="") as csvfile:
                        writer = csv.writer(csvfile)
                        for sample in o.samples.query.all():
                            sample_output_directory = os.path.join(samples_output_directory, sample.display_key)
                            utils.make_directory(sample_output_directory)

                            writer.writerow(
                                    [
                                        sample.display_key,
                                        os.path.join(utils.get_path("sample_data", "hpc"), sample.display_key),
                                        sample_output_directory,
                                        "EMPTY INFORMATION"
                                    ]
                            )

            # Build variables string
            vstring = ""
            for value in current_module_instance.option_values:
                marker = value.option.parameter_name
                result = value.value
                vstring += marker + "=\"" + result + "\","

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

            local_modules_output_directory = os.path.join(local_pipeline_directory, "modules_output")
            local_module_directory = os.path.join(local_modules_output_directory, current_module_instance.module.name)

            remote_modules_output_directory = os.path.join(remote_pipeline_directory, "modules_output")
            remote_module_directory = os.path.join(remote_modules_output_directory, current_module_instance.module.name)

            pipeline_scripts = os.path.join(os.path.join(utils.get_path("scripts", "hpc"), "pipelines"),
                                            pipeline_instance.pipeline.name)
            pipeline_output_directory = os.path.join(local_pipeline_directory, "pipeline_output")
            with open(os.path.join(local_module_directory,
                                   current_module_instance.module.name + "_hpc_submission_out.txt"), "wb") as out, \
                    open(os.path.join(local_module_directory,
                                      current_module_instance.module.name + "_hpc_submission_error.txt"), "wb") as err:
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
                            "-w=" + remote_pipeline_directory,
                            "-mo=" + remote_module_directory,
                            "-po=" + pipeline_output_directory,
                            "-i=" + csv_path,
                            "-v=" + vstring,
                            "-c=" + cleanup_script_path,
                            "-s=" + server
                        ],
                        stdout=out,
                        stderr=err
                )

            pipeline_instance.update(current_execution_status="RUNNING")
            return

    except Exception as e:
        print("There was an exception when executing the current pipeline: " + str(e))
        return


@async
def finish_pipeline_instance(app, pid="", oid=""):
    try:
        with app.app_context():
            # Check that our objects exist
            pipeline_instance = PipelineInstance.query.filter_by(display_key=pid).first()
            if pipeline_instance is None:
                return

            # Build the output directory path
            output_directory = os.path.join(
                    os.path.join(utils.get_path("pipeline_data", "webserver"), pipeline_instance.display_key),
                    "samples_output")

            # Create a sample group
            group = SampleGroup.create(name="Sample Group from Pipeline: " + pipeline_instance.pipeline.name,
                                       creator=pipeline_instance.user,
                                       group=pipeline_instance.user.group,
                                       pipeline=pipeline_instance)

            # Get the submission if relevant
            submission = None
            if pipeline_instance.pipeline.type == "I":
                submission = Submission.query.filter_by(display_key=oid).first()
                if submission is None:
                    print(
                            "Could not locate the submission for this pipeline. The pipeline outcome will not be submitted.")

            # Look in the output directory for folders - these will be our sample names
            filepaths = next(os.walk(output_directory))
            for file in filepaths[1]:
                try:
                    if submission is not None:
                        # Create the sample
                        s = Sample.create(name=file, submission=submission, pipeline=pipeline_instance)

                        # Make the directory
                        utils.make_directory(
                                os.path.join(utils.get_path("sample_data", "webserver"), s.display_key))
                    else:
                        # Find the sample
                        s = Sample.query.filter_by(display_key=file).first()
                        if s is None:
                            continue

                        output_directory = os.path.join(output_directory, s.display_key)

                    # Update the sample group
                    group.samples.append(s)
                    group.save()

                    # Make the data directory and update the sample
                    source = os.path.join(output_directory, file)
                    data_path = os.path.join(os.path.join(utils.get_path("sample_data", "webserver"), s.display_key),
                                             pipeline_instance.display_key)
                    utils.make_directory(data_path)
                    s.pipeline_runs.append(pipeline_instance)
                    s.save()

                    # Transfer the data using an sh
                    script_path = os.path.join(utils.get_path("scripts", "webserver"), "io")
                    script_path = os.path.join(script_path, "move.sh")
                    subprocess.Popen(
                            [
                                "sudo",
                                script_path,
                                "-s=" + source,
                                "-t=" + data_path
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                    ).wait()

                except Exception as e:
                    print("There was an exception when executing the current pipeline: " + str(e))
                    pass

            # Prevent modification for record keeping
            group.update(modifiable=False)

            # Fix up the pipeline status
            data_source = pipeline_instance.data_consigner
            data_source.update(running_pipeline=None)

    except Exception as e:
        print("There was an exception when executing the current pipeline: " + str(e))
        return
