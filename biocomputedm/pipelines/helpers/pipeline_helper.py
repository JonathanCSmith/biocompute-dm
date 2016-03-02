import json
import os
import subprocess

import codecs

import jsonschema
from biocomputedm import utils
from biocomputedm.decorators import async
from biocomputedm.manage.models import Submission, SampleGroup
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
                      "enum": ["boolean", "string", "library", "file", "enum"]
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
        "modules"
      ]
    }
    '''


def validate(file):
    json_instance = json.loads(codecs.open(file, mode="r", encoding="utf-8").read())

    try:
        jsonschema.validate(json_instance, json.loads(pipeline))

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
    pipeline = Pipeline.query.filter_by(name=name, description=description, author=author, version=version,
                                        type=type).first()
    if pipeline is not None:
        pipeline.update(executable=True)
        return False

    pipeline = Pipeline.create(name=name, description=description, author=author, version=version, type=type)
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
def execute_pipeline_instance(app, pid="", oid=""):
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
            remote_pipeline_directory = os.path.join(utils.get_path("pipeline_data", "hpc"), pipeline_instance.display_key)
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
                        writer = csv.writer(csvfile)
                        writer.writerow(
                                [
                                    oid,
                                    os.path.join(utils.get_path("submission_data", "hpc"), oid),
                                    samples_output_directory
                                ])

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
                                        sample_output_directory
                                    ]
                            )

            # Build variables string
            vstring = ""
            for value in current_module_instance.option_values:
                marker = value.option.parameter_name
                if value.option.user_interaction_type == "library":
                    # TODO: input actual library path for webserver
                    result = value.value
                    vstring += marker + "=\"" + result + "\","

                else:
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
            utils.make_directory(local_module_directory)

            remote_modules_output_directory = os.path.join(remote_pipeline_directory, "modules_output")
            remote_module_directory = os.path.join(remote_modules_output_directory, current_module_instance.module.name)

            with open(os.path.join(local_module_directory, current_module_instance.module.name + "_hpc_submission_out.log"), "wb") as out, \
                    open(os.path.join(local_module_directory, current_module_instance.module.name + "_hpc_submission_error.log"), "wb") as err:
                subprocess.Popen(
                        [
                            shell_path,
                            "-u=" + current_app.config["HPC_USERNAME"],
                            "-h=" + current_app.config["NETWORK_PATH_TO_HPC_FROM_WEBSERVER"],
                            "-m=" + current_module_instance.module.name,
                            "-t=" + current_module_instance.display_key,
                            "-e=" + current_module_instance.module.executor,
                            "-l=" + local_module_directory,
                            "-w=" + remote_pipeline_directory,
                            "-o=" + remote_module_directory,
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
        flash("There was an exception when executing the current pipeline: " + str(e), "error")
        return
