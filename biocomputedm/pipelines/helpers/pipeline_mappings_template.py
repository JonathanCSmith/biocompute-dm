import json

import jsonschema
from biocomputedm import utils
from biocomputedm.pipelines.models import Pipeline, PipelineModule, PipelineModuleOption
from flask import flash

pipeline = \
    '''
    {
    "$schema": "http://json-schema.org/draft-04/schema#",
      "id": "/",
      "type": "object",
      "required": [
        "name",
        "description",
        "pipeline_type",
        "author",
        "version",
        "modules"
      ],
      "properties": {
        "name": {
          "id": "name",
          "type": "string"
        },
        "description": {
          "id": "description",
          "type": "string"
        },
        "pipeline_type": {
          "id": "pipeline_type",
          "type": {
            "enum": [
              "I",
              "II",
              "III"
            ]
          }
        },
        "author": {
          "id": "author",
          "type": "string"
        },
        "version": {
          "id": "version",
          "type": "string"
        },
        "modules": {
          "id": "modules",
          "type": "array",
          "items": {
            "id": "1",
            "type": "object",
            "required": [
              "name",
              "description",
              "executor",
              "options"
            ],
            "properties": {
              "name": {
                "id": "name",
                "type": "string"
              },
              "description": {
                "id": "description",
                "type": "string"
              },
              "executor": {
                "id": "executor",
                "type": "string"
              },
              "options:": {
                "id": "options:",
                "type": "array",
                "items": {
                  "id": "1",
                  "type": "object",
                  "required": [
                    "display_name",
                    "description",
                    "parameter_name",
                    "default_value",
                    "user_interaction_type",
                    "necessary"
                  ],
                  "properties": {
                    "display_name": {
                      "id": "display_name",
                      "type": "string"
                    },
                    "description": {
                      "id": "description",
                      "type": "string"
                    },
                    "parameter_name": {
                      "id": "parameter_name",
                      "type": "string"
                    },
                    "default_value": {
                      "id": "default_value",
                      "type": "string"
                    },
                    "user_interaction_type": {
                      "id": "user_interaction_type",
                      "type": {
                        "enum": [
                          "file",
                          "string",
                          "boolean",
                          "library"
                        ]
                      }
                    },
                    "necessary": {
                      "id": "necessary",
                      "type": "boolean"
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    '''


def validate(file):
    json_instance = json.loads(open(file).read())

    try:
        jsonschema.validate(json_instance, json.loads(pipeline))

    except jsonschema.ValidationError as e:
        flash(e.message + " for file: " + file, "error")
        return False

    except jsonschema.SchemaError as e:
        flash(e.message + " for file: " + file, "error")
        return False

    return True


def build(file):
    json_instance = json.loads(open(file).read())
    name = json_instance.get("name")
    description = json_instance.get("description")
    author = json_instance.get("author")
    version = json_instance.get("version")
    type = json_instance.get("pipeline_type")
    pipeline = Pipeline.query.filter_by(name=name, description=description, author=author, version=version, type=type).first()
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
