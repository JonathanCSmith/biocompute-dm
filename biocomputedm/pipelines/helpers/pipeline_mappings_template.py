import json
import jsonschema

from biocomputedm import utils
from biocomputedm.extensions import db
from biocomputedm.pipelines.models import create_module
from biocomputedm.pipelines.models import create_pipeline, create_option
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
          "enum": [
            "I",
            "II",
            "III"
          ]
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
              "index_in_execution_order",
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
              "index_in_execution_order": {
                "id": "index_in_execution_order",
                "type": "integer"
              },
              "options:": {
                "id": "options:",
                "type": "array",
                "items": {
                  "id": "1",
                  "type": "object",
                  "required": [
                    "display_name",
                    "parameter_name",
                    "default_value",
                    "user_interaction_type"
                  ],
                  "properties": {
                    "display_name": {
                      "id": "display_name",
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
                      "enum": [
                        "string",
                        "boolean",
                        "library"
                      ]
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
    pipeline = utils.get_allowed_pipeline(name, description, author, version, type)
    if pipeline is not None:
        return False

    pipeline = create_pipeline(name, description, author, version, type)
    for module in json_instance.get("modules"):
        mod = create_module(module.get("name"), module.get("description"), module.get("executor"), module.get("index_in_execution_order"), pipeline)

        for option in module.get("options"):
            opt = create_option(option.get("display_name"), option.get("parameter_name"), option.get("default_value"), option.get("user_interaction_type"), mod)

    return True