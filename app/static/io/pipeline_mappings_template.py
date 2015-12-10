import json

import jsonschema
from app import utils, db
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
                    "option_key",
                    "default",
                    "type"
                  ],
                  "properties": {
                    "option_key": {
                      "id": "option_key",
                      "type": "string"
                    }

                    "default": {
                      "id": "default",
                      "type": "string"
                    }

                    "type": {
                      "id": "type",
                      "enum": [
                        "boolean",
                        "file",
                        "string"
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
    pipeline_instance = utils.get_allowed_pipeline(name, description, author, version)
    if pipeline_instance is not None:
        return False

    pipeline_instance = utils.create_pipeline(name, description, author, version)
    for module in json_instance.get("modules"):
        mod = utils.create_module(module.get("name"), module.get("description"), module.get("executor"), module.get("index_in_execution_order"), pipeline_instance)

        for option in module.get("options"):
            opt = utils.create_option(option.get("option_key"), option.get("default"), mod)

    db.session.add(pipeline_instance)
    db.session.commit()

    return True