import codecs
import json
import jsonschema
from biocomputedm.manage.models import ReferenceData

from flask import flash

reference_data = '''
    {
      "type": "object",
      "properties": {
        "name": {
          "type": "string"
        },
        "description": {
          "type": "string"
        },
        "version": {
          "type": "string"
        }
      },
      "additionalProperties": false,
      "required": [
        "name",
        "description",
        "version"
      ]
    }
'''


def validate(file):
    json_instance = json.loads(codecs.open(file, mode="r", encoding="utf-8").read())

    try:
        jsonschema.validate(json_instance, json.loads(reference_data))

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
    version = json_instance.get("version")

    reference_data_instance = ReferenceData.query.filter_by(name=name, description=description, version=version).first()
    if reference_data_instance is not None:
        reference_data_instance.update(current=True)
        return False

    reference_data_instance = ReferenceData.create(name=name, description=description, version=version)
    reference_data_instance.update(current=True)
    return True
