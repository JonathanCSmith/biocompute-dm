import os

import datetime

from flask import current_app
from flask import flash

__author__ = 'jon'


# Function to retrieve a path variable based on the target and the environment
def get_path(target, environment):
    if target == "scripts":
        target_path = current_app.config["MANAGEMENT_SCRIPTS_PATH_AFTER_RELATIVE_ROOT"]
    elif target == "pipeline_data":
        target_path = current_app.config["PIPELINE_DATA_PATH_AFTER_RELATIVE_ROOT"]
    elif target == "submission_data":
        target_path = current_app.config["SUBMISSION_DATA_PATH_AFTER_RELATIVE_ROOT"]
    elif target == "sample_data":
        target_path = current_app.config["SAMPLE_DATA_PATH_AFTER_RELATIVE_ROOT"]
    elif target == "project_data":
        target_path = current_app.config["PROJECT_DATA_PATH_AFTER_RELATIVE_ROOT"]
    elif target == "reference_data":
        target_path = current_app.config["REFERENCE_DATA_PATH_AFTER_RELATIVE_ROOT"]
    else:
        raise ValueError(
                "Valid inputs for target are: scripts, pipeline_data, submission_data, sample_data, project_data and reference data")

    if environment == "webserver":
        environment_path = current_app.config["WEBSERVER_ROOT_PATH"]
    elif environment == "hpc":
        environment_path = current_app.config["HPC_ROOT_PATH"]
    elif environment == "serve":
        return os.path.join("serve", target_path.split("/")[-1])
    else:
        raise ValueError("environment must either be: serve, webserver or hpc")

    return os.path.join(environment_path, target_path)


# Function to display errors
def flash_errors(form):
    for field, errors in form.errors.items():
        for error in errors:
            flash(u"Error in the %s field - %s" % (
                getattr(form, field).label.text,
                error
            ), "error")


def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def get_current_date():
    return datetime.datetime.now()
