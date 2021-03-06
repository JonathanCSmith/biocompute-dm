import os
import subprocess

from biocomputedm import utils
from flask import current_app


def create_group_directory(group, realpass):
    path = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "admin"), "add_group.sh")
    process_out = subprocess.Popen(
            [
                "sudo",
                path,
                "-g=" + group.name,
                "-p=" + realpass,
                "-r=" + current_app.config["SFTP_USER_ROOT_PATH"],
                "-t=" + current_app.config["TEMP_DIRECTORY"]
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE
    )

    print(process_out.communicate())
    return


def create_user_directory(user, realpass):
    path = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "admin"), "add_user.sh")
    process_out = subprocess.Popen(
            [
                "sudo",
                path,
                "-u=" + user.username,
                "-p=" + realpass,
                "-r=" + current_app.config["SFTP_USER_ROOT_PATH"],
                "-d=" + user.display_key,
                "-t=" + current_app.config["TEMP_DIRECTORY"]
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE
    )

    print(process_out.communicate())
    return


def change_password(name, realpass):
    path = os.path.join(os.path.join(utils.get_path("scripts", "webserver"), "admin"), "change_password.sh")
    process_out = subprocess.Popen(
        [
            "sudo",
            path,
            name,
            realpass
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE
    )

    print(process_out.communicate())
    return