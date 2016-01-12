import os
import subprocess

from biocomputedm import utils
from flask import current_app


def create_group_directory(group):
    utils.make_directory(group.directory)


def create_user_directory(user, realpass):
    path = os.path.join(current_app.config["SFTP_SCRIPTS_PATH"], "add_user.sh")
    cmd = "sudo " + path
    process_out = subprocess.Popen(
            [
                cmd,
                "-u=" + user.name,
                "-p=" + realpass,
                "-r=" + current_app.config["RAW_DATA_PATH_ON_WEBSERVER"],
                "-d=" + user.display_name
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE
    ).stdout

    # Print out our logs
    while True:
        lines = process_out.readline()
        if lines == "" and process_out.poll() is not None:
            break

        if lines:
            print(lines.strip())

    return
