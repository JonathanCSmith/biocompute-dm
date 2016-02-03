import os
import subprocess

import getpass

from biocomputedm import utils
from flask import current_app


def create_group_directory(group):
    # TODO: Currently group members don't have write access
    utils.make_directory(os.path.join(current_app.config["SFTP_USER_ROOT_PATH"], group.name))


def create_user_directory(user, realpass):
    # TODO: Currently users are not members of their group!
    path = os.path.join(current_app.config["SFTP_SCRIPTS_PATH"], "add_user.sh")
    process_out = subprocess.Popen(
            [
                "sudo",
                path,
                "-u=" + user.username,
                "-p=" + realpass,
                "-r=" + current_app.config["SFTP_USER_ROOT_PATH"],
                "-d=" + user.display_key
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE
    )

    print(process_out.communicate())
    return
