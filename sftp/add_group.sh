#!/usr/bin/env bash

for i in "$@"
do
case $i in
    # Group
    -g=*|-group=*)
    GROUP_NAME="${i#*=}"
    shift
    ;;

    -r=*|-root=*)
    SFTP_ROOT="${i#*=}"
    shift
    ;;

    *)
    ;;
esac
done

GROUP_SFTP_DIRECTORY="${SFTP_ROOT}/${GROUP_NAME}"
LANDING_DIRECTORY="${GROUP_SFTP_DIRECTORY}/landing_zone"

mkdir "${GROUP_SFTP_DIRECTORY}"
chown root:sftpusers "${GROUP_SFTP_DIRECTORY}"
chmod 750 "${GROUP_SFTP_DIRECTORY}"

mkdir "${LANDING_DIRECTORY}"
chown "${root:sftpusers}" "${LANDING_DIRECTORY}"
chmod 755 "${LANDING_DIRECTORY}"

exit

