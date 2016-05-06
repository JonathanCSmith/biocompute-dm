#!/usr/bin/env bash

for i in "$@"
do
case $i in
    # Group
    -g=*|-group=*)
    GROUP_NAME="${i#*=}"
    shift
    ;;

    # Job name argument
    -p=*|-password=*)
    PASSWORD="${i#*=}"
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
LANDING_DIRECTORY="${GROUP_SFTP_DIRECTORY}/staged_files"

# Add the user
useradd "biocompute-DM_group_${GROUP_NAME}" -g sftpusers -d "${GROUP_SFTP_DIRECTORY}" -s /sbin/nologin
echo "biocompute-DM_group_${GROUP_NAME}":"${PASSWORD}" | chpasswd

mkdir "${GROUP_SFTP_DIRECTORY}"
chown root:sftpusers "${GROUP_SFTP_DIRECTORY}"
chmod 750 "${GROUP_SFTP_DIRECTORY}"

mkdir "${LANDING_DIRECTORY}"
chown "biocompute-DM_group_${GROUP_NAME}":sftpusers "${LANDING_DIRECTORY}"
chmod 775 "${LANDING_DIRECTORY}"

exit

