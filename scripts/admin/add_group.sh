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

    # Temp directory
    -t=*)
    TEMP_DIRECTORY="${i#*]}"
    shift
    ;;

    *)
    ;;
esac
done

TMP_GROUP_DIRECTORY="${TEMP_DIRECTORY}/${GROUP_NAME}"
TMP_LANDING_DIRECTORY="${TMP_GROUP_DIRECTORY}/staged_files"
GROUP_SFTP_DIRECTORY="${SFTP_ROOT}/${GROUP_NAME}"
#LANDING_DIRECTORY="${GROUP_SFTP_DIRECTORY}/staged_files"

# Make everything in the temporary directory
mkdir "${TMP_GROUP_DIRECTORY}"
chown root:sftpusers "${TMP_GROUP_DIRECTORY}"
chmod 750 "${TMP_GROUP_DIRECTORY}"
mkdir "${TMP_LANDING_DIRECTORY}"
chown "biocompute-DM_group_${GROUP_NAME}":sftpusers "${TMP_LANDING_DIRECTORY}"
chmod 755 "${TMP_LANDING_DIRECTORY}"

# Deploy to actual directory - handles the case where sftp dir is a mount of a remote dir
mv "${TMP_GROUP_DIRECTORY}" "${SFTP_ROOT}"

# Add the user
useradd "biocompute-DM_group_${GROUP_NAME}" -g sftpusers -d "${GROUP_SFTP_DIRECTORY}" -s /sbin/nologin
echo "biocompute-DM_group_${GROUP_NAME}":"${PASSWORD}" | chpasswd

#mkdir "${GROUP_SFTP_DIRECTORY}"
#chown root:sftpusers "${GROUP_SFTP_DIRECTORY}"
#chmod 750 "${GROUP_SFTP_DIRECTORY}"
#
#mkdir "${LANDING_DIRECTORY}"
#chown "biocompute-DM_group_${GROUP_NAME}":sftpusers "${LANDING_DIRECTORY}"
#chmod 775 "${LANDING_DIRECTORY}"

exit

