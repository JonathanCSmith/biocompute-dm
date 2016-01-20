#!/usr/bin/env bash
# Note: for SFTP to work in the manner implied below the following conditions must be met
#
# 1) Openssh sftp is allowed and configured correctly using the /etc/ssh/sshd_config file
#   At the start of the file:
#       Sybsystem sftp internal-sftp
#   At the end of the file:
#       Match Group sftpusers
#           ChrootDirectory %h
#           ForceCommand internal-sftp
#           AllowTcpForwarding no
#
# 2) The sftpusers group exists in the system

# Take all script input arguments and build strings for them
for i in "$@"
do
case $i in
    # Ticket for transactions
    -u=*|-user=*)
    USERNAME="biocompute-DM_user_${i#*=}"
    shift
    ;;

    # Job name argument
    -p=*|-password=*)
    PASSWORD="${i#*=}"
    shift
    ;;

    # SFTP Root
    -r=*|-root=*)
    SFTP_ROOT="${i#*=}"
    shift
    ;;

    # Directory Name
    -d=*|-directory=*)
    USER_DIRECTORY_NAME="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

# Add the user, to the sftp group, with password and home directory
SFTP_DIRECTORY=${SFTP_ROOT}"/"${USER_DIRECTORY_NAME}
USER_DIRECTORY=${SFTP_DIRECTORY}"/landing_zone"
ENCRYPTED_PASS=$(mkpasswd -m sha-512 ${PASSWORD})
useradd -g sftpusers -p ${ENCRYPTED_PASS} -d ${USER_DIRECTORY} -s /sbin/nologin ${USERNAME}
exit
