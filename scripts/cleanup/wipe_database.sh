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
    # Username
    -u=*)
    USERNAME="${i#*=}"
    shift
    ;;

    # Job name argument
    -p=*)
    PASSWORD="${i#*=}"
    shift
    ;;

    # SFTP Root
    -l=*)
    LOCATION="${i#*=}"
    shift
    ;;

    # Directory Name
    -n=*)
    NAME="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

QUERY="DROP DATABASE ${NAME};CREATE DATABASE ${NAME};\q"
mysql -u "${USERNAME}" --password="${PASSWORD}" -h "${LOCATION}" -e "${QUERY}" -D "${NAME}"

exit
