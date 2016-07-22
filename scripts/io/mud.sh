#!/usr/bin/env bash
# Handle args
for i in "$@"
do
case ${i} in
    # Files to move
    -s=*)
    SOURCE="${i#*=}"
    shift
    ;;

    -d=*)
    DIRECTORY="${i#*=}"
    shift
    ;;

    # Submission job number
    -i=*)
    SUBMISSION_ID="${i#*=}"
    shift
    ;;

    # Webserver ip from hpc
    -p=*)
    PORT="${i#*=}"
    shift
    ;;

    # Should unpack
    -u=*)
    UNPACK="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

### NOTE REMOVED UNPACK FOR NOW

IFS=',' read -r -a ARRAY <<< "${SOURCE}"

# Move (and remove original), unpack and delete the newly moved files
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
for FILE in "${ARRAY[@]}"
do
    NEW_FILE="$(basename "${FILE}")"
    mv "${FILE}" "${DIRECTORY}/${NEW_FILE}"

    if [ "$UNPACK" == "True" ]; then
        if (cd "${DIRECTORY}"; "${DIR}"/unpack.sh "${NEW_FILE}"); then
            rm -f "${DIRECTORY}/${NEW_FILE}"
        else
            mv "${DIRECTORY}/${NEW_FILE}" "${FILE}"
            exit
        fi
    fi

    #chown -R biocompute-dm:sftpusers "${NEW_FILE}"
    #chmod -R 660 "${NEW_FILE}"
done

# Post message
curl --form "event=subcomplete" 127.0.0.1:"${PORT}"/message/manage\|"${SUBMISSION_ID}"
exit
