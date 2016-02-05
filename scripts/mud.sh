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

    # Unknown
    *)
    ;;
esac
done

IFS=',' read -r -a ARRAY <<< "${SOURCE}"

# Move (and remove original), unpack and delete the newly moved files
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
for FILE in "${ARRAY[@]}"
do
    mv "${FILE}" "${DIRECTORY}"
    NEW_FILE="${DIRECTORY}/$(basename "${FILE}")"
    (cd "${DIRECTORY}"; "${DIR}"/unpack.sh "${NEW_FILE}")
    rm -f "${NEW_FILE}"
done

# Post message
curl --form "submission=${SUBMISSION_ID}" 127.0.0.1:5000/manage_message
exit
