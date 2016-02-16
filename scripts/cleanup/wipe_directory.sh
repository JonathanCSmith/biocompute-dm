#!/usr/bin/env bash
for i in "$@"
do
case ${i} in
    # Username
    -p=*)
    DIRECTORY="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

# Delete the user
find "${DIRECTORY}" -mindepth 1 -delete
exit

