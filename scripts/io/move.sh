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

    -t=*)
    TARGET="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

# We want to move all files within the source directory (doesn't catch . files) into the destination directory
mv "${SOURCE}"/. "${TARGET}"
exit
