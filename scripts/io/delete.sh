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

    # Unknown
    *)
    ;;
esac
done

# We want to move all files within the source directory (doesn't catch . files) into the destination directory
echo "Deleting: ${SOURCE}"
rm -rf "${SOURCE}"
exit
