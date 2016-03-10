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

# Link dir for static file serving
echo "Linking: ${SOURCE} to ${TARGET}"
rm -rf "${TARGET}"
ln -s "${SOURCE}/" "${TARGET}"
exit
