#!/usr/bin/env bash

for i in "$@"
do
case ${i} in
    # Files to move
    -s=*)
    SOURCE="${i#*=}"
    shift
    ;;

    -d=*)
    TARGET="${i#*=}"
    shift
    ;;

    -f=*)
    FILE="${i*-}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

if [ -z "$FILE" ]; then
    cp -a "${SOURCE}/." "${TARGET}/"

else
    cp -rf "${SOURCE}" "${TARGET}"

fi